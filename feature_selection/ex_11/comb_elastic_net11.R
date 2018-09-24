setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_11/")

library(tictoc)
library(limma)
library(glmnet)

set.seed(12)

tic()

# Combined protein and gene probe set
df.gene.all <- read.csv("../../data/ex_11/gene_train_body.csv", header=TRUE, row.names = 1)  #Gene test/train set
df.prot.all <- read.csv("../../data/ex_11/prot_train_body.csv", header=TRUE, row.names = 1)  # Protein train set
df.gp.all <- cbind(df.prot.all, df.gene.all)
df.meta <- read.csv("../../data/ex_11/gp_train_meta.csv", header=TRUE, row.names = 1)
df.meta$group <- as.character(df.meta$group)

# Selected features from DE

sel.gene.tb_od <- read.csv("../../data/ex_11/feat_sel_1_2/gene_tb_od_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
#sel.prot.tb_od <- read.csv("../../data/ex_11/feat_sel/prot_tb_od_BH_LFC_sig_factors.csv", header=TRUE, row.names = 1)
#sel.gp.tb_od <- rbind(sel.prot.tb_od, sel.gene.tb_od)


sel.gene.tb_ltbi <- read.csv("../../data/ex_11/feat_sel_1_2/gene_tb_od_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
#sel.prot.tb_ltbi <- read.csv("../../data/ex_11/feat_sel/prot_tb_od_BH_LFC_sig_factors.csv", header=TRUE, row.names = 1)
#sel.gp.tb_ltbi <- rbind(sel.prot.tb_ltbi, sel.gene.tb_ltbi)


sel.gene.tb_nontb <- read.csv("../../data/ex_11/feat_sel_1_2/gene_tb_nontb_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
#sel.prot.tb_nontb <- read.csv("../../data/ex_11/feat_sel/prot_tb_nontb_BH_LFC_sig_factors.csv", header=TRUE, row.names = 1)
#sel.gp.tb_nontb <- rbind(sel.prot.tb_nontb, sel.gene.tb_nontb)


# Create meta variant for easy tb vs non-tb identification w/in hiv status group
df.meta.tb_nontb <- df.meta
df.meta.tb_nontb$group <- as.numeric(df.meta.tb_nontb$group)
df.meta.tb_nontb$group[df.meta.tb_nontb$group == 3] <- 7
df.meta.tb_nontb$group[df.meta.tb_nontb$group == 6] <- 7
df.meta.tb_nontb$group[df.meta.tb_nontb$group == 4] <- 8
df.meta.tb_nontb$group[df.meta.tb_nontb$group == 5] <- 8
df.meta.tb_nontb$group <- as.character(df.meta.tb_nontb$group)

# To direct to the correct folder
date <- "2018-08-20/"
ex_dir <- "ex_11/"
feat_sel <- "feat_sel_1_2/"

# Parameters
alpha = 0.5 # Elastic net
K = 20

# comparison sets

comps <- list(
  list(sel.gene.tb_od, df.meta, 1, 6, "TB", "OD", "HIV-", "gp_tb_od_hivneg")
  ,
  list(sel.gene.tb_ltbi, df.meta, 1, 3, "TB", "LTBI", "HIV-", "gp_tb_ltbi_hivneg")
  ,
  list(sel.gene.tb_nontb, df.meta.tb_nontb, 1, 7, "TB", "non-TB", "HIV-", "gp_tb_nontb_hivneg")
)

for (comp in comps){
  comp.sel.feats <- comp[[1]]
  comp.meta <- comp[[2]]
  comp.g1 <- comp[[3]][1]
  comp.g2 <- comp[[4]][1]
  comp.verb.g1 <- comp[[5]][1]
  comp.verb.g2 <- comp[[6]][1]
  comp.hivstat <- comp[[7]][1]
  comp.abbrv <- comp[[8]][1]
  
  # Subset for correct patients and fetures for comparison
  
  ind.comp <- c()
  
  for (i in 1:nrow(comp.meta)){
    if((comp.meta$group[i] == comp.g1) || (comp.meta$group[i] == comp.g2)){
      ind.comp <- c(ind.comp, i)
    }
  }
  
  df.meta.comp <- comp.meta[ind.comp,]
  comp.sel.data.noprot <- df.gene.all[ind.comp, match(comp.sel.feats$features, colnames(df.gene.all))]
  
  # Add proteins
  comp.sel.data <- cbind(df.prot.all[ind.comp,], comp.sel.data.noprot)
  
  ##############################################
  # Use limma just to get logFC and beta coeff #
  ##############################################
  
  # Make factors for analysis to consider
  
  fac.sex <- factor(df.meta.comp$sex)
  fac.site <- factor(df.meta.comp$site) # Correct for site? Maybe use as intercept
  fac.tb <- factor(df.meta.comp$tb.status)
  
  design <- model.matrix(~fac.site + fac.sex + fac.tb)
  
  fit <- lmFit(t(comp.sel.data), design)
  fit <- eBayes(fit, trend=TRUE, robust=TRUE)
  results <- decideTests(fit)
  tab.res <- topTable(fit, coef=4, n=ncol(comp.sel.data))
  
   ####################
  ## EN selection ##
  ####################
  
  # Fit to EN
  
  # Cross-validated analysis of coefficients
  foldid <- sample(1:K,size=length(data.matrix(as.numeric(df.meta.comp$group))),replace=TRUE)
  
  cvfit <- cv.glmnet(data.matrix(comp.sel.data), 
                     data.matrix(as.numeric(df.meta.comp$group)), 
                     family="gaussian",
                     foldid = foldid,
                     alpha=alpha
  )
  
  png(paste("../../img/", ex_dir, date, comp.abbrv, "_BH_LFC_EN_cv_glmnet_coeff.png", sep=""),
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8        # smaller font size
  )
  
  plot(cvfit, main=paste("Cross-validated coefficient plot for HIV-", comp.verb.g1, "vs", comp.verb.g2, "gene and protein data", sep=" "))
  dev.off()
  
  # Check which are selected for in second round with elastic net with lambda + 1se
  
  #Convert to dataframe for easier sorting
  coef.cvfit.lambda1se <- coef(cvfit, s=cvfit$lambda.1se)
  
  features <- coef.cvfit.lambda1se@Dimnames[[1]][which(coef.cvfit.lambda1se != 0 ) ]
  
  df.coef.cvfit.lambda1se <- data.frame(
    features, 
    coefs <- coef.cvfit.lambda1se[which(coef.cvfit.lambda1se != 0 ) ] 
  )
  colnames(df.coef.cvfit.lambda1se) <- c("features", "coefficients")
  
  #Remove intercept
  if(df.coef.cvfit.lambda1se$features[1] == "(Intercept)"){
    df.coef.cvfit.lambda1se <- df.coef.cvfit.lambda1se[-1,]
    features <- features[-1]
  }
  
  # Attach beta coefficients and log fold changes
  beta <- tab.res[match(df.coef.cvfit.lambda1se$features, rownames(tab.res)),]$B
  logFC <- tab.res[match(df.coef.cvfit.lambda1se$features, rownames(tab.res)),]$logFC
  df.coef.cvfit.lambda1se <- cbind(df.coef.cvfit.lambda1se, data.frame(logFC, beta))
  
  # Add label desribing direction of expression
  reg_dir <- c()
  
  for (l in 1:length(logFC)){
    if (logFC[l] < 0){
      reg_dir <- c(reg_dir, "down")
    } else {
      reg_dir <- c(reg_dir, "up")
    }
  }
  
  df.coef.cvfit.lambda1se <- cbind(df.coef.cvfit.lambda1se, reg_dir)
  
  write.csv(df.coef.cvfit.lambda1se, paste("../../data/", ex_dir, feat_sel, comp.abbrv, "_BH_LFC_EN_sig_factors.csv", sep=""), row.names=1:nrow(df.coef.cvfit.lambda1se))
  
}


toc()