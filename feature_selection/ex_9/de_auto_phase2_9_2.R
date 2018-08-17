setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_9/")

library(tictoc)
library(limma)
library(glmnet)

set.seed(12)

tic()

# Combined protein and gene probe set
#df.gp.all <- read.csv("../../data/ex_9/gp_train_body.csv", header=TRUE, row.names = 1)  #Gene test/train set

df.meta <- read.csv("../../data/ex_9/gp_train_meta.csv", header=TRUE, row.names = 1)
df.meta$group <- as.character(df.meta$group)

# Create meta variant for easy tb vs non-tb identification w/in hiv status group
df.meta.tb_nontb <- df.meta
df.meta.tb_nontb$group <- as.numeric(df.meta.tb_nontb$group)
df.meta.tb_nontb$group[df.meta.tb_nontb$group == 3] <- 7
df.meta.tb_nontb$group[df.meta.tb_nontb$group == 6] <- 7
df.meta.tb_nontb$group[df.meta.tb_nontb$group == 4] <- 8
df.meta.tb_nontb$group[df.meta.tb_nontb$group == 5] <- 8
df.meta.tb_nontb$group <- as.character(df.meta.tb_nontb$group)

# Selected features for tb vs od

sel.gp.tb_od <- read.csv("../../data/ex_9/feat_sel/gp_tb_od_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)

# Selected features for tb vs ltbi

sel.gp.tb_ltbi <- read.csv("../../data/ex_9/feat_sel/gp_tb_ltbi_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)

# Selected features for tb vs non-TB

sel.gp.tb_nontb <- read.csv("../../data/ex_9/feat_sel/gp_tb_nontb_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)

# To direct to the correct folder
date <- "2018-08-09/"
ex_dir <- "ex_9/"
feat_sel <- "feat_sel_2/"

# Parameters
alpha = 1 # lasso
#alpha = 0.5 # elastic net
K = 20

phase1s = list(
  list(sel.gp.tb_od, df.meta, 1, 6, "TB", "OD", "HIV-", "tb_od_hivneg")
  ,
  list(sel.gp.tb_ltbi, df.meta, 1, 3, "TB", "LTBI", "HIV-", "tb_ltbi_hivneg")
  ,
  list(sel.gp.tb_nontb, df.meta.tb_nontb, 1, 7, "TB", "non-TB", "HIV-", "tb_nontb_hivneg")
)

for (phase1 in phase1s){
  phase1.sel_feat <- phase1[[1]]
  phase1.meta <- phase1[[2]]
  phase1.g1 <- phase1[[3]][1]
  phase1.g2 <- phase1[[4]][1]
  phase1.verb.g1 <- phase1[[5]][1]
  phase1.verb.g2 <- phase1[[6]][1]
  phase1.hivstat <- phase1[[7]][1]
  phase1.abbrv <- phase1[[8]][1]
  
  
  # Select for correct patients for comparison
  
  ind.comp <- c()
  
  for (i in 1:nrow(df.gp.all)){
    if((phase1.meta$group[i] == phase1.g1) || (phase1.meta$group[i] == phase1.g2)){
      ind.comp <- c(ind.comp, i)
    }
  }
  
  df.meta.comp <- phase1.meta[ind.comp,]
  df.data.hivneg <- df.gp.all[ind.comp,]
  
  set.sel.data <- df.data.hivneg[, match(phase1.sel_feat$features, colnames(df.data.hivneg))]
  
  #############################
  ## Limma-based DE analysis ##
  #############################
  
  # Make factors for analysis to consider
  
  fac.sex <- factor(df.meta.comp$sex)
  fac.site <- factor(df.meta.comp$site) # Correct for site? Maybe use as intercept
  fac.tb <- factor(df.meta.comp$tb.status)
  
  design <- model.matrix(~fac.site + fac.sex + fac.tb)
  
  fit <- lmFit(t(set.sel.data), design)
  fit <- eBayes(fit, trend=TRUE, robust=TRUE)
  results <- decideTests(fit)
  print(summary(results))
  tab.res <- topTable(fit, coef=4, n=ncol(set.sel.data))
  print(tab.res)
  
  png(paste("../../img/", ex_dir, date, phase1.abbrv, "_phase2_meandiff.png", sep=""),
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8        # smaller font size
  )
  
  plotMD(fit, coef=4, status=results[,4], values=c(1,-1), hl.col=c("red","blue"), main=paste("Differential expresson of gene and protein data in", phase1.hivstat, "patients,", phase1.verb.g1, "vs", phase1.verb.g2, sep=" "))
  dev.off()
  
  # Get factors by sig |LFC| and BH-corrected p-value
  
  ind.dif_ex <- c()
  
  for (k in 1:length(tab.res$logFC)){
    if ((abs(tab.res$logFC[k]) > 0.5) & (tab.res$adj.P.Val[k] < 0.05)){
      ind.dif_ex <- c(ind.dif_ex, k)
    }
  }
  
  sig_P = tab.res$adj.P.Val[ind.dif_ex]
  sig_factor = rownames(tab.res)[ind.dif_ex]
  sig_rows = tab.res[ind.dif_ex,]
  write.csv(sig_rows, paste("../../data/", ex_dir, feat_sel , phase1.abbrv, "_BH_LFC_phase2_sig_factors.csv", row.namessep=""))
  
  #Select significant features
  sel.comp.data <- df.data.hivneg[match(sig_factor, colnames(df.data.hivneg))]
 
  ###########################
  ## Lasso selection ##
  ###########################
  
  # Fit to elastic net
  # Feed in selected factors
  
  #fit.glm <- glmnet(as.matrix(sel.comp.data), 
  #                  df.meta.comp$group, 
  #                  family="gaussian", 
  #                  alpha=alpha
  #)
  
  #png(paste("../../img/", ex_dir, date, phase1.abbrv, "_BH_LFC_phase2_lasso_glmnet_coeff.png", sep=""),
  #    width = 5*300,        # 5 x 300 pixels
  #    height = 5*300,
  #    res = 300,            # 300 pixels per inch
  #    pointsize = 8        # smaller font size
  #)
  
  #plot(fit.glm, xvar = "lambda", label=TRUE)
  #dev.off()
  
  # Cross-validated analysis of coefficients
  
  foldid <- sample(1:K,size=length(data.matrix(as.numeric(df.meta.comp$group))),replace=TRUE)
  
  cvfit <- cv.glmnet(data.matrix(sel.comp.data), 
                     data.matrix(as.numeric(df.meta.comp$group)), 
                     family="gaussian",
                     foldid = foldid,
                     alpha=alpha
  )
  
  png(paste("../../img/", ex_dir, date, phase1.abbrv, "_BH_LFC_phase2_lasso_cv_glmnet_coeff.png", sep=""),
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8        # smaller font size
  )
  
  plot(cvfit, main=paste("Cross-validated coefficient plot for gene protein data for", phase1.hivstat, "patients,", phase1.verb.g1, "vs", phase1.verb.g2, sep=" "))
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
  
  # Attach beta coefficients
  beta <- tab.res[match(df.coef.cvfit.lambda1se$features, rownames(tab.res)),]$B
  logFC <- tab.res[match(df.coef.cvfit.lambda1se$features, rownames(tab.res)),]$logFC
  df.coef.cvfit.lambda1se <- cbind(df.coef.cvfit.lambda1se, beta)
  
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
  
  write.csv(df.coef.cvfit.lambda1se, paste("../../data/", ex_dir, feat_sel, phase1.abbrv, "_BH_LFC_lasso_phase2_sig_factors.csv", sep=""), row.names=1:nrow(df.coef.cvfit.lambda1se))
  
}
