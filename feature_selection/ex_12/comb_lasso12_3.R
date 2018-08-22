setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_12/")

library(limma)
library(heatmap3)
library(SNFtool)
library(glmnet)
library(stringr)

set.seed(12)

df.gene.body <- read.csv("../../data/ex_12/gene_train_body.csv", header=TRUE, row.names = 1)  #Gene test/train set

df.prot.body <- read.csv("../../data/ex_12/prot_train_body.csv", header=TRUE, row.names = 1)  # Protein train set

df.meta <- read.csv("../../data/ex_12/gp_train_meta.csv", header=TRUE, row.names = 1)
df.meta$group <- as.character(df.meta$group)

# To direct to the correct folder
date <- "2018-08-21/"
ex_dir <- "ex_12/"

# Parameters
alpha = 1 # lasso
K = 20

# Select HIV- patients and tb status comaprisons

ind.hiv_neg.tb_ltbi <- c()
ind.hiv_neg.tb_od <- c()
ind.hiv_neg.tb_nontb <- c()

for (i in 1:nrow(df.meta)){
  if((df.meta$group[i] == 1) || (df.meta$group[i] == 3)){
    ind.hiv_neg.tb_ltbi <- c(ind.hiv_neg.tb_ltbi, i)
  }
}

df.gene.body.hiv_neg.tb_ltbi <- df.gene.body[ind.hiv_neg.tb_ltbi,]
df.prot.body.hiv_neg.tb_ltbi <- df.prot.body[ind.hiv_neg.tb_ltbi,]
df.meta.hiv_neg.tb_ltbi <- df.meta[ind.hiv_neg.tb_ltbi,]

# Select HIV- TB vs OD patients

for (i in 1:nrow(df.meta)){
  if((df.meta$group[i] == 1) || (df.meta$group[i] == 6)){
    ind.hiv_neg.tb_od <- c(ind.hiv_neg.tb_od, i)
  }
}

df.gene.body.hiv_neg.tb_od <- df.gene.body[ind.hiv_neg.tb_od,]
df.prot.body.hiv_neg.tb_od <- df.prot.body[ind.hiv_neg.tb_od,]
df.meta.hiv_neg.tb_od <- df.meta[ind.hiv_neg.tb_od,]

# Select HIV- TB vs nontb patients

for (i in 1:nrow(df.meta)){
  if((df.meta$group[i] == 1) || (df.meta$group[i] == 3) || (df.meta$group[i] == 6)){
    ind.hiv_neg.tb_nontb <- c(ind.hiv_neg.tb_nontb, i)
  }
}

df.gene.body.hiv_neg.tb_nontb <- df.gene.body[ind.hiv_neg.tb_nontb,]
df.prot.body.hiv_neg.tb_nontb <- df.prot.body[ind.hiv_neg.tb_nontb,]
df.meta.hiv_neg.tb_nontb <- df.meta[ind.hiv_neg.tb_nontb,]

#Change group of LTBI and OD to 7
df.meta.hiv_neg.tb_nontb$group[df.meta.hiv_neg.tb_nontb$group == 3] <- 7
df.meta.hiv_neg.tb_nontb$group[df.meta.hiv_neg.tb_nontb$group == 6] <- 7

comparisons <- list(
  list(df.gene.body.hiv_neg.tb_ltbi, df.prot.body.hiv_neg.tb_ltbi, df.meta.hiv_neg.tb_ltbi, "TB vs LTBI", "tb_ltbi")
  ,
  list(df.gene.body.hiv_neg.tb_od, df.prot.body.hiv_neg.tb_od, df.meta.hiv_neg.tb_od, "TB vs OD", "tb_od")
  ,
  list(df.gene.body.hiv_neg.tb_nontb, df.prot.body.hiv_neg.tb_nontb, df.meta.hiv_neg.tb_nontb, "TB vs nontb", "tb_nontb")
)

for (comparison in comparisons){
  comp.gene <- comparison[[1]]
  comp.prot <- comparison[[2]]
  comp.meta <- comparison[[3]]
  comp.verbose <- comparison[[4]][[1]]
  comp.abrv <- comparison[[5]][[1]]
  
  #############################
  ## Limma-based DE analysis ##
  #############################
  
  # Make factors for analysis to consider
  
  fac.sex <- factor(comp.meta$sex)
  fac.site <- factor(comp.meta$site) # Correct for site? Maybe use as intercept
  fac.tb <- factor(comp.meta$tb.status)
  
  design <- model.matrix(~fac.site + fac.sex + fac.tb)
  
  # run analysis for genes and select based on analysis
  
  fit.gene <- lmFit(t(comp.gene), design)
  fit.gene <- eBayes(fit.gene, trend=TRUE, robust=TRUE)
  results.gene <- decideTests(fit.gene)
  #print(summary(results.gene))
  tab.res.gene <- topTable(fit.gene, coef="fac.tbTB", n=ncol(comp.gene))
  #print(tab.res.gene)
  
  png(paste("../../img/", ex_dir, date, "gene_", comp.abrv, "_meandiff.png", sep=""),
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8        # smaller font size
  )
  
  plotMD(fit.gene, coef=4, status=results.gene["fac.tbTB"], values=c(1,-1), hl.col=c("red","blue"))
  dev.off()
  
  # Get factors by sig |LFC| and BH-corrected p-value
  
  ind.dif_ex.gene <- c()
  
  for (k in 1:nrow(tab.res.gene)){
    if ((abs(tab.res.gene$logFC[k]) > 0.5) & (tab.res.gene$adj.P.Val[k] < 0.02)){
      ind.dif_ex.gene <- c(ind.dif_ex.gene, k)
    }
  }
  
  sig_P.gene = tab.res.gene$adj.P.Val[ind.dif_ex.gene]
  sig_factor.gene = rownames(tab.res.gene)[ind.dif_ex.gene]
  sig_rows.gene = tab.res.gene[ind.dif_ex.gene,]
  write.csv(sig_rows.gene, paste("../../data/", ex_dir, "feat_sel/gene_", comp.abrv, "_BH_LFC_sig_factors.csv", sep=""))
  #print(comp.verbose)
  #print(nrow(sig_rows.gene))
  #Select significant features
  sel.comp.gene.body <- comp.gene[,match(sig_factor.gene, colnames(comp.gene))]
  
  # run analysis for proteins and select based on analysis
  
  fit.prot <- lmFit(t(comp.prot), design)
  fit.prot <- eBayes(fit.prot, trend=TRUE, robust=TRUE)
  results.prot <- decideTests(fit.prot)
  #print(summary(results.prot))
  tab.res.prot <- topTable(fit.prot, coef="fac.tbTB", n=ncol(comp.prot))
  #print(tab.res.prot)
  
  png(paste("../../img/", ex_dir, date, "prot_", comp.abrv, "_meandiff.png", sep=""),
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8        # smaller font size
  )
  
  plotMD(fit.prot, coef=4, status=results.prot["fac.tbTB"], values=c(1,-1), hl.col=c("red","blue"))
  dev.off()
  
  # Get factors by sig |LFC| and BH-corrected p-value
  
  ind.dif_ex.prot <- c()
  
  for (k in 1:nrow(tab.res.prot)){
    if ((abs(tab.res.prot$logFC[k]) > 0.5) & (tab.res.prot$adj.P.Val[k] < 0.05)){
      ind.dif_ex.prot <- c(ind.dif_ex.prot, k)
    }
  }
  
  sig_P.prot = tab.res.prot$adj.P.Val[ind.dif_ex.prot]
  sig_factor.prot = rownames(tab.res.prot)[ind.dif_ex.prot]
  sig_rows.prot = tab.res.prot[ind.dif_ex.prot,]
  write.csv(sig_rows.prot, paste("../../data/", ex_dir, "feat_sel/gene_", comp.abrv, "_BH_LFC_sig_factors.csv", sep=""))
  #print(comp.verbose)
  #print(nrow(sig_rows.prot))
  #Select significant features
  sel.comp.prot.body <- comp.prot[,match(sig_factor.prot, colnames(comp.prot))]
  
  
  # Concatenate prot and seleted gene diff ex tables and body data
  tab.res <- rbind(sig_rows.prot, sig_rows.gene)
  sel.comp.data <- cbind(sel.comp.prot.body, sel.comp.gene.body)
  
  #Run DE to get logFC and beta for ROC curves
  
  #fit.gp <- lmFit(t(sel.comp.data), design)
  #fit.gp <- eBayes(fit.gp, trend=TRUE, robust=TRUE)
  #results.gp <- decideTests(fit.gp)
  #print(summary(results.gp))
  #tab.res.gp <- topTable(fit.gp, coef="fac.tbTB", n=ncol(sel.comp.data))
  #print(tab.res.gp)
  
  #####################
  ## Lasso selection ##
  #####################
  
  foldid <- sample(1:20,size=length(data.matrix(as.numeric(comp.meta$group))), replace=TRUE)
  
  cvfit <- cv.glmnet(data.matrix(sel.comp.data), 
                     data.matrix(as.numeric(comp.meta$group)), 
                     family="gaussian",
                     foldid = foldid,
                     alpha=alpha
  )
  
  png(paste("../../img/", ex_dir, date, "gp_", comp.abrv, "_BH_LFC_lasso_cv_glmnet_coeff.png", sep=""),
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8        # smaller font size
  )
  
  plot(cvfit, main=paste("Cross-validated coefficient plot for HIV-", comp.verbose, "protein data", sep=" "))
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
  
  write.csv(df.coef.cvfit.lambda1se, paste("../../data/", ex_dir, "feat_sel_1_2/", "gp_", comp.abrv, "_BH_LFC_lasso_sig_factors.csv", sep=""), row.names=1:nrow(df.coef.cvfit.lambda1se))
  
  #Display number of selected factors
  
  num.sig.factors <- length(df.coef.cvfit.lambda1se$features)
  
  sig.factors.genecount <- 0
  ind.factors.gene <- c()
  
  for (i in 1:num.sig.factors){
    if (str_detect(toString(df.coef.cvfit.lambda1se$features[i]), "ILMN_")){
      sig.factors.genecount <- sig.factors.genecount + 1
      ind.factors.gene <- c(ind.factors.gene, i)
    }
  }
  
  sig.factors.protcount <- num.sig.factors - sig.factors.genecount
  
  #Write CSVs for individual features
  write.csv(df.coef.cvfit.lambda1se[ind.factors.gene,], paste("../../data/", ex_dir, "feat_sel_2/", "gene_", comp.abrv, "_BH_LFC_lasso_sig_factors.csv", sep=""), row.names=1:nrow(df.coef.cvfit.lambda1se[ind.factors.gene,]))
  write.csv(df.coef.cvfit.lambda1se[-ind.factors.gene,], paste("../../data/", ex_dir, "feat_sel_2/", "prot_", comp.abrv, "_BH_LFC_lasso_sig_factors.csv", sep=""), row.names=1:nrow(df.coef.cvfit.lambda1se[-ind.factors.gene,]))
  
  
  print(paste("For", comp.verbose, "after DE analysis,", nrow(sig_rows.gene), "gene probes selected and", nrow(sig_rows.gene), "proteins. After lasso,", num.sig.factors, "factors selected.", sig.factors.genecount, "genes and", sig.factors.protcount, "proteins.", sep = " "))
  
}
