setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_10/")

library(tictoc)
library(limma)
library(heatmap3)
library(SNFtool)
library(glmnet)

set.seed(12)

tic()

df.gene.body <- read.csv("../../data/ex_9/gene_train_body.csv", header=TRUE, row.names = 1)  #Gene test/train set

df.prot.body <- read.csv("../../data/ex_9/prot_train_body.csv", header=TRUE, row.names = 1)  # Protein train set

df.meta <- read.csv("../../data/ex_9/gp_train_meta.csv", header=TRUE, row.names = 1)
df.meta$group <- as.character(df.meta$group)

# To direct to the correct folder
date <- "2018-08-21/"
ex_dir <- "ex_12/"
feat_sel <- "feat_sel_1_2/"

# Parameters
alpha = 1 # lasso
K = 20

datasets = list(
  list(df.gene.body, df.meta, "gene", "gene")
  #,
  #list(df.prot.body, df.meta, "protein", "prot")
)

for (dataset in datasets){
  set.data <- dataset[[1]]
  set.meta <- dataset[[2]]
  set.verbose <- dataset[[3]][[1]]
  set.abrv <- dataset[[4]][[1]]
  
  # Select HIV- TB vs LTBI patients
  
  ind.hiv_neg.tb_ltbi <- c()
  ind.hiv_neg.tb_od <- c()
  ind.hiv_neg.tb_nontb <- c()
  
  for (i in 1:nrow(set.data)){
    if((set.meta$group[i] == 1) || (set.meta$group[i] == 3)){
      ind.hiv_neg.tb_ltbi <- c(ind.hiv_neg.tb_ltbi, i)
    }
  }
  
  set.data.hiv_neg.tb_ltbi <- set.data[ind.hiv_neg.tb_ltbi,]
  set.meta.hiv_neg.tb_ltbi <- set.meta[ind.hiv_neg.tb_ltbi,]
  
  # Select HIV- TB vs OD patients
  
  for (i in 1:nrow(set.data)){
    if((set.meta$group[i] == 1) || (set.meta$group[i] == 6)){
      ind.hiv_neg.tb_od <- c(ind.hiv_neg.tb_od, i)
    }
  }
  
  set.data.hiv_neg.tb_od <- set.data[ind.hiv_neg.tb_od,]
  set.meta.hiv_neg.tb_od <- set.meta[ind.hiv_neg.tb_od,]
  
  # Select HIV- TB vs nontb patients
  
  for (i in 1:nrow(set.data)){
    if((set.meta$group[i] == 1) || (set.meta$group[i] == 3) || (set.meta$group[i] == 6)){
      ind.hiv_neg.tb_nontb <- c(ind.hiv_neg.tb_nontb, i)
    }
  }
  
  set.data.hiv_neg.tb_nontb <- set.data[ind.hiv_neg.tb_nontb,]
  set.meta.hiv_neg.tb_nontb <- set.meta[ind.hiv_neg.tb_nontb,]
  
  #Change group of LTBI and OD to 7
  set.meta.hiv_neg.tb_nontb$group[set.meta.hiv_neg.tb_nontb$group == 3] <- 7
  set.meta.hiv_neg.tb_nontb$group[set.meta.hiv_neg.tb_nontb$group == 6] <- 7
  
  comparisons <- list(
    list(set.data.hiv_neg.tb_ltbi, set.meta.hiv_neg.tb_ltbi, "TB vs LTBI", "tb_ltbi")
    ,
    list(set.data.hiv_neg.tb_od, set.meta.hiv_neg.tb_od, "TB vs OD", "tb_od")
    ,
    list(set.data.hiv_neg.tb_nontb, set.meta.hiv_neg.tb_nontb, "TB vs nontb", "tb_nontb")
  )
  
  for (j in 1:length(comparisons)){
    comp.data <- comparisons[[j]][[1]]
    comp.meta <- comparisons[[j]][[2]]
    comp.verbose <- comparisons[[j]][[3]]
    comp.abrv<- comparisons[[j]][[4]]
    
    #############################
    ## Limma-based DE analysis ##
    #############################
    
    # Make factors for analysis to consider
    
    fac.sex <- factor(comp.meta$sex)
    fac.site <- factor(comp.meta$site) # Correct for site? Maybe use as intercept
    fac.tb <- factor(comp.meta$tb.status)
    
    design <- model.matrix(~fac.site + fac.sex + fac.tb)
    
    fit <- lmFit(t(comp.data), design)
    fit <- eBayes(fit, trend=TRUE, robust=TRUE)
    results <- decideTests(fit)
    print(summary(results))
    tab.res <- topTable(fit, coef=4, n=ncol(comp.data))
    print(tab.res)
    
    png(paste("../../img/", ex_dir, date, set.abrv, "_", comp.abrv, "_meandiff.png", sep=""),
        width = 5*300,        # 5 x 300 pixels
        height = 5*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8        # smaller font size
    )
    
    plotMD(fit, coef=4, status=results[,4], values=c(1,-1), hl.col=c("red","blue"))
    dev.off()
    
    # Get factors by sig |LFC| and BH-corrected p-value
    
    ind.dif_ex <- c()
    
    
    for (k in 1:length(tab.res$logFC)){
      if ((abs(tab.res$logFC[k]) > 0.5) & (tab.res$adj.P.Val[k] < 0.03)){
    ind.dif_ex <- c(ind.dif_ex, k)
      }
    }
    
    sig_P = tab.res$adj.P.Val[ind.dif_ex]
    sig_factor = rownames(tab.res)[ind.dif_ex]
    sig_rows = tab.res[ind.dif_ex,]
    write.csv(sig_rows, paste("../../data/", ex_dir, feat_sel, set.abrv, "_", comp.abrv, "_BH_LFC_sig_factors.csv", sep=""))
    
    #Select significant features
    sel.comp.data <- comp.data[,match(sig_factor, colnames(comp.data))]
    
    # If only one feature, stop here
    
    if (length(colnames(sel.comp.data))==1){
      beta <- tab.res[match(colnames(sel.comp.data), rownames(tab.res)),]$B
      logFC <- tab.res[match(colnames(sel.comp.data), rownames(tab.res)),]$logFC
      
      # Add label desribing direction of expression
      reg_dir <- c()
      
      for (l in 1:length(logFC)){
        if (logFC[l] < 0){
          reg_dir <- c(reg_dir, "down")
        } else {
          reg_dir <- c(reg_dir, "up")
        }
      }
      
      features <- colnames(sel.comp.data)
      
      df.sel_feat <- data.frame(features, beta, logFC, reg_dir)
      
      write.csv(df.sel_feat, paste("../../data/", ex_dir, "feat_sel/", set.abrv, "_", comp.abrv, "_BH_LFC_lasso_sig_factors.csv", sep=""), row.names=1:nrow(df.sel_feat))
      
    } else {
    
    
      ###########################
      ## Lasso selection ##
      ###########################
      
      # Fit to lasso model
      # Feed in BH selected factors and all proteins
      
      # Cross-validated analysis of coefficients
      foldid <- sample(1:K,size=length(data.matrix(as.numeric(comp.meta$group))),replace=TRUE)
      
      cvfit <- cv.glmnet(data.matrix(gp.sel.comp.data), 
                         data.matrix(as.numeric(comp.meta$group)), 
                         family="gaussian",
                         foldid = foldid,
                         alpha=alpha
      )
      
      png(paste("../../img/", ex_dir, date, set.abrv, "_", comp.abrv, "_BH_LFC_lasso_cv_glmnet_coeff.png", sep=""),
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
      
      write.csv(df.coef.cvfit.lambda1se, paste("../../data/", ex_dir, "feat_sel/", set.abrv, "_", comp.abrv, "_BH_LFC_lasso_sig_factors.csv", sep=""), row.names=1:nrow(df.coef.cvfit.lambda1se))
      
      
      }
    }
}

# Display time elapsed
toc()