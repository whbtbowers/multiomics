setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_9/")

library(limma)
library(heatmap3)
library(SNFtool)
library(glmnet)

set.seed(12)

#df.gene.body <- read.csv("../../data/ex_9/gene_train_body.csv", header=TRUE, row.names = 1)  #Gene test/train set
df.prot.body <- read.csv("../../data/ex_9/prot_train_body.csv", header=TRUE, row.names = 1)  # Protein train set
#df.gp.body <- cbind(df.prot.body, df.gene.body)

df.meta <- read.csv("../../data/ex_9/gp_train_meta.csv", header=TRUE, row.names = 1)
df.meta$group <- as.character(df.meta$group)

df.meta.tb_nontb <- df.meta
df.meta.tb_nontb$group <- as.numeric(df.meta.tb_nontb$group)
df.meta.tb_nontb$group[df.meta.tb_nontb$group == 3] <- 7
df.meta.tb_nontb$group[df.meta.tb_nontb$group == 6] <- 7
df.meta.tb_nontb$group[df.meta.tb_nontb$group == 4] <- 8
df.meta.tb_nontb$group[df.meta.tb_nontb$group == 5] <- 8
df.meta.tb_nontb$group <- as.character(df.meta.tb_nontb$group)

# To direct to the correct folder
date <- "2018-08-15/"
ex_dir <- "ex_10/"
feat_sel <- "feat_sel/"

# Parameters
alpha = 1 # lasso
K = 20 # for KNN

datasets <- list(
  #list(df.gene.body, "gene", "gene"),
  list(df.prot.body, "prot", "prot")
)

comps <- list(
  list(df.meta, "tb_od", "TB", "OD", 1, 6),
  list(df.meta, "tb_ltbi", "TB", "LTBI", 1, 3),
  list(df.meta.tb_nontb, "tb_nontb", "TB", "non-TB", 1, 7)
)

for (dataset in datasets){
  set.data <- dataset[[1]]
  set.abrv <- dataset[[2]][[1]]
  set.verbose <- dataset[[3]][[1]]
  
  for (comp in comps){
    full.meta <- comp[[1]]
    comp.abrv <- comp[[2]][[1]]
    comp.verb.g1 <- comp[[3]][[1]]
    comp.verb.g2 <- comp[[4]][[1]]
    comp.g1 <- comp[[5]][[1]]
    comp.g2 <- comp[[6]][[1]]
    
    #Get patients and metadata for comparison
    
    ind.comp <- c()
    
    for (i in 1:nrow(set.data)){
      if((full.meta$group[i] == comp.g1) || (full.meta$group[i] == comp.g2)){
        ind.comp <- c(ind.comp, i)
      }
    }
    
    comp.data <- set.data[ind.comp, ]
    comp.meta <- full.meta[ind.comp, ]
    
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
      if ((abs(tab.res$logFC[k]) > 0.5) & (tab.res$adj.P.Val[k] < 0.05)){
        ind.dif_ex <- c(ind.dif_ex, k)
      }
    }
    
    sig_P = tab.res$adj.P.Val[ind.dif_ex]
    sig_factor = rownames(tab.res)[ind.dif_ex]
    sig_rows = tab.res[ind.dif_ex,]
    #write.csv(sig_rows, paste("../../data/", ex_dir, "feat_sel/",set.abrv, "_", comp.abrv, "_BH_LFC_sig_factors.csv", sep=""))
    
    #Select significant features
    sel.comp.data <- comp.data[match(sig_factor, colnames(comp.data))]
    
    ###########################
    ## Elastic net selection ##
    ###########################
    
    # Fit to elastic net
    # Feed in BH selected factors
    
    fit.glm <- glmnet(as.matrix(sel.comp.data), 
                      comp.meta$group, 
                      family="gaussian", 
                      alpha=alpha
    )
    
    png(paste("../../img/", ex_dir, date, set.abrv, "_", comp.abrv, "_BH_LFC_lasso_glmnet_coeff.png", sep=""),
        width = 5*300,        # 5 x 300 pixels
        height = 5*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8        # smaller font size
    )
    
    plot(fit.glm, xvar = "lambda", label=TRUE)
    dev.off()
    
    # Cross-validated analysis of coefficients
    foldid <- sample(1:K,size=length(data.matrix(as.numeric(comp.meta$group))),replace=TRUE)
    
    cvfit <- cv.glmnet(data.matrix(sel.comp.data), 
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
    
    write.csv(df.coef.cvfit.lambda1se, paste("../../data/", ex_dir, "feat_sel/",set.abrv, "_", comp.abrv, "_BH_LFC_lasso_sig_factors.csv", sep=""), row.names=1:nrow(df.coef.cvfit.lambda1se))
    
  }
}
