setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_9/")

library(limma)
library(heatmap3)
library(SNFtool)
library(glmnet)

set.seed(12)

df.gene.body <- read.csv("../../data/ex_9/gene_train_body.csv", header=TRUE, row.names = 1)  # Protein test/train set

df.prot.body <- read.csv("../../data/ex_9/prot_train_body.csv", header=TRUE, row.names = 1)  # Protein train set

df.meta <- read.csv("../../data/ex_9/gp_train_meta.csv", header=TRUE, row.names = 1)
df.meta$group <- as.character(df.meta$group)

# To direct to the correct folder
date <- "2018-07-31/"
ex_dir <- "ex_9/"

# Parameters
#alphas = c(0, 0.5, 1)
alpha = 0.5
K = 20

datasets = list(
  #list(df.gene.body, df.meta, "gene", "gene")
  #,
  list(df.prot.body, df.meta, "protein", "prot")
)

for (i in 1:length(datasets)){
  set.data <- datasets[[i]][[1]]
  set.meta <- datasets[[i]][[2]]
  set.verbose <- datasets[[i]][[3]]
  set.abrv <- datasets[[i]][[4]]
  
  #Select HIV- TB vs LTBI patients
  
  ind.hiv_neg.tb_ltbi <- c()
  ind.hiv_neg.tb_od <- c()
  
  for (i in 1:nrow(set.data)){
    if((set.meta$group[i] == 1) || (set.meta$group[i] == 3)){
      ind.hiv_neg.tb_ltbi <- c(ind.hiv_neg.tb_ltbi, i)
    }
  }
  
  for (i in 1:nrow(set.data)){
    if((set.meta$group[i] == 1) || (set.meta$group[i] == 5)){
      ind.hiv_neg.tb_od <- c(ind.hiv_neg.tb_od, i)
    }
  }
  
  set.data.hiv_neg.tb_ltbi <- set.data[ind.hiv_neg.tb_ltbi,]
  set.meta.hiv_neg.tb_ltbi <- set.meta[ind.hiv_neg.tb_ltbi,]
  
  set.data.hiv_neg.tb_od <- set.data[ind.hiv_neg.tb_od,]
  set.meta.hiv_neg.tb_od <- set.meta[ind.hiv_neg.tb_od,]
  
  comparisons <- list(
    list(set.data.hiv_neg.tb_ltbi, set.meta.hiv_neg.tb_ltbi, "TB vs LTBI", "tb_ltbi")
    ,
    list(set.data.hiv_neg.tb_od, set.meta.hiv_neg.tb_od, "TB vs OD", "tb_od")
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
    tab.res <- topTable(fit, coef=4, n=22)
    print(tab.res)
    
    png(paste("../../img/", ex_dir, date, set.abrv, "_", comp.abrv, "_meandiff.png", sep=""),
        width = 5*300,        # 5 x 300 pixels
        height = 5*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8        # smaller font size
    )
    
    plotMD(fit, coef=4, status=results[,4], values=c(1,-1), hl.col=c("red","blue"))
    dev.off()
    
    # Get significant BH corrected values
    
    ind.dif_ex <- c()
    
    
    for (i in 1:length(tab.res$adj.P.Val)){
      if (tab.res$adj.P.Val[i] < 0.05){
    ind.dif_ex <- c(ind.dif_ex, i)
      }
    }
    
    sig_P = tab.res$adj.P.Val[ind.dif_ex]
    sig_factor = rownames(tab.res)[ind.dif_ex]
    sig_rows = tab.res[ind.dif_ex,]
    write.csv(sig_rows, paste("../../data/", ex_dir, "feat_sel/",set.abrv, "_", comp.abrv, "_BH_sig_factors.csv", sep=""))
    
    ###########################
    ## Elastic net selection ##
    ###########################
    
    # Fit to elastic net
    
    fit.glm <- glmnet(as.matrix(comp.data), 
                      comp.meta$group, 
                      family="gaussian", 
                      alpha=alpha
    )
    
    png(paste("../../img/", ex_dir, date, set.abrv, "_", comp.abrv, "_glmnet_coeff.png", sep=""),
        width = 5*300,        # 5 x 300 pixels
        height = 5*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8        # smaller font size
    )
    
    plot(fit.glm, label=TRUE)
    dev.off()
    
    # Cross-validated analysis of coefficients
    
    cvfit <- cv.glmnet(data.matrix(comp.data), 
                       data.matrix(as.numeric(comp.meta$group)), 
                       family="gaussian",
                       alpha=alpha
    )
    
    png(paste("../../img/", ex_dir, date, set.abrv, "_", comp.abrv, "_cv_glmnet_coeff.png", sep=""),
        width = 5*300,        # 5 x 300 pixels
        height = 5*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8        # smaller font size
    )
    
    plot(cvfit, main=paste("Cross-validated coefficient plot for HIV-", comp.verbose, "protein data", sep=" "))
    dev.off()
    
    ######################################################
    ## Out of curiosity, comparing rige, EMN, and lasso ##
    ######################################################
    
    foldid=sample(1:K,size=length(data.matrix(as.numeric(comp.meta$group))),replace=TRUE)
    cv1=cv.glmnet(data.matrix(comp.data),data.matrix(as.numeric(comp.meta$group)),foldid=foldid,alpha=1)
    cv.5=cv.glmnet(data.matrix(comp.data),data.matrix(as.numeric(comp.meta$group)),foldid=foldid,alpha=.5)
    cv0=cv.glmnet(data.matrix(comp.data),data.matrix(as.numeric(comp.meta$group)),foldid=foldid,alpha=0)
    
    png(paste("../../img/", ex_dir, date, set.abrv, "_", comp.abrv, "_alpha_comp.png", sep=""),
        width = 5*300,        # 5 x 300 pixels
        height = 5*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8        # smaller font size
    )
    
    par(mfrow=c(2,2))
    plot(cv1);plot(cv.5);plot(cv0)
    plot(log(cv1$lambda),cv1$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=cv1$name)
    points(log(cv.5$lambda),cv.5$cvm,pch=19,col="grey")
    points(log(cv0$lambda),cv0$cvm,pch=19,col="blue")
    legend("topleft",legend=c("alpha= 1","alpha= .5","alpha 0"),pch=19,col=c("red","grey","blue"))
    dev.off()
  }
}
