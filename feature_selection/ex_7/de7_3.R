setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_7/")

library(limma)
library(heatmap3)
library(SNFtool)
library(glmnet)

set.seed(12)

df.prot.tt <- read.csv("../../data/ex_7/prot_testtrain_body.csv", header=TRUE, row.names = 1)  # Protein test/train set

df.gene.tt <- read.csv("../../data/ex_7/gene_testtrain_body.csv", header=TRUE, row.names = 1)  # Gene test/train set

df.all.meta.tt <- read.csv("../../data/ex_7/gene_prot_testtrain_meta.csv", header=TRUE, row.names = 1)
df.all.meta.tt$group <- as.character(df.all.meta.tt$group)

# To direct to the correct folder
date <- "2018-07-20/"
ex_dir <- "ex_7/"

# Parameters
#alphas = c(0, 0.5, 1)
alpha = 0.5
K = 20

datasets = list(
  list(df.gene.tt, df.all.meta.tt, "gene", "gene")
  ,
  list(df.prot.tt, df.all.meta.tt, "protein", "prot")
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
    
    diff.comp.data <- as.matrix(dist(comp.data))
    
    # Get affinity heatmap for each regularisation method
    
    
    affMat.comp.data <- affinityMatrix(diff.comp.data, K, alpha)
    
    png(paste("../../img/", ex_dir, date, set.abrv, "_affinity_", comp.abrv, "_alpha", alpha,".png", sep=""),
        width = 3000,        # 5 x 300 pixels
        height = 3000,
        res = 300,            # 300 pixels per inch
        pointsize = 8)        # smaller font size
    
    heatmap3(affMat.comp.data,
             margins=c(10,5),
             main= paste(set.verbose, "data, HIV-,",  comp.verbose, sep=" ")
    )
    
    dev.off()
      
    
    
    # Analyse for differential expression
    tb_status.factor <- factor(comp.meta$tb.status)
    site.factor <- factor(comp.meta$site)
    sex.factor <- factor(comp.meta$sex)
    
    mat.design <- model.matrix(~1 + tb_status.factor + sex.factor + site.factor)
    
    fit <- lmFit(t(comp.data), mat.design)
    fit <- eBayes(fit, trend=TRUE, robust=TRUE)
    results <- decideTests(fit)
    cat("\n\n\n", "HIV-", set.verbose, "HIV-", comp.verbose, "summary table\n\n")
    print(summary(results))
    
    # Plot fold-changes
    png(paste("../../img/", ex_dir, date, set.abrv, "_", comp.abrv, "_diff_ex.png", sep=""),
        width = 5*300,        # 5 x 300 pixels
        height = 5*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8        # smaller font size
    )
    
    plotMD(fit,status=results[,4],values=c(1,-1),hl.col=c("red","blue"))
    dev.off()
    
    #tab.res <- topTable(fit, coef="TB", n=ncol(comp.data))
    tab.res <- topTable(fit, coef = "tb_status.factorTB", n=ncol(comp.data))
    #print(tab.res)
    
    # Select by BH corrected P-values
    ind.dif_ex <- c()
    
    
    for (i in 1:length(tab.res$adj.P.Val)){
      if (tab.res$adj.P.Val[i] > 0.05){
        ind.dif_ex <- c(ind.dif_ex, i)
      }
    }
    
    sig_P = tab.res$adj.P.Val[ind.dif_ex]
    sig_factor = rownames(tab.res)[ind.dif_ex]
    
    sel.tab.res <- topTable(fit, coef = "tb_status.factorTB", n=length(sig_P))
    
    #cat("\n\n\n", "HIV-", set.verbose, "HIV-", comp.verbose, "factors with significant p-values\n\n")
    print(sel.tab.res)
    cat("\n\n")
    print(paste("Number of significant factors:", length(sig_factor), sep=" "))
    
    
      
    fit.glm <- glmnet(as.matrix(comp.data), 
                           comp.meta$group, 
                           family="gaussian", 
                           alpha=alpha
    )
    
    png(paste("../../img/", ex_dir, date, set.abrv, "_", comp.abrv,  "_coefficients_alpha", alpha, ".png", sep=""),
        width = 5*300,        # 5 x 300 pixels
        height = 5*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8        # smaller font size
    )
    plot(fit.glm, label=TRUE)
    dev.off()
    
    # Cross validate to find lamda cutoff
    
    cvfit <- cv.glmnet(data.matrix(comp.data), 
                       data.matrix(as.numeric(comp.meta$group)), 
                       family="gaussian",
                       alpha=alpha
    )
    
    #Get lamda values
    lamda.min <- cvfit$lambda.min
    lamda.1se <- cvfit$lambda.1se
    lamda.0.5se <- (lamda.1se + lamda.min)/2
    
    png(paste("../../img/", ex_dir, date, set.abrv, "_", comp.abrv, "_cv_lm_alpha", alpha, ".png", sep=""),
        width = 5*300,        # 5 x 300 pixels
        height = 5*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8         # smaller font size
    )
    
    plot(cvfit)
    dev.off()
      
    
  }
}


