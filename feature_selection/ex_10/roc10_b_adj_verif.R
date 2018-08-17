setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_9/")
#setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/ex_9/")

library(pROC)
library(ggplot2)

set.seed(12)

# To direct to the correct folder
date <- "2018-08-18/"
ex_dir <- "ex_10/"

# Selected features for tb vs od
sel.gene.tb_od <- read.csv("../../data/ex_10/feat_sel_2/gene_tb_od_hivneg_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.prot.tb_od <- read.csv("../../data/ex_10/feat_sel_2/prot_tb_od_hivneg_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.gp.tb_od <- read.csv("../../data/ex_10/feat_sel_1_2/gp_tb_od_hivneg_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)


# Selected features for tb vs ltbi
sel.gene.tb_ltbi <- read.csv("../../data/ex_10/feat_sel_2/gene_tb_ltbi_hivneg_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.prot.tb_ltbi <- read.csv("../../data/ex_10/feat_sel_2/prot_tb_ltbi_hivneg_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.gp.tb_ltbi <- read.csv("../../data/ex_10/feat_sel_1_2/gp_tb_ltbi_hivneg_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)


# Selected features for tb vs non-tb
sel.gene.tb_nontb <- read.csv("../../data/ex_10/feat_sel_2/gene_tb_nontb_hivneg_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.prot.tb_nontb <- read.csv("../../data/ex_10/feat_sel_2/prot_tb_nontb_hivneg_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.gp.tb_nontb <- read.csv("../../data/ex_10/feat_sel_1_2/gp_tb_nontb_hivneg_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)

# Complete datasets
df.gene.val <- read.csv("../../data/ex_9/gene_validation_body.csv", header=TRUE, row.names = 1)  #Gene validation set
df.prot.val <- read.csv("../../data/ex_9/prot_validation_body.csv", header=TRUE, row.names = 1)  #Gene validation set
df.gp.val <- cbind(df.prot.val, df.gene.val)
df.meta.val <- read.csv("../../data/ex_9/gp_validation_meta.csv", header=TRUE, row.names = 1)
df.meta.val$group <- as.character(df.meta.val$group)

# Create meta variant for easy tb vs non-tb identification w/in hiv status group
df.meta.val.tb_nontb <- df.meta.val
df.meta.val.tb_nontb$group <- as.numeric(df.meta.val.tb_nontb$group)
df.meta.val.tb_nontb$group[df.meta.val.tb_nontb$group == 3] <- 7
df.meta.val.tb_nontb$group[df.meta.val.tb_nontb$group == 6] <- 7
df.meta.val.tb_nontb$group[df.meta.val.tb_nontb$group == 4] <- 8
df.meta.val.tb_nontb$group[df.meta.val.tb_nontb$group == 5] <- 8
df.meta.val.tb_nontb$group <- as.character(df.meta.val.tb_nontb$group)

# Create meta variant for easy tb vs other factor identification between hiv status groups
df.meta.val.allhiv <- df.meta.val
df.meta.val.allhiv$group <- as.numeric(df.meta.val.allhiv$group)
df.meta.val.allhiv$group[df.meta.val.allhiv$group == 1] <- 9
df.meta.val.allhiv$group[df.meta.val.allhiv$group == 2] <- 9
df.meta.val.allhiv$group[df.meta.val.allhiv$group == 3] <- 10
df.meta.val.allhiv$group[df.meta.val.allhiv$group == 4] <- 10
df.meta.val.allhiv$group[df.meta.val.allhiv$group == 5] <- 11
df.meta.val.allhiv$group[df.meta.val.allhiv$group == 6] <- 11
df.meta.val.allhiv$group <- as.character(df.meta.val.allhiv$group)

# Create meta variant for easy tb vs non-tb identification between hiv status groups
df.meta.val.allhiv.tb_nontb <- df.meta.val
df.meta.val.allhiv.tb_nontb$group <- as.numeric(df.meta.val.allhiv.tb_nontb$group)
df.meta.val.allhiv.tb_nontb$group[df.meta.val.allhiv.tb_nontb$group == 1] <- 9
df.meta.val.allhiv.tb_nontb$group[df.meta.val.allhiv.tb_nontb$group == 2] <- 9
df.meta.val.allhiv.tb_nontb$group[df.meta.val.allhiv.tb_nontb$group == 3] <- 12
df.meta.val.allhiv.tb_nontb$group[df.meta.val.allhiv.tb_nontb$group == 4] <- 12
df.meta.val.allhiv.tb_nontb$group[df.meta.val.allhiv.tb_nontb$group == 5] <- 12
df.meta.val.allhiv.tb_nontb$group[df.meta.val.allhiv.tb_nontb$group == 6] <- 12
df.meta.val.allhiv.tb_nontb$group <- as.character(df.meta.val.allhiv.tb_nontb$group)

comps <- list(
  list(sel.gene.tb_od, sel.prot.tb_od, sel.gp.tb_od, df.meta.val, 1, 6, "TB", "OD", "tb_od")
  ,
  list(sel.gene.tb_ltbi, sel.prot.tb_ltbi, sel.gp.tb_ltbi, df.meta.val, 1, 3, "TB", "LTBI", "tb_ltbi")
  ,
  list(sel.gene.tb_nontb, sel.prot.tb_nontb, sel.gp.tb_nontb, df.meta.val.tb_nontb, 1, 7, "TB", "non-TB", "tb_nontb")
) 

for (comp in comps){
  comp.gene <- comp[[1]]
  comp.prot <- comp[[2]]
  comp.all <- comp[[3]]
  comp.meta <- comp[[4]]
  comp.g1 <- comp[[5]][[1]]
  comp.g2 <- comp[[6]][[1]]
  comp.verb.g1 <- comp[[7]][[1]]
  comp.verb.g2 <- comp[[8]][[1]]
  comp.abbrv <- comp[[9]][[1]]
  
  #Get patient in comparison group
  
  ind.comp <- c()
  
  for (i in 1:nrow(comp.meta)){
    if ((comp.meta$group[i] == comp.g1) || (comp.meta$group[i] == comp.g2)){
      ind.comp <- c(ind.comp, i)
      
    }
  }
  
  df.comp.gene <- df.gene.val[ind.comp,]
  df.comp.prot <- df.prot.val[ind.comp,]
  df.comp.gp <- df.gp.val[ind.comp,]
  df.comp.meta <- comp.meta[ind.comp,]
  
  #Get ROC curves for all proteins and individual proteins
  
  for (i in 1:length(comp.prot$features)){
    feat <- df.comp.prot[,match(sel.prot.tb_nontb$features[i], colnames(df.comp.prot))]
    
    feat.roc <- roc(df.comp.meta$group, feat, auc=TRUE)
    
    #Finish tomorrow
    if (feat.roc$auc == 1){}
    
    feat.roc.smooth <- smooth(feat.roc)
    
    png(paste("../../img/", ex_dir, date, comp.abbrv, "_prot_indv_roc", i, ".png", sep=""),
        width = 5*300,        # 5 x 300 pixels
        height = 5*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8) # smaller font size
    
    plot(feat.roc,
         col="red",
         lwd=1,
         main=paste(comp.prot$features[i], "," , comp.verb.g1,  "vs", comp.verb.g2),
         legacy.axes = TRUE
    )
    
    lines.roc(feat.roc.smooth, col="blue", lwd=3)
    legend("bottomright",
           title = paste("AUC =", format(round(feat.roc$auc, 2), nsmall = 2)),
           legend = c("Empirical data", "Fit"),
           col = c("red", "blue"),
           lwd = c(1, 3)
    )
    
    dev.off()
    
  }
  
}



for (i in 1:length(sel.prot.tb_nontb$features)){
  feat <- df.prot.val[,match(comp.prot$features[i], colnames(df.prot.val))]
  
  feat.roc <- roc(df.comp.meta$group, feat, auc=TRUE)
  feat.roc.smooth <- smooth(feat.roc)
  
  png(paste("../../img/", ex_dir, date, comp.abbrv, "_prot_indv_roc", i, ".png", sep=""),
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8) # smaller font size
  
  plot(drs.roc,
       col="red",
       lwd=1,
       main=paste(comp.prot$features[i], "," , comp.verb.g1,  "vs", comp.verb.g2),
       legacy.axes = TRUE
  )
  
  lines.roc(drs.roc.smooth, col="blue", lwd=3)
  legend("bottomright",
         title = paste("AUC =", format(round(drs.roc$auc, 2), nsmall = 2)),
         legend = c("Empirical data", "Fit"),
         col = c("red", "blue"),
         lwd = c(1, 3)
  )
  
  dev.off()
  
  print(feat)
}
