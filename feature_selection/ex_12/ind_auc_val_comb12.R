setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_12/")
#setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/ex_9/")

library(pROC)
library(ggplot2)
library(stringr)

set.seed(12)

# To direct to the correct folder
date <- "2018-08-30/"
ex_dir <- "ex_12/"

# Selected features for tb vs od
sel.gene.tb_od <- read.csv("../../data/ex_12/feat_sel_2/gene_tb_od_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.prot.tb_od <- read.csv("../../data/ex_12/feat_sel_2/prot_tb_od_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.gp.tb_od <- read.csv("../../data/ex_12/feat_sel_1_2/gp_tb_od_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)


# Selected features for tb vs ltbi
sel.gene.tb_ltbi <- read.csv("../../data/ex_12/feat_sel_2/gene_tb_ltbi_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.prot.tb_ltbi <- read.csv("../../data/ex_12/feat_sel_2/prot_tb_ltbi_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.gp.tb_ltbi <- read.csv("../../data/ex_12/feat_sel_1_2/gp_tb_ltbi_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)


# Selected features for tb vs non-tb
sel.gene.tb_nontb <- read.csv("../../data/ex_12/feat_sel_2/gene_tb_nontb_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.prot.tb_nontb <- read.csv("../../data/ex_12/feat_sel_2/prot_tb_nontb_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.gp.tb_nontb <- read.csv("../../data/ex_12/feat_sel_1_2/gp_tb_nontb_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)

# Complete datasets
df.gene.val <- read.csv("../../data/ex_12/gene_validation_body.csv", header=TRUE, row.names = 1)  #Gene validation set
df.prot.val <- read.csv("../../data/ex_12/prot_validation_body.csv", header=TRUE, row.names = 1)  #Gene validation set
df.gp <- cbind(df.prot.val, df.gene.val)
df.meta <- read.csv("../../data/ex_12/gp_validation_meta.csv", header=TRUE, row.names = 1)
df.meta$group <- as.character(df.meta$group)

# Create meta variant for easy tb vs non-tb identification w/in hiv status group
df.meta.tb_nontb <- df.meta
df.meta.tb_nontb$group <- as.numeric(df.meta.tb_nontb$group)
df.meta.tb_nontb$group[df.meta.tb_nontb$group == 3] <- 7
df.meta.tb_nontb$group[df.meta.tb_nontb$group == 6] <- 7
df.meta.tb_nontb$group[df.meta.tb_nontb$group == 4] <- 8
df.meta.tb_nontb$group[df.meta.tb_nontb$group == 5] <- 8
df.meta.tb_nontb$group <- as.character(df.meta.tb_nontb$group)

# Create meta variant for easy tb vs other factor identification between hiv status groups
df.meta.allhiv <- df.meta
df.meta.allhiv$group <- as.numeric(df.meta.allhiv$group)
df.meta.allhiv$group[df.meta.allhiv$group == 1] <- 9
df.meta.allhiv$group[df.meta.allhiv$group == 2] <- 9
df.meta.allhiv$group[df.meta.allhiv$group == 3] <- 10
df.meta.allhiv$group[df.meta.allhiv$group == 4] <- 10
df.meta.allhiv$group[df.meta.allhiv$group == 5] <- 11
df.meta.allhiv$group[df.meta.allhiv$group == 6] <- 11
df.meta.allhiv$group <- as.character(df.meta.allhiv$group)

# Create meta variant for easy tb vs non-tb identification between hiv status groups
df.meta.allhiv.tb_nontb <- df.meta
df.meta.allhiv.tb_nontb$group <- as.numeric(df.meta.allhiv.tb_nontb$group)
df.meta.allhiv.tb_nontb$group[df.meta.allhiv.tb_nontb$group == 1] <- 9
df.meta.allhiv.tb_nontb$group[df.meta.allhiv.tb_nontb$group == 2] <- 9
df.meta.allhiv.tb_nontb$group[df.meta.allhiv.tb_nontb$group == 3] <- 12
df.meta.allhiv.tb_nontb$group[df.meta.allhiv.tb_nontb$group == 4] <- 12
df.meta.allhiv.tb_nontb$group[df.meta.allhiv.tb_nontb$group == 5] <- 12
df.meta.allhiv.tb_nontb$group[df.meta.allhiv.tb_nontb$group == 6] <- 12
df.meta.allhiv.tb_nontb$group <- as.character(df.meta.allhiv.tb_nontb$group)

comps <- list(
  list(sel.gene.tb_od, sel.prot.tb_od, sel.gp.tb_od, df.meta, 1, 6, "TB", "OD", "tb_od")
  ,
  list(sel.gene.tb_ltbi, sel.prot.tb_ltbi, sel.gp.tb_ltbi, df.meta, 1, 3, "TB", "LTBI", "tb_ltbi")
  ,
  list(sel.gene.tb_nontb, sel.prot.tb_nontb, sel.gp.tb_nontb, df.meta.tb_nontb, 1, 7, "TB", "non-TB", "tb_nontb")
) 

gene_rocs <- list()

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
  df.comp.gp <- df.gp[ind.comp,]
  df.comp.meta <- comp.meta[ind.comp,]
  
  #Individual protein ROCs
  
  prot.aucs <- c()
  prot.down.cis <- c()
  prot.up.cis <- c()
  
  for (i in 1:length(comp.prot$features)){
    feat <- df.comp.prot[,match(comp.prot$features[i], colnames(df.comp.prot))]
    
    feat.roc <- roc(df.comp.meta$group, feat, auc=TRUE, ci=TRUE)
    prot.aucs <- c(prot.aucs, feat.roc$auc)
    prot.down.cis <- c(prot.down.cis, feat.roc$ci[1])
    prot.up.cis <- c(prot.up.cis, feat.roc$ci[3])
    
  }
  
  df.aucs.prots <- cbind(comp.prot[-5], data.frame(prot.aucs), data.frame(prot.down.cis), data.frame(prot.up.cis))
  colnames(df.aucs.prots) <- c("features", "coefficients", "beta", "reg_dir", "aucs", "lower.ci", "upper.ci")
  
  write.csv(df.aucs.prots, paste("../../data/", ex_dir, "feat_sel_2/", "prot_", comp.abbrv, "_BH_LFC_lasso_sig_factors_withaucs.csv", sep=""), row.names=TRUE)
  
  # Get individual aucs for gene features
  
  gene.aucs <- c()
  gene.down.cis <- c()
  gene.up.cis <- c()
  
  for (i in 1:length(comp.gene$features)){
    feat <- df.comp.gene[,match(sel.gene.tb_nontb$features[i], colnames(df.comp.gene))]
    
    feat.roc <- roc(df.comp.meta$group, feat, auc=TRUE, ci=TRUE)
    gene.aucs <- c(gene.aucs, feat.roc$auc)
    gene.down.cis <- c(gene.down.cis, feat.roc$ci[1])
    gene.up.cis <- c(gene.up.cis, feat.roc$ci[3])
  }
  
  df.aucs.genes <-  cbind(comp.gene, data.frame(gene.aucs), data.frame(gene.down.cis), data.frame(gene.up.cis))
  colnames(df.aucs.genes) <- c("features", "coefficients", "beta", "reg_dir", "aucs", "lower.ci", "upper.ci")
  
  write.csv(df.aucs.genes, paste("../../data/", ex_dir, "feat_sel_2/", "gene_", comp.abbrv, "_BH_LFC_lasso_sig_factors_withaucs.csv", sep=""), row.names=TRUE)
  
  df.aucs.gp <- rbind(df.aucs.prots, df.aucs.genes)
  df.aucs.gp <- df.aucs.gp[order(-df.aucs.gp$aucs),]
  
  write.csv(df.aucs.gp, paste("../../data/", ex_dir, "feat_sel_2/", "gp_", comp.abbrv, "_BH_LFC_lasso_sig_factors_withaucs.csv", sep=""), row.names=TRUE)
 
} 
  
  # Seperately compare protein and  gene, then both combined by DRS 
 # sets <- list(
 #   list(df.comp.gene, comp.gene, "gene", "gene"),
 #   list(df.comp.prot, comp.prot, "protein", "prot"),
 #   list(df.comp.gp, comp.all, "combined gene and protein", "gp")
 # )
  
  # Set empty variables to become ROC objects for comparison
 # roc.gene <- 0
 # roc.prot <- 0
 # roc.gp <- 0
 # roc.meta <- 0
  
#for (set in sets){
   # set.data <- set[[1]]
   # set.sel_features <- set[[2]]
   # set.verbose <- set[[3]][1]
    #set.abbrv <- set[[4]][1]
    
    ####################################################################################
    ## Get DRS ROC curves - First for protein and gene seperately, then both combined ##
    ####################################################################################
    
   # upreg.features <- c()
   # downreg.features <- c()
    
   # for (i in 1:nrow(set.sel_features)){
     # if (set.sel_features$reg_dir[i] == "up"){
     #   upreg.features <- c(upreg.features, as.character(set.sel_features$features[i]))
        #print(paste("UP:", comp.sel.features$features[i]))
     # } else {
       # downreg.features <- c(downreg.features, as.character(set.sel_features$features[i]))
        #print(paste("DOWN:", comp.sel.features$features[i]))
    #  }
   # }
    
    #df.upreg <- set.data[,match(upreg.features, colnames(set.data))]
    #df.downreg <- set.data[,match(downreg.features, colnames(set.data))]
    
    # Get DRS values of patients
    
    #drs <- c()
    
   # if (length(upreg.features) == 0){
    #  if (length(downreg.features) == 1){
        
    #    for (i in 1:nrow(set.data)){
    #      drs.i <-  - sum(df.downreg[i])
    #      drs <- c(drs, drs.i)
   #     }
        
   #   } 
   # } else if (length(upreg.features) == 1){
    #  if (length(downreg.features) == 1){
        
    #    for (i in 1:nrow(set.data)){
    #      drs.i <- df.upreg[i] - df.downreg[i]
     #     drs <- c(drs, drs.i)
    #    }
        
    #  } else if (length(downreg.features) == 0){
        
   #     for (i in 1:nrow(set.data)){
   #       drs.i <- sum(df.upreg[i])
   #       drs <- c(drs, drs.i)
   #     }
   #     
   #   }
   # } else if (length(upreg.features) > 1) {
      
    #  if (length(downreg.features) == 1){
        
    #    for (i in 1:nrow(set.data)){
     #     drs.i <- sum(df.upreg[i,]) - df.downreg[i]
     #     drs <- c(drs, drs.i)
     #   }

   #   } else {
      
   #     for (i in 1:nrow(set.data)){
    #      drs.i <- sum(df.upreg[i,]) - sum(df.downreg[i,])
   #       drs <- c(drs, drs.i)
  #      }
  #    }
   # }
  
    # Select DRS by TB status
    
  #  ind.g1 <- c()
 #   ind.g2 <- c()
 #   ind.g1_g2 <- c()
    
    
 #   for (j in 1:nrow(df.comp.meta)){
 #     if (df.comp.meta$group[j] == comp.g1){
 #       ind.g1 <- c(ind.g1, j)
  #    }
  #  }
    
   # for (j in 1:nrow(df.comp.meta)){
   #   if (df.comp.meta$group[j] == comp.g2){
   #     ind.g2 <- c(ind.g2, j)
   #   }
  #  }
    
   # for (j in 1:nrow(df.comp.meta)){
  #    if ((df.comp.meta$group[j] == comp.g1) || (df.comp.meta$group[j] == comp.g2)){
 #       ind.g1_g2 <- c(ind.g1_g2, j)
  #    }
  #  }
    
   # comp.meta.g1_g2 <- df.comp.meta[ind.g1_g2,]
    
  #  drs.g1 <- drs[ind.g1]
  #  drs.g2 <- drs[ind.g2]
 #   drs.g1_g2 <- drs[ind.g1_g2]
    
    #Get means and sds of drs
    
 #   mean.drs.g1 <- mean(drs.g1)
 #   mean.drs.g2 <- mean(drs.g2)
  #  sd.drs.g1 <- sd(drs.g1)
  #  sd.drs.g2 <- sd(drs.g2)
    
    # Get threshold and get boxplot
 #   if (set.abbrv == "gp"){
      
  #    melt.gp <- data.frame(comp.meta.g1_g2$group, drs.g1_g2)
  #    colnames(melt.gp) <- c("group", "gp.drs")
  #    drs.threshold.gp <- ((mean.drs.g1/sd.drs.g1)+(mean.drs.g2/sd.drs.g2))/((1/sd.drs.g1)+(1/sd.drs.g2))
      
  #    bplot.gp <- ggplot(data=melt.gp, aes(x=group, y=gp.drs)) +
   #     geom_boxplot(fill="purple", alpha=0.8) +
   #     geom_jitter(aes(x=group,y=gp.drs),
     #                         position=position_jitter(width=0.1,height=0),
      #                        alpha=0.6,
      #                        size=3) +
      #  geom_hline(yintercept = drs.threshold.gp) +
      #  xlab("Group") +
      #  scale_x_discrete(labels=c(comp.verb.g1, comp.verb.g2)) +
      #  ylab("DRS") +
      #  theme_minimal() +
      #  theme(axis.text.x  = element_text(size = 10),
      #        axis.text.y = element_text(size = 10),
     #         axis.title.x = element_text(size = 20),
     #         axis.title.y = element_text(size = 20)
      #  )
      
     # ggsave(paste("../../img/", ex_dir, date, "val/", set.abbrv, "_", comp.abbrv, "_val_drs_thresh_boxplot.png", sep=""),
    #         bplot.gp,
    #         width=5, 
    #         height=15, 
    #         unit="cm")
      
    #} if (set.abbrv == "prot"){
      
   #   drs.threshold.prot <- ((mean.drs.g1/sd.drs.g1)+(mean.drs.g2/sd.drs.g2))/((1/sd.drs.g1)+(1/sd.drs.g2))
      
  #    melt.prot <- data.frame(comp.meta.g1_g2$group, drs.g1_g2)
   #   colnames(melt.prot) <- c("group", "prot.drs")
   #   drs.threshold.prot <- ((mean.drs.g1/sd.drs.g1)+(mean.drs.g2/sd.drs.g2))/((1/sd.drs.g1)+(1/sd.drs.g2))
      
   #   bplot.prot <- ggplot(data=melt.prot, aes(x=group, y=prot.drs)) +
   #     geom_boxplot(fill="red", alpha=0.8) +
   #     geom_jitter(aes(x=group,y=prot.drs),
    #                position=position_jitter(width=0.1,height=0),
    #                alpha=0.6,
    #                size=3) +
    #    geom_hline(yintercept = drs.threshold.prot) +
    #    xlab("Group") +
    #    scale_x_discrete(labels=c(comp.verb.g1, comp.verb.g2)) +
    #    ylab("DRS") +
   #     theme_minimal() +
   #     theme(axis.text.x  = element_text(size = 10),
   #           axis.text.y = element_text(size = 10),
   #           axis.title.x = element_text(size = 20),
   #           axis.title.y = element_text(size = 20)
    #    )
      
   #   ggsave(paste("../../img/", ex_dir, date, "val/", set.abbrv, "_", comp.abbrv, "_val_drs_thresh_boxplot.png", sep=""),
  #           bplot.prot,
  #           width=5, 
  #           height=15, 
  #           unit="cm")
      
   # } if (set.abbrv == "gene"){
      
   #   drs.threshold.gene <- ((mean.drs.g1/sd.drs.g1)+(mean.drs.g2/sd.drs.g2))/((1/sd.drs.g1)+(1/sd.drs.g2))
      
   #   melt.gene <- data.frame(comp.meta.g1_g2$group, drs.g1_g2)
   #   colnames(melt.gene) <- c("group", "gene.drs")
   #   drs.threshold.gene <- ((mean.drs.g1/sd.drs.g1)+(mean.drs.g2/sd.drs.g2))/((1/sd.drs.g1)+(1/sd.drs.g2))
      
   #   bplot.gene <- ggplot(data=melt.gene, aes(x=group, y=gene.drs)) +
     #   geom_boxplot(fill="blue", alpha=0.8) +
       # geom_jitter(aes(x=group,y=gene.drs),
       #             position=position_jitter(width=0.1,height=0),
       #             alpha=0.6,
       #             size=3) +
       # geom_hline(yintercept = drs.threshold.gene) +
       # xlab("Group") +
       # scale_x_discrete(labels=c(comp.verb.g1, comp.verb.g2)) +
       # ylab("DRS") +
        #theme_minimal() +
        #theme(axis.text.x  = element_text(size = 10),
        #      axis.text.y = element_text(size = 10),
        #      axis.title.x = element_text(size = 20),
        #      axis.title.y = element_text(size = 20)
        #)
      
     # ggsave(paste("../../img/", ex_dir, date, "val/", set.abbrv, "_", comp.abbrv, "_val_drs_thresh_boxplot.png", sep=""),
     #        bplot.gene,
     #        width=5, 
     #        height=15, 
    #         unit="cm")
    #  
    #}
    # Get ROC
    
    #drs.roc <- roc(comp.meta.g1_g2$group, drs.g1_g2, auc=TRUE)
    
    #roc.groups <- comp.meta.g1_g2$group
    
    #if (set.abbrv == "prot"){
    ##  roc.prot <- drs.roc
    #} else if (set.abbrv == "gene"){
    #  roc.gene <- drs.roc
    #} else if (set.abbrv == "gp"){
    #  roc.gp <- drs.roc
    #}
    
    #if (set.verbose == "gene"){
    #  my_rocs_list <- list(my_rocs_list, list(drs.roc, comp.abbrv))
    #}
    
    #if (drs.roc$auc < 1){
      
      #drs.roc.smooth <- smooth(drs.roc)
      
      #png(paste("../../img/", ex_dir, date, "val/", set.abbrv, "_", comp.abbrv, "_val_drs_roc.png", sep=""),
       #   width = 5*300,        # 5 x 300 pixels
      #    height = 5*300,
      #    res = 300,            # 300 pixels per inch
      #    pointsize = 8) # smaller font size
      
      #plot(drs.roc,
     #      col="red",
      #     lwd=1,
      #     main=paste("ROC curve for selected features for", set.verbose, "data,", "HIV-patients," , comp.verb.g1,  "vs", comp.verb.g2),
      #     legacy.axes = TRUE
      #)
      
      #lines.roc(drs.roc.smooth, col="blue", lwd=3)
     # legend("bottomright",
     #        title = paste("AUC =", format(round(drs.roc$auc, 2), nsmall = 2)),
     #        #legend = c("Empirical data", "Fit"),
     #        #col = c("red", "blue"),
     #        #lwd = c(1, 3)
     #        legend = c("Empirical data"),
     #        col = c("red"),
      #       lwd = c(1)
      #)
      
      #dev.off()
      
    #} else {
      
      #png(paste("../../img/", ex_dir, date, "val/", set.abbrv, "_", comp.abbrv, "_val_drs_roc.png", sep=""),
      #    width = 5*300,        # 5 x 300 pixels
      #    height = 5*300,
     #     res = 300,            # 300 pixels per inch
      #    pointsize = 8) # smaller font size
      
     # plot(drs.roc,
      #     col="red",
     #      lwd=1,
     #      main=paste("ROC curve for selected features for", set.verbose, "data,", "HIV-patients," , comp.verb.g1,  "vs", comp.verb.g2),
     #      legacy.axes = TRUE
      #)
      
      #legend("bottomright",
       #      title = paste("AUC =", format(round(drs.roc$auc, 2), nsmall = 2)),
      #       legend = c("Empirical data"),
      #       col = c("red"),
     #        lwd = c(1, 3)
     #)
      
      #dev.off()
      
    #}
    
    # Create dataframe and ggplot boxplot
    #inf.group <- as.character(comp.meta.g1_g2$group)
    
    #data = data.frame(inf.group, drs.g1_g2)
    
    #ggplot(data, aes(x=inf.group, y=drs.g1_g2, group=inf.group)) + 
      #geom_boxplot() +
      #labs(x = "TB status", y="Disease risk score", title =paste("Distribution of DRS values for", set.verbose, "data,", "HIV-patients," , comp.verb.g1,  "vs", comp.verb.g2)) +
      #scale_x_discrete(labels=c(comp.verb.g1, comp.verb.g2))
    
    #ggsave(paste("../../img/", ex_dir, date, "val/", set.abbrv, "_", comp.abbrv, "_val_drs_boxplot.png", sep=""))
     
 # }
  
  #gene_rocs <- list(gene_rocs, roc.gene)
  
  # Plot superimposed ROCs
  #png(paste("../../img/", ex_dir, date, "val/", comp.abbrv, "_val_superimposed_drs_rocs.png", sep=""),
  #    width = 5*300,        # 5 x 300 pixels
  #    height = 5*300,
  #    res = 300,            # 300 pixels per inch
  #    pointsize = 8) # smaller font size
  
  #plot(roc.gp,
  #     col = alpha("purple", 0.7),
  #     lwd=4,
  #     #main=paste("Superimposed ROC curve for selected features for seperate and combined data types, HIV-" , comp.verb.g1,  "vs", comp.verb.g2),
  #     legacy.axes = TRUE,
  #     cex.lab=3,
  #     xlab="",
  #     ylab=""
  #     #par(mai=c(10,10,10,10))
  #)
  
  #par(mar=c(10,10,10,10))
  
  #title(main="", xlab="1 - Specificity", ylab="Sensitivity", cex.lab=3, pos=1)
  
  #lines.roc(roc.gene,
  #          col = alpha("blue", 0.7)
  #)
#
  #lines.roc(roc.prot,
  #          col = alpha("red", 0.7)
  #)
  
  #legend("bottomright",
  #       #title = paste("AUC =", format(round(drs.roc$auc, 2), nsmall = 2)),
  #       legend = c("Gene", "Protein", "Combined"),
  #       col = c("blue", "red", "purple"),
  #       lwd = c(2, 2, 4)
  #)
  
  #dev.off()
  
#}

# Run roc11_drs_combo first to get kaforou rocs
#test_rocs_tb_od <- roc.test(gene_rocs[[1]][[1]][[2]], kaforou_rocs_list[[1]][[1]][[2]][[1]])
#test_rocs_tb_ltbi <- roc.test(gene_rocs[[1]][[2]], kaforou_rocs_list[[1]][[2]][[1]])
#test_rocs_tb_nontb <- roc.test(gene_rocs[[2]], kaforou_rocs_list[[2]][[1]])

#test_rocs_tb_od
#test_rocs_tb_ltbi
#test_rocs_tb_nontb
