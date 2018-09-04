#setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_12/")
setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/ex_9/")

library(pROC)
library(ggplot2)
library(stringr)

set.seed(12)

# To direct to the correct folder
date <- "2018-09-03/"
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

df.gene.train <- read.csv("../../data/ex_12/gene_train_body.csv", header=TRUE, row.names = 1)  #Gene validation set
df.prot.train <- read.csv("../../data/ex_12/prot_train_body.csv", header=TRUE, row.names = 1)  #Protein validation set
df.gp.train <- cbind(df.prot.train, df.gene.train)
df.meta.train <- read.csv("../../data/ex_12/gp_train_meta.csv", header=TRUE, row.names = 1)
df.meta.train$group <- as.character(df.meta.train$group)

df.gene.val <- read.csv("../../data/ex_12/gene_validation_body.csv", header=TRUE, row.names = 1)  #Gene validation set
df.prot.val <- read.csv("../../data/ex_12/prot_validation_body.csv", header=TRUE, row.names = 1)  #Protein validation set
df.gp <- cbind(df.prot.val, df.gene.val)
df.meta <- read.csv("../../data/ex_12/gp_validation_meta.csv", header=TRUE, row.names = 1)
df.meta$group <- as.character(df.meta$group)

testrains <- list(
  list(),
  list()
)

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
  
  # Seperately compare protein and  gene, then both combined by DRS 
  sets <- list(
    list(df.comp.gene, comp.gene, "gene", "gene"),
    list(df.comp.prot, comp.prot, "protein", "prot"),
    list(df.comp.gp, comp.all, "combined gene and protein", "gp")
  )
  
  # Set empty variables to become ROC objects for comparison
  roc.gene <- 0
  roc.prot <- 0
  roc.gp <- 0
  roc.meta <- 0
  
  #Get melt tables for combined boxplot
  melt.gene <- 0
  melt.prot <- 0
  melt.gp <- 0
  drs.threshold <- 0
    
  
  
  for (set in sets){
    set.data <- set[[1]]
    set.sel_features <- set[[2]]
    set.verbose <- set[[3]][1]
    set.abbrv <- set[[4]][1]
    
    ####################################################################################
    ## Get DRS ROC curves - First for protein and gene seperately, then both combined ##
    ####################################################################################
    
    upreg.features <- c()
    downreg.features <- c()
    
    for (i in 1:nrow(set.sel_features)){
      if (set.sel_features$reg_dir[i] == "up"){
        upreg.features <- c(upreg.features, as.character(set.sel_features$features[i]))
        #print(paste("UP:", comp.sel.features$features[i]))
      } else {
        downreg.features <- c(downreg.features, as.character(set.sel_features$features[i]))
        #print(paste("DOWN:", comp.sel.features$features[i]))
      }
    }
    
    df.upreg <- set.data[,match(upreg.features, colnames(set.data))]
    df.downreg <- set.data[,match(downreg.features, colnames(set.data))]
    
    # Get DRS values of patients
    
    drs <- c()
    
    if (length(upreg.features) == 0){
      if (length(downreg.features) == 1){
        
        for (i in 1:nrow(set.data)){
          drs.i <-  - sum(df.downreg[i])
          drs <- c(drs, drs.i)
        }
        
      } 
    } else if (length(upreg.features) == 1){
      if (length(downreg.features) == 1){
        
        for (i in 1:nrow(set.data)){
          drs.i <- df.upreg[i] - df.downreg[i]
          drs <- c(drs, drs.i)
        }
        
      } else if (length(downreg.features) == 0){
        
        for (i in 1:nrow(set.data)){
          drs.i <- sum(df.upreg[i])
          drs <- c(drs, drs.i)
        }
        
      }
    } else if (length(upreg.features) > 1) {
      
      if (length(downreg.features) == 1){
        
        for (i in 1:nrow(set.data)){
          drs.i <- sum(df.upreg[i,]) - df.downreg[i]
          drs <- c(drs, drs.i)
        }

      } else {
      
        for (i in 1:nrow(set.data)){
          drs.i <- sum(df.upreg[i,]) - sum(df.downreg[i,])
          drs <- c(drs, drs.i)
        }
      }
    }
  
    # Select DRS by TB status
    
    ind.g1 <- c()
    ind.g2 <- c()
    ind.g1_g2 <- c()
    
    
    for (j in 1:nrow(df.comp.meta)){
      if (df.comp.meta$group[j] == comp.g1){
        ind.g1 <- c(ind.g1, j)
      }
    }
    
    for (j in 1:nrow(df.comp.meta)){
      if (df.comp.meta$group[j] == comp.g2){
        ind.g2 <- c(ind.g2, j)
      }
    }
    
    for (j in 1:nrow(df.comp.meta)){
      if ((df.comp.meta$group[j] == comp.g1) || (df.comp.meta$group[j] == comp.g2)){
        ind.g1_g2 <- c(ind.g1_g2, j)
      }
    }
    
    comp.meta.g1_g2 <- df.comp.meta[ind.g1_g2,]
    
    drs.g1 <- drs[ind.g1]
    drs.g2 <- drs[ind.g2]
    drs.g1_g2 <- drs[ind.g1_g2]
    
    #Get means and sds of drs
    
    mean.drs.g1 <- mean(drs.g1)
    mean.drs.g2 <- mean(drs.g2)
    sd.drs.g1 <- sd(drs.g1)
    sd.drs.g2 <- sd(drs.g2)
    
    # Get threshold and get boxplot
    if (set.abbrv == "gp"){
      
      melt.gp <- data.frame(comp.meta.g1_g2$group, drs.g1_g2)
      colnames(melt.gp) <- c("group", "gp.drs")
      drs.threshold.gp <- ((mean.drs.g1/sd.drs.g1)+(mean.drs.g2/sd.drs.g2))/((1/sd.drs.g1)+(1/sd.drs.g2))
      
      bplot.gp <- ggplot(data=melt.gp, aes(x=group, y=gp.drs)) +
        geom_boxplot(fill="#FF0447", alpha=0.8) +
        geom_jitter(aes(x=group,y=gp.drs),
                    position=position_jitter(width=0.1,height=0),
                    alpha=0.6,
                    size=3) +
        geom_hline(yintercept = drs.threshold.gp) +
        xlab("Group") +
        scale_x_discrete(labels=c(comp.verb.g1, comp.verb.g2)) +
        ylab("DRS") +
        theme_minimal() +
        theme(axis.text.x  = element_text(size = 10),
              axis.text.y = element_text(size = 10),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size = 20)
        )
      
      ggsave(paste("../../img/", ex_dir, date, "val/", set.abbrv, "_", comp.abbrv, "_val_drs_thresh_boxplot.png", sep=""),
             bplot.gp,
             width=5, 
             height=15, 
             unit="cm")
      
    } else if (set.abbrv == "prot"){
      
      drs.threshold.prot <- ((mean.drs.g1/sd.drs.g1)+(mean.drs.g2/sd.drs.g2))/((1/sd.drs.g1)+(1/sd.drs.g2))
      
      melt.prot <- data.frame(comp.meta.g1_g2$group, drs.g1_g2)
      colnames(melt.prot) <- c("group", "prot.drs")
      drs.threshold.prot <- ((mean.drs.g1/sd.drs.g1)+(mean.drs.g2/sd.drs.g2))/((1/sd.drs.g1)+(1/sd.drs.g2))
      
      bplot.prot <- ggplot(data=melt.prot, aes(x=group, y=prot.drs)) +
        geom_boxplot(fill="#FFDE04", alpha=0.8) +
        geom_jitter(aes(x=group,y=prot.drs),
                    position=position_jitter(width=0.1,height=0),
                    alpha=0.6,
                    size=3) +
        geom_hline(yintercept = drs.threshold.prot) +
        xlab("Group") +
        scale_x_discrete(labels=c(comp.verb.g1, comp.verb.g2)) +
        ylab("DRS") +
        theme_minimal() +
        theme(axis.text.x  = element_text(size = 10),
              axis.text.y = element_text(size = 10),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size = 20)
        )
      
      ggsave(paste("../../img/", ex_dir, date, "val/", set.abbrv, "_", comp.abbrv, "_val_drs_thresh_boxplot.png", sep=""),
             bplot.prot,
             width=5, 
             height=15, 
             unit="cm")
      
    } else if (set.abbrv == "gene"){
      
      drs.threshold.gene <- ((mean.drs.g1/sd.drs.g1)+(mean.drs.g2/sd.drs.g2))/((1/sd.drs.g1)+(1/sd.drs.g2))
      
      melt.gene <- data.frame(comp.meta.g1_g2$group, drs.g1_g2)
      colnames(melt.gene) <- c("group", "gene.drs")
      drs.threshold.gene <- ((mean.drs.g1/sd.drs.g1)+(mean.drs.g2/sd.drs.g2))/((1/sd.drs.g1)+(1/sd.drs.g2))
      
      bplot.gene <- ggplot(data=melt.gene, aes(x=group, y=gene.drs)) +
        geom_boxplot(fill="#048AFF", alpha=0.8) +
        geom_jitter(aes(x=group,y=gene.drs),
                    position=position_jitter(width=0.1,height=0),
                    alpha=0.6,
                    size=3) +
        geom_hline(yintercept = drs.threshold.gene) +
        xlab("Group") +
        scale_x_discrete(labels=c(comp.verb.g1, comp.verb.g2)) +
        ylab("DRS") +
        theme_minimal() +
        theme(axis.text.x  = element_text(size = 10),
              axis.text.y = element_text(size = 10),
              axis.title.x = element_blank(),
              axis.title.y = element_text(size = 20)
        )
      
      ggsave(paste("../../img/", ex_dir, date, "val/", set.abbrv, "_", comp.abbrv, "_val_drs_thresh_boxplot.png", sep=""),
             bplot.gene,
             width=5, 
             height=15, 
             unit="cm")
      
    }
    
    # Get ROC
    
    drs.roc <- roc(comp.meta.g1_g2$group, drs.g1_g2, auc=TRUE, ci=TRUE)
    
    roc.groups <- comp.meta.g1_g2$group
    
    if (set.abbrv == "prot"){
      roc.prot <- drs.roc
      
     
    } else if (set.abbrv == "gene"){
      roc.gene <- drs.roc
      
     
      
    } else if (set.abbrv == "gp"){
      roc.gp <- drs.roc
      
      
      
     
    }
  }
  
  gene_rocs <- list(gene_rocs, roc.gene)
  
  #Display spec, sens, and ci of roc
  print(comp.abbrv)
  print(roc.gene$thresholds)
  print(roc.gene$sensitivities)
  print(roc.gene$specificities)
  print(roc.gene$ci)
  print(ci.thresholds(roc.gene))
  print(roc.prot$thresholds)
  print(roc.prot$sensitivities)
  print(roc.prot$specificities)
  print(roc.prot$ci)
  print(ci.thresholds(roc.prot))
  print(roc.gp$thresholds)
  print(roc.gp$sensitivities)
  print(roc.gp$specificities)
  print(roc.gp$ci)
  print(ci.thresholds(roc.gp))
  
  
  # Plot superimposed ROCs
  png(paste("../../img/", ex_dir, date, "val/", comp.abbrv, "_val_superimposed_drs_rocs.png", sep=""),
      width = 3000,        # 5 x 300 pixels
      height = 3000,
      res = 300,            # 300 pixels per inch
      pointsize = 8) # smaller font size
  
  plot(roc.gp,
       col = alpha("#FF0447", 0.7),
       lwd=32,
       #main=paste("Superimposed ROC curve for selected features for seperate and combined data types, HIV-" , comp.verb.g1,  "vs", comp.verb.g2),
       legacy.axes = TRUE,
       #cex.lab=4,
       cex.axis=2,
       xlab="",
       ylab="",
       asp=NA
  )
  
  #axis(side=1, cex.axis=1.5)
  #axis(side=2, cex.axis=1.5)
  
  par(mar=c(5,5,5,5))
  
  #title(main="", xlab="1 - Specificity", ylab="Sensitivity", cex.lab=3, pos=1)
  
  lines.roc(roc.gene,
            col = alpha("#048AFF", 1),
            lwd=8
  )

  lines.roc(roc.prot,
            col = alpha("#FFDE04", 1),
            lwd=8
  )
  
  legend("bottomright",
         legend = c(paste("Gene", "(AUC =", format(round(roc.gene$auc, 2), nsmall = 2), ")", sep=" "),
                    paste("Protein", "(AUC =", format(round(roc.prot$auc, 2), nsmall = 2), ")", sep=" "),
                    paste("Combined", "(AUC =", format(round(roc.gp$auc, 2), nsmall = 2), ")", sep=" ")
                    ),
         col = c("#048AFF", "#FFDE04", "#FF0447"),
         lwd = c(8, 8, 16),
         cex=3
  )
  
  
  
  dev.off()
  
  # Plot side-by side boxplots w/threshold
  
  #melt.drs <- rbind(melt.prot, melt.gene, melt.gp)
  
  #ggplot(data=melt.gp, aes(x=tb_stat, y=drs.value)) +
    #geom_boxplot(aes(fill=omtype)) +
   # geom_jitter(aes(x=tb_stat,y=drs.value),
      #          position=position_jitter(width=0.1,height=0),
   #             alpha=0.6,
      #          size=3) +
    #geom_hline(yintercept=drs.threshold)
  
  
}

# Run roc11_drs_combo first to get kaforou rocs
test_rocs_tb_od <- roc.test(gene_rocs[[1]][[1]][[2]], kaforou_rocs_list[[1]][[1]][[2]][[1]])
test_rocs_tb_ltbi <- roc.test(gene_rocs[[1]][[2]], kaforou_rocs_list[[1]][[2]][[1]])
test_rocs_tb_nontb <- roc.test(gene_rocs[[2]], kaforou_rocs_list[[2]][[1]])

test_rocs_tb_od
test_rocs_tb_ltbi
test_rocs_tb_nontb
