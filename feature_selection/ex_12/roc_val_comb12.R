setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_12/")
#setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/ex_9/")

library(pROC)
library(ggplot2)
library(stringr)

set.seed(12)

# To direct to the correct folder
date <- "2018-08-24/"
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
  
  for (i in 1:length(comp.prot$features)){
    feat <- df.comp.prot[,match(comp.prot$features[i], colnames(df.comp.prot))]
    
    feat.roc <- roc(df.comp.meta$group, feat, auc=TRUE)
    
    # Not all ROCs smoothable; deal with them
    if ((feat.roc$auc == 1) || (feat.roc$auc == -1)){
      
      pdf(paste("../../img/", ex_dir, date, "val/", comp.abbrv, "_prot_val_indv_roc", comp.prot$gr_features[i], ".pdf", sep=""),
          width = 5*300,        # 5 x 300 pixels
          height = 5*300,
          res = 300,            # 300 pixels per inch
          pointsize = 8) # smaller font size
      
      plot(feat.roc,
           col="red",
           lwd=1,
           #main=paste(comp.prot$features[i], "," , comp.verb.g1,  "vs", comp.verb.g2),
           legacy.axes = TRUE
      )
      
      legend("bottomright",
             title = paste("AUC =", format(round(feat.roc$auc, 2), nsmall = 2)),
             #legend = c("Empirical data", "Fit"),
             #col = c("red", "blue"),
             #lwd = c(1, 3)
             legend = c("Empirical data"),
             col = c("red"),
             lwd = c(1)
      )
      
    } else {
    
      #feat.roc.smooth <- smooth(feat.roc)
      
      pdf(paste("../../img/", ex_dir, date, "val/", comp.abbrv, "_prot_val_indv_roc", comp.prot$gr_features[i], ".pdf", sep=""),
          width = 5*300,        # 5 x 300 pixels
          height = 5*300,
          res = 300,            # 300 pixels per inch
          pointsize = 8) # smaller font size
      
      plot(feat.roc,
           col="red",
           lwd=1,
           #main=paste(comp.prot$features[i], "," , comp.verb.g1,  "vs", comp.verb.g2),
           legacy.axes = TRUE
      )
      
      #lines.roc(feat.roc.smooth, col="blue", lwd=3)
      legend("bottomright",
             title = paste("AUC =", format(round(feat.roc$auc, 2), nsmall = 2)),
             #legend = c("Empirical data", "Fit"),
             #col = c("red", "blue"),
             #lwd = c(1, 3)
             legend = c("Empirical data"),
             col = c("red"),
             lwd = c(1)
      )
      
      dev.off()
    }
    
    # Create dataframe and ggplot boxplot
    inf.group <- as.character(df.comp.meta$group)
    
    data = data.frame(inf.group, feat)
    
    ggplot(data, aes(x=inf.group, y=feat, group=inf.group)) + 
      geom_boxplot() +
      labs(x = "TB status",
           y = as.character(str_replace_all(comp.prot$gr_features[i], "_", " "))
           #title =paste("Distribution of DRS values for", comp.prot$features[i], "HIV-patients," , comp.verb.g1,  "vs", comp.verb.g2)
           ) +
      scale_x_discrete(labels=c(comp.verb.g1, comp.verb.g2))
    
    ggsave(paste("../../img/", ex_dir, date, "val/", comp.abbrv, "_prot_val_indv_boxplot", comp.prot$gr_features[i], ".pdf", sep=""))
    
  }
  
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
    
    # Get threshold
    if (set.abbrv == "gp"){
      drs.threshold <- ((mean.drs.g1/sd.drs.g1)+(mean.drs.g2/sd.drs.g2))/((1/sd.drs.g1)+(1/sd.drs.g2))
    
      print(paste("Threshold for ", comp.verb.g1, " vs ", comp.verb.g2, ": ", drs.threshold, sep=""))
    }
    # Get ROC
    
    drs.roc <- roc(comp.meta.g1_g2$group, drs.g1_g2, auc=TRUE)
    
    roc.groups <- comp.meta.g1_g2$group
    
    if (set.abbrv == "prot"){
      roc.prot <- drs.roc
    } else if (set.abbrv == "gene"){
      roc.gene <- drs.roc
    } else if (set.abbrv == "gp"){
      roc.gp <- drs.roc
    }
    
    #if (set.verbose == "gene"){
    #  my_rocs_list <- list(my_rocs_list, list(drs.roc, comp.abbrv))
    #}
    
    if (drs.roc$auc < 1){
      
      #drs.roc.smooth <- smooth(drs.roc)
      
      pdf(paste("../../img/", ex_dir, date, "val/", set.abbrv, "_", comp.abbrv, "_val_drs_roc.pdf", sep=""),
          width = 5*300,        # 5 x 300 pixels
          height = 5*300,
          res = 300,            # 300 pixels per inch
          pointsize = 8) # smaller font size
      
      plot(drs.roc,
           col="red",
           lwd=1,
           main=paste("ROC curve for selected features for", set.verbose, "data,", "HIV-patients," , comp.verb.g1,  "vs", comp.verb.g2),
           legacy.axes = TRUE
      )
      
      #lines.roc(drs.roc.smooth, col="blue", lwd=3)
      legend("bottomright",
             title = paste("AUC =", format(round(drs.roc$auc, 2), nsmall = 2)),
             #legend = c("Empirical data", "Fit"),
             #col = c("red", "blue"),
             #lwd = c(1, 3)
             legend = c("Empirical data"),
             col = c("red"),
             lwd = c(1)
      )
      
      dev.off()
      
    } else {
      
      pdf(paste("../../img/", ex_dir, date, "val/", set.abbrv, "_", comp.abbrv, "_val_drs_roc.pdf", sep=""),
          width = 5*300,        # 5 x 300 pixels
          height = 5*300,
          res = 300,            # 300 pixels per inch
          pointsize = 8) # smaller font size
      
      plot(drs.roc,
           col="red",
           lwd=1,
           main=paste("ROC curve for selected features for", set.verbose, "data,", "HIV-patients," , comp.verb.g1,  "vs", comp.verb.g2),
           legacy.axes = TRUE
      )
      
      legend("bottomright",
             title = paste("AUC =", format(round(drs.roc$auc, 2), nsmall = 2)),
             legend = c("Empirical data"),
             col = c("red"),
             lwd = c(1, 3)
      )
      
      dev.off()
      
    }
    
    # Create dataframe and ggplot boxplot
    inf.group <- as.character(comp.meta.g1_g2$group)
    
    data = data.frame(inf.group, drs.g1_g2)
    
    ggplot(data, aes(x=inf.group, y=drs.g1_g2, group=inf.group)) + 
      geom_boxplot() +
      labs(x = "TB status", y="Disease risk score", title =paste("Distribution of DRS values for", set.verbose, "data,", "HIV-patients," , comp.verb.g1,  "vs", comp.verb.g2)) +
      scale_x_discrete(labels=c(comp.verb.g1, comp.verb.g2))
    
    ggsave(paste("../../img/", ex_dir, date, "val/", set.abbrv, "_", comp.abbrv, "_val_drs_boxplot.pdf", sep=""))
     
  }
  
  gene_rocs <- list(gene_rocs, roc.gene)
  
  # Plot superimposed ROCs
  pdf(paste("../../img/", ex_dir, date, "val/", comp.abbrv, "_val_superimposed_drs_rocs.pdf", sep=""),
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8) # smaller font size
  
  plot(roc.gp,
       col = alpha("purple", 0.7),
       lwd=4,
       #main=paste("Superimposed ROC curve for selected features for seperate and combined data types, HIV-" , comp.verb.g1,  "vs", comp.verb.g2),
       legacy.axes = TRUE,
       cex.lab=3,
       xlab="",
       ylab=""
       #par(mai=c(10,10,10,10))
  )
  
  par(mar=c(10,10,10,10))
  
  #title(main="", xlab="1 - Specificity", ylab="Sensitivity", cex.lab=3, pos=1)
  
  lines.roc(roc.gene,
            col = alpha("blue", 0.7)
  )

  lines.roc(roc.prot,
            col = alpha("red", 0.7)
  )
  
  legend("bottomright",
         #title = paste("AUC =", format(round(drs.roc$auc, 2), nsmall = 2)),
         legend = c("Gene", "Protein", "Combined"),
         col = c("blue", "red", "purple"),
         lwd = c(2, 2, 4)
  )
  
  dev.off()
  
}

# Run roc11_drs_combo first to get kaforou rocs
test_rocs_tb_od <- roc.test(gene_rocs[[1]][[1]][[2]], kaforou_rocs_list[[1]][[1]][[2]][[1]])
test_rocs_tb_ltbi <- roc.test(gene_rocs[[1]][[2]], kaforou_rocs_list[[1]][[2]][[1]])
test_rocs_tb_nontb <- roc.test(gene_rocs[[2]], kaforou_rocs_list[[2]][[1]])

test_rocs_tb_od
test_rocs_tb_ltbi
test_rocs_tb_nontb
