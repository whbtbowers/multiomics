setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_9/")
#setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/ex_9/")

library(pROC)
library(ggplot2)

set.seed(12)

# To direct to the correct folder
date <- "2018-08-17/"
ex_dir <- "ex_10/"

# Features selected in Kaforou 2013
sel.gene.kaforou.tb_od <- read.csv("../../data/kaforou_2013/gene_tb_od_kaforou_2013.csv", header=TRUE, row.names = 1)
sel.gene.kaforou.tb_ltbi <- read.csv("../../data/kaforou_2013/gene_tb_ltbi_kaforou_2013.csv", header=TRUE, row.names = 1)
sel.gene.kaforou.tb_nontb <- read.csv("../../data/kaforou_2013/gene_tb_nontb_kaforou_2013.csv", header=TRUE, row.names = 1)

# Selected features for tb vs od
sel.gene.tb_od <- read.csv("../../data/ex_10/feat_sel_1_2/gene_tb_od_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.prot.tb_od <- read.csv("../../data/ex_10/feat_sel_1_2/prot_tb_od_hivneg_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.gp.tb_od <- rbind(sel.prot.tb_od, sel.gene.tb_od)

# Selected features for tb vs ltbi
sel.gene.tb_ltbi <- read.csv("../../data/ex_10/feat_sel_1_2/gene_tb_ltbi_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.prot.tb_ltbi <- read.csv("../../data/ex_10/feat_sel_1_2/prot_tb_ltbi_hivneg_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.gp.tb_ltbi <- rbind(sel.prot.tb_ltbi, sel.gene.tb_ltbi)

# Selected features for tb vs non-tb
sel.gene.tb_nontb <- read.csv("../../data/ex_10/feat_sel_1_2/gene_tb_nontb_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.prot.tb_nontb <- read.csv("../../data/ex_10/feat_sel_1_2/prot_tb_nontb_hivneg_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.gp.nontb <- rbind(sel.prot.tb_nontb, sel.gene.tb_nontb)

# Selected features for combined lasso

sel.p1.gp.tb_od <- read.csv("../../data/ex_10/feat_sel_1_2/gp_tb_od_hivneg_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)

sel.p1.gp.tb_ltbi <- read.csv("../../data/ex_10/feat_sel_1_2/gp_tb_ltbi_hivneg_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)

sel.p1.gp.tb_nontb <- read.csv("../../data/ex_10/feat_sel_1_2/gp_tb_nontb_hivneg_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)


# Selected features for phase 2 comparisons

sel.phase2.gp.tb_od <- read.csv("../../data/ex_9/feat_sel_2/tb_od_hivneg_BH_LFC_lasso_phase2_sig_factors.csv", header=TRUE, row.names = 1)

sel.phase2.gp.tb_ltbi <- read.csv("../../data/ex_9/feat_sel_2/tb_ltbi_hivneg_BH_LFC_lasso_phase2_sig_factors.csv", header=TRUE, row.names = 1)

sel.phase2.gp.nontb <- read.csv("../../data/ex_9/feat_sel_2/tb_nontb_hivneg_BH_LFC_lasso_phase2_sig_factors.csv", header=TRUE, row.names = 1)

# Complete datasets
#df.gene.all <- read.csv("../../data/ex_9/gene_train_body.csv", header=TRUE, row.names = 1)
df.prot.all <- read.csv("../../data/ex_9/prot_train_body.csv", header=TRUE, row.names = 1)
df.gp.all <- cbind(df.prot.all, df.gene.all)
#write.csv(df.gp.all, "../../data/ex_9/gp_train_body.csv", row.names=TRUE)
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

# Reconstitute probe ids so kaforou stuff can be searched by id
all.probe.ids <- c()

for (i in (1+ncol(df.prot.all)):length(df.gp.all)){
  id.parts <- strsplit(colnames(df.gp.all)[i], split="_")
  recon.probe.id <- paste(id.parts[[1]][1], "_",id.parts[[1]][2], sep = "")
  all.probe.ids <- c(all.probe.ids, recon.probe.id)
}

# lists for comparing roc curves

list.kaforou.rocs <- list()
list.my.rocs <- list()

# Comparisons lists
comps <- list(
  #list(sel.gene.tb_od, "", df.meta, "my", 1, 6, "TB", "OD", "HIV-","tb_od_hivneg")
  #,
  #list(sel.gene.tb_ltbi, "", df.meta, "my", 1, 3, "TB", "LTBI", "HIV-","tb_ltbi_hivneg")
  #,
  #list(sel.gene.tb_nontb, "", df.meta.tb_nontb, "my", 1, 7, "TB", "non-TB", "HIV-","tb_nontb_hivneg")
  #,
  list(sel.gene.tb_od, sel.prot.tb_ltbi, df.meta, "my", 1, 6, "TB", "OD", "HIV-","tb_od_hivneg")
  ,
  list(sel.gene.tb_ltbi, sel.prot.tb_ltbi, df.meta, "my", 1, 3, "TB", "LTBI", "HIV-","tb_ltbi_hivneg")
  ,
  list(sel.gene.tb_nontb, sel.prot.tb_nontb, df.meta.tb_nontb, "my", 1, 7, "TB", "non-TB", "HIV-","tb_nontb_hivneg")
  ,
  list(sel.p1.gp.tb_od, "", df.meta, "my_phase2", 1, 6, "TB", "OD", "HIV-", "gp_tb_od_hivneg")
  ,
  list(sel.p1.gp.tb_ltbi, "", df.meta, "my_phase2", 1, 3, "TB", "LTBI", "HIV-", "gp_tb_ltbi_hivneg")
  ,
  list(sel.p1.gp.tb_nontb, "", df.meta.tb_nontb, "my_phase2", 1, 7, "TB", "non-TB", "HIV-", "gp_tb_nontb_hivneg")
  ,
  #list(sel.phase2.gp.tb_od, "", df.meta.tb_od, "my_phase2", 1, 6, "TB", "OD", "HIV-", "od_nontb_hivneg_phase2")
  #,
  #list(sel.phase2.gp.tb_ltbi, "", df.meta, "my_phase2", 1, 3, "TB", "LTBI", "HIV-", "tb_ltbi_hivneg_phase2")
  #,
  #list(sel.phase2.gp.tb_nontb, "", df.meta.tb_nontb, "my_phase2", 1, 7, "TB", "non-TB", "HIV-", "tb_nontb_hivneg_phase2")
  #,
  list(sel.gene.kaforou.tb_od, "", df.meta.allhiv, "kaforou_2013", 9, 11, "TB", "OD", "HIV+/-", "tb_od_allhiv")
  ,
  list(sel.gene.kaforou.tb_ltbi, "", df.meta.allhiv, "kaforou_2013", 9, 10, "TB", "LTBI", "HIV+/-", "tb_ltbi_allhiv")
  ,
  list(sel.gene.kaforou.tb_nontb, "", df.meta.allhiv.tb_nontb, "kaforou_2013", 9, 12, "TB", "non-TB", "HIV+/-", "tb_nontb_allhiv")
)

for (comp in comps){
  comp.sel.genes <- comp[[1]]
  comp.sel.prots <- comp[[2]]
  comp.meta <- comp[[3]]
  comp.origin <- comp[[4]][1]
  comp.g1 <- comp[[5]][1]
  comp.g2 <- comp[[6]][1]
  comp.verb.g1 <- comp[[7]][1]
  comp.verb.g2 <- comp[[8]][1]
  comp.hivstat <- comp[[9]][1]
  comp.abbrv <- comp[[10]][1]
  
  if (comp.origin == "kaforou_2013"){
    
    ######################################
    ######################################
    ## Kaforou 2013 specific processing ##
    ######################################
    ######################################
    
    upreg.features <- c()
    downreg.features <- c()
    
    for (i in 1:nrow(comp.sel.genes)){
      if (comp.sel.genes$regulation_direction[i] == "Up"){
        upreg.features <- c(upreg.features, as.character(comp.sel.genes$probe_id[i]))
        #print(paste("UP:", sel.gene.kaforou.tb_od$probe_id[i]))
      } else {
        downreg.features <- c(downreg.features, as.character(comp.sel.genes$probe_id[i]))
        #print(paste("DOWN:", sel.gene.kaforou.tb_od$probe_id[i]))
      }
    }
    
    df.upreg <- df.gene.all[,match(upreg.features, all.probe.ids)]
    df.downreg <- df.gene.all[,match(downreg.features, all.probe.ids)]
    
    # Get DRS values of patients
    
    drs <- c()
    
    for (i in 1:nrow(df.gp.all)){
      drs.i <- sum(df.upreg[i,]) - sum(df.downreg[i,])
      drs <- c(drs, drs.i)
    }
    
    # Select DRS by TB status
    
    ind.g1 <- c()
    ind.g2 <- c()
    ind.g1_g2 <- c()
    
    
    for (j in 1:nrow(comp.meta)){
      if (comp.meta$group[j] == comp.g1){
        ind.g1 <- c(ind.g1, j)
      }
    }
    
    for (j in 1:nrow(comp.meta)){
      if (comp.meta$group[j] == comp.g2){
        ind.g2 <- c(ind.g2, j)
      }
    }
    
    for (j in 1:nrow(comp.meta)){
      if ((comp.meta$group[j] == comp.g1) || (comp.meta$group[j] == comp.g2)){
        ind.g1_g2 <- c(ind.g1_g2, j)
      }
    }
    
    comp.meta.g1_g2 <- comp.meta[ind.g1_g2,]
    
    drs.g1 <- drs[ind.g1]
    drs.g2 <- drs[ind.g2]
    drs.g1_g2 <- drs[ind.g1_g2]
    
    #Get means and sds of drs
    
    mean.drs.g1 <- mean(drs.g1)
    mean.drs.g2 <- mean(drs.g2)
    sd.drs.g1 <- sd(drs.g1)
    sd.drs.g2 <- sd(drs.g2)
    
    # Get threshold
    
    drs.threshold <- ((mean.drs.g1/sd.drs.g1)+(mean.drs.g2/sd.drs.g2))/((1/sd.drs.g1)+(1/sd.drs.g2))
    
    # Get ROC
    
    drs.roc <- roc(comp.meta.g1_g2$group, drs.g1_g2, auc=TRUE)
    list.kaforou.rocs <- list(list.kaforou.rocs, list(drs.roc, comp.abbrv))
    drs.roc.smooth <- smooth(drs.roc)
    
    png(paste("../../img/", ex_dir, date, comp.origin, "_", comp.abbrv, "_drs_roc.png", sep=""),
        width = 5*300,        # 5 x 300 pixels
        height = 5*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8) # smaller font size
    
    plot(drs.roc,
         col="red",
         lwd=1,
         main=paste("ROC curve for selected gene probes for", comp.hivstat, "patients," , comp.verb.g1,  "vs", comp.verb.g2, "(Kaforou et al 2013), (DRS)"),
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
    
    # Create dataframe and ggplot boxplot
    inf.group <- as.character(comp.meta.g1_g2$group)
    
    data = data.frame(inf.group, drs.g1_g2)
    
    ggplot(data, aes(x=inf.group, y=drs.g1_g2, group=inf.group)) + 
      geom_boxplot() +
      labs(x = "TB status", y="Disease risk score", title =paste("Distribution of DRS values for", comp.hivstat, "patients," , comp.verb.g1,  "vs", comp.verb.g2, "(Kaforou et al 2013)")) +
      scale_x_discrete(labels=c(comp.g1 = comp.verb.g1, comp.g2 = comp.verb.g2))
    
    ggsave(paste("../../img/", ex_dir, date, comp.origin, "_", comp.abbrv, "_drs_boxplot.png", sep=""))
    
  } else {
    
    ##############################################
    ##############################################
    ## My feature selection-specific processing ##
    ##############################################
    ##############################################
    
    # Select for correct patients for comparison
    
    ind.comp <- c()
    
    for (i in 1:nrow(df.gp.all)){
      if((comp.meta$group[i] == comp.g1) || (comp.meta$group[i] == comp.g2)){
        ind.comp <- c(ind.comp, i)
      }
    }
    
    df.meta.comp <- comp.meta[ind.comp,]
    
    if (comp.origin == "my_phase2"){
      
      ######################################
      # Process phase 2/ combination lasso #
      ######################################
      
      # Get data for selected features and betas
      
      set.sel.data <- df.gp.all[, match(comp.sel.genes$features, colnames(df.gp.all))]
      betas <- comp.sel.genes$beta
      
      set.sel.data.hivneg <- set.sel.data[ind.comp,]
      
      # Calculate beta adjustments
      
      b_adj <- c()
      
      for (i in 1:nrow(set.sel.data.hivneg)){
        b_adj.i <- c()
        for (j in 1:length(betas)){
          b_adj.i.j <- set.sel.data.hivneg[i, j] * betas[j]
          b_adj.i <- c(b_adj.i, b_adj.i.j)
        }
        
        b_adj <- c(b_adj, sum(b_adj.i))
      }
      
      # Get ROC
      
      b_adj.roc <- roc(df.meta.comp$group, b_adj, auc=TRUE)
      
      if (b_adj.roc$auc < 1){
        
        b_adj.roc.smooth <- smooth(b_adj.roc)
        
        png(paste("../../img/", ex_dir, date, comp.abbrv, "_b_adj_roc.png", sep=""),
            width = 5*300,        # 5 x 300 pixels
            height = 5*300,
            res = 300,            # 300 pixels per inch
            pointsize = 8) # smaller font size
        
        plot(b_adj.roc,
             col="red",
             lwd=1,
             main=paste("ROC curve for selected features for combined data,", comp.hivstat, "patients," , comp.verb.g1,  "vs", comp.verb.g2, "(beta-coefficient)"),
             legacy.axes = TRUE
        )
        
        lines.roc(b_adj.roc.smooth, col="blue", lwd=3)
        legend("bottomright",
               title = paste("AUC =", format(round(b_adj.roc$auc, 2), nsmall = 2)),
               legend = c("Empirical data", "Fit"),
               col = c("red", "blue"),
               lwd = c(1, 3)
        )
        
        dev.off()
        
      } else {
        
        png(paste("../../img/", ex_dir, date, set.abbrv, "_b_adj_roc.png", sep=""),
            width = 5*300,        # 5 x 300 pixels
            height = 5*300,
            res = 300,            # 300 pixels per inch
            pointsize = 8) # smaller font size
        
        plot(b_adj.roc,
             col="red",
             lwd=1,
             main=paste("ROC curve for selected features for combined data,", comp.hivstat, "patients," , comp.verb.g1,  "vs", comp.verb.g2, "(beta-coefficient)"),
             legacy.axes = TRUE
        )
        
        legend("bottomright",
               title = paste("AUC =", format(round(b_adj.roc$auc, 2), nsmall = 2)),
               legend = c("Empirical data"),
               col = c("red"),
               lwd = c(1, 3)
        )
        
        dev.off()
        
      }
      
      # Create dataframe and ggplot boxplot
      inf.group <- as.character(df.meta.comp$group)
      
      data = data.frame(inf.group, b_adj)
      
      ggplot(data, aes(x=inf.group, y=b_adj, group=inf.group)) + 
        geom_boxplot() +
        labs(x = "TB status", y="Beta adjustment score", title =paste("Distribution of Beta adjustment scores for combined data,", comp.hivstat, "patients," , comp.verb.g1,  "vs", comp.verb.g2, "(beta-coefficient)")) +
        scale_x_discrete(labels=c(comp.verb.g1, comp.verb.g2))
      
      ggsave(paste("../../img/", ex_dir, date, comp.abbrv, "_b_adj_boxplot.png", sep=""))
      
    } else {
      
      ###################
      # Process phase 1 #
      ###################
    
      sets <- list(
        list(comp.sel.genes, df.gene.all, "gene", "gene")
        #,
        #list(comp.sel.prots, df.prot.all, "protein", "prot")
      )
      
      for (set in sets){
        set.sel_features <- set[[1]]
        set.data <- set[[2]]
        set.verbose <- set[[3]][1]
        set.abbrv <- set[[4]][1]
        
        # Get data for selected features and betas
        
        set.sel.data <- set.data[, match(set.sel_features$features, colnames(set.data))]
        betas <- set.sel_features$beta
        
        set.sel.data.hivneg <- set.sel.data[ind.comp,]
  
        # Calculate beta adjustments
        
        b_adj <- c()
        
        for (i in 1:nrow(set.sel.data.hivneg)){
          b_adj.i <- c()
          for (j in 1:length(betas)){
            b_adj.i.j <- set.sel.data.hivneg[i, j] * betas[j]
            b_adj.i <- c(b_adj.i, b_adj.i.j)
          }
          
          b_adj <- c(b_adj, sum(b_adj.i))
        }
        
        # Get ROC
        
        b_adj.roc <- roc(df.meta.comp$group, b_adj, auc=TRUE)
        if (set.verbose == "gene"){
          list.my.rocs <- list(list.my.rocs, list(b_adj.roc, comp.abbrv))
        }
        
        if (b_adj.roc$auc < 1){
          
          b_adj.roc.smooth <- smooth(b_adj.roc)
          
          png(paste("../../img/", ex_dir, date, set.abbrv, "_", comp.abbrv, "_b_adj_roc.png", sep=""),
              width = 5*300,        # 5 x 300 pixels
              height = 5*300,
              res = 300,            # 300 pixels per inch
              pointsize = 8) # smaller font size
          
          plot(b_adj.roc,
               col="red",
               lwd=1,
               main=paste("ROC curve for selected features for", set.verbose, "data,", comp.hivstat, "patients," , comp.verb.g1,  "vs", comp.verb.g2, "(beta-coefficient)"),
               legacy.axes = TRUE
          )
          
          lines.roc(b_adj.roc.smooth, col="blue", lwd=3)
          legend("bottomright",
                 title = paste("AUC =", format(round(b_adj.roc$auc, 2), nsmall = 2)),
                 legend = c("Empirical data", "Fit"),
                 col = c("red", "blue"),
                 lwd = c(1, 3)
          )
          
          dev.off()
          
        } else {
          
          png(paste("../../img/", ex_dir, date, set.abbrv, "_", comp.abbrv, "_b_adj_roc.png", sep=""),
              width = 5*300,        # 5 x 300 pixels
              height = 5*300,
              res = 300,            # 300 pixels per inch
              pointsize = 8) # smaller font size
          
          plot(b_adj.roc,
               col="red",
               lwd=1,
               main=paste("ROC curve for selected features for", set.verbose, "data,", comp.hivstat, "patients," , comp.verb.g1,  "vs", comp.verb.g2, "(beta-coefficient)"),
               legacy.axes = TRUE
          )
          
          legend("bottomright",
                 title = paste("AUC =", format(round(b_adj.roc$auc, 2), nsmall = 2)),
                 legend = c("Empirical data"),
                 col = c("red"),
                 lwd = c(1, 3)
          )
          
          dev.off()
          
        }
        
        # Create dataframe and ggplot boxplot
        inf.group <- as.character(df.meta.comp$group)
        
        data = data.frame(inf.group, b_adj)
        
        ggplot(data, aes(x=inf.group, y=b_adj, group=inf.group)) + 
          geom_boxplot() +
          labs(x = "TB status", y="Beta adjustment score", title =paste("Distribution of Beta adjustment scores for", set.verbose, "data,", comp.hivstat, "patients," , comp.verb.g1,  "vs", comp.verb.g2, "(beta-coefficient)")) +
          scale_x_discrete(labels=c(comp.verb.g1, comp.verb.g2))
        
        ggsave(paste("../../img/", ex_dir, date, set.abbrv, "_", comp.abbrv, "_b_adj_boxplot.png", sep=""))
      }
    }
  }
}

# ROC curve comparisons
test.rocs.tb_od <- roc.test(list.my.rocs[[1]][[1]][[2]][[1]], list.kaforou.rocs[[1]][[1]][[2]][[1]])
test.rocs.tb_ltbi <- roc.test(list.my.rocs[[1]][[2]][[1]], list.kaforou.rocs[[1]][[2]][[1]])
test.rocs.tb_nontb <- roc.test(list.my.rocs[[2]][[1]], list.kaforou.rocs[[2]][[1]])

test.rocs.tb_od
test.rocs.tb_ltbi
test.rocs.tb_nontb

