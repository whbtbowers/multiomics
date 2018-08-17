setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_9/")

library(ggplot2)
library(stats)
library(ggfortify)

set.seed(12)

tic()

# Combined protein and gene probe set
df.gp.all <- read.csv("../../data/ex_9/gp_train_body.csv", header=TRUE, row.names = 1)  #Gene test/train set

df.meta <- read.csv("../../data/ex_9/gp_train_meta.csv", header=TRUE, row.names = 1)
df.meta$group <- as.character(df.meta$group)

# Selected features for tb vs od

sel.gp.tb_od <- read.csv("../../data/ex_9/feat_sel/gp_tb_od_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)

# Selected features for tb vs ltbi

sel.gp.tb_ltbi <- read.csv("../../data/ex_9/feat_sel/gp_tb_ltbi_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)

# To direct to the correct folder
date <- "2018-08-07/"
ex_dir <- "ex_9/"

phase1s = list(
  list(sel.gp.tb_od, 1, 6, "TB", "OD", "HIV-", "tb_od_hivneg")
  ,
  list(sel.gp.tb_ltbi, 1, 3, "TB", "LTBI", "HIV-", "tb_ltbi_hivneg")
)

#for (phase1 in phase1s){
#  phase1.sel_feat <- phase1[[1]][1]
#  phase1.g1 <- phase1[[2]][1]
#  phase1.g2 <- phase1[[3]][1]
#  phase1.verb.g1 <- phase1[[4]][1]
#  phase1.verb.g2 <- phase1[[5]][1]
#  phase1.hivstat <- phase1[[6]][1]
#  phase1.abbrv <- phase1[[7]][1]
  
  # Extract selected features from dataset
  
#  df.sel.feat <- df.gp.all[match(phase1.sel_feat$features, colnames(df.gp.all))]

# Select HIV- patients

ind.hiv_neg <- c()

for (i in 1:nrow(df.gp.all)){
  if((df.meta$group[i] == 1) || (df.meta$group[i] == 3) || (df.meta$group[i] == 6)){
    ind.hiv_neg <- c(ind.hiv_neg, i)
  }
}

df.data.hiv_neg <- df.gp.all[ind.hiv_neg,]
df.meta.hiv_neg <- df.meta[ind.hiv_neg,]

# Select HIV- TB vs OD patients, select appropriate features

ind.hiv_neg.tb_od <- c()

for (i in 1:nrow(df.gp.all)){
  if((df.meta$group[i] == 1) || (df.meta$group[i] == 6)){
    ind.hiv_neg.tb_od <- c(ind.hiv_neg.tb_od, i)
  }
}

df.sel.feat.hiv_neg.tb_od <- df.gp.all[,match(sel.gp.tb_od$features, colnames(df.gp.all))]

df.data.hiv_neg.tb_od <- df.sel.feat.hiv_neg.tb_od[ind.hiv_neg.tb_od,]
df.meta.hiv_neg.tb_od <- df.meta[ind.hiv_neg.tb_od,]

# Select HIV- TB vs LTBI patients, select appropriate features

ind.hiv_neg.tb_ltbi <- c()

for (i in 1:nrow(df.gp.all)){
  if((df.meta$group[i] == 1) || (df.meta$group[i] == 3)){
    ind.hiv_neg.tb_ltbi <- c(ind.hiv_neg.tb_ltbi, i)
  }
}

df.sel.feat.hiv_neg.tb_ltbi <- df.gp.all[,match(sel.gp.tb_ltbi$features, colnames(df.gp.all))]

df.data.hiv_neg.tb_ltbi <- df.sel.feat.hiv_neg.tb_ltbi[ind.hiv_neg.tb_ltbi,]
df.meta.hiv_neg.tb_ltbi <- df.meta[ind.hiv_neg.tb_ltbi,]

# Lazily select data and meta of choice here

#pca.data <- df.data.hiv_neg
#pca.meta <- df.meta.hiv_neg
#abbrv <- "gp_comb"
#verbose = "HIV- combined gene and protein"

#pca.data <- df.data.hiv_neg.tb_od
#pca.meta <- df.meta.hiv_neg.tb_od
#abbrv <- "sel_hivneg_gp_tb_od"
#verbose = "Selected features from HIV- TB vs OD, combined gene and protein"

pca.data <- df.data.hiv_neg.tb_ltbi
pca.meta <- df.meta.hiv_neg.tb_ltbi
abbrv <- "sel_hivneg_gp_tb_ltbi"
verbose = "Selected features from HIV- TB vs LTBI, combined gene and protein"

# Principal components 1 & 2 by location

png(paste("../../img/", ex_dir, date,   abbrv, "_pca_pc1_pc2_bysite.png", sep=""),
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size


autoplot(prcomp(pca.data),
         data =pca.meta,
         colour = "site",
         main=paste( verbose, "data by location", sep=" ")
)

dev.off()

# Principal components 3 & 4 by by location

png(paste("../../img/", ex_dir, date,   abbrv, "_pca_pc3_pc4_bysite.png", sep=""),
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(pca.data),
         x=3, 
         y=4,
         data = pca.meta,
         colour = "site",
         main=paste( verbose, "data by location", sep=" ")
)

dev.off()

# Principal components 1 & 2 by infection status

png(paste("../../img/", ex_dir, date,   abbrv, "_pca_pc1_pc2_byinf.png", sep=""),
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size


autoplot(prcomp(pca.data),
         data = pca.meta,
         colour = "group",
         main=paste( verbose, "data by TB status", sep=" ")
)

dev.off()

# Principal components 3 & 4 by infection status

png(paste("../../img/", ex_dir, date,   abbrv, "_pca_pc3_pc4_byinf.png", sep=""),
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(pca.data),
         x=3, 
         y=4,
         data = pca.meta,
         colour = "group",
         main=paste( verbose, "data by TB status", sep=" ")
)

dev.off()

