setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_9/")
#setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/ex_9/")

library(pROC)
library(ggplot2)


set.seed(12)

# To direct to the correct folder
date <- "2018-08-07/"
ex_dir <- "ex_9/"

# Features selected in Kaforou 2013
sel.gene.kaforou.tb_od <- read.csv("../../data/kaforou_2013/gene_tb_od_kaforou_2013.csv", header=TRUE, row.names = 1)

# Complete datasets
#df.gene.all <- read.csv("../../data/ex_9/gene_train_body.csv", header=TRUE, row.names = 1)
df.prot.all <- read.csv("../../data/ex_9/prot_train_body.csv", header=TRUE, row.names = 1)
df.gp.all <- cbind(df.prot.all, df.gene.all)
df.meta <- read.csv("../../data/ex_9/gp_train_meta.csv", header=TRUE, row.names = 1)
df.meta$group <- as.character(df.meta$group)

# Selected features for tb vs od
sel.gene.tb_od <- read.csv("../../data/ex_9/feat_sel/gene_tb_od_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.prot.tb_od <- read.csv("../../data/ex_9/feat_sel/prot_tb_od_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.gp.tb_od <- rbind(sel.prot.tb_od, sel.gene.tb_od)

# Reconstitute probe ids so kaforou stuff can be searched by id
all.probe.ids <- c()

for (i in 1:length(df.gene.all)){
  id.parts <- strsplit(colnames(df.gene.all)[i], split="_")
  recon.probe.id <- paste(id.parts[[1]][1], "_",id.parts[[1]][2], sep = "")
  all.probe.ids <- c(all.probe.ids, recon.probe.id)
}

#Get upreg and downreg factors for mine

my.upreg.factors <- c()
my.downreg.factors <- c()

for (i in 1:nrow(sel.gp.tb_od)){
  if (sel.gp.tb_od$reg_dir[i] == "up"){
    my.upreg.factors <- c(my.upreg.factors, as.character(sel.gp.tb_od$features[i]))
    #print(paste("UP:", sel.gp.tb_od$features[i]))
  } else {
    my.downreg.factors <- c(my.downreg.factors, as.character(sel.gp.tb_od$features[i]))
    #print(paste("DOWN:", sel.gp.tb_od$features[i]))
  }
}

df.upreg.my.tb_od <- df.gp.all[match(my.upreg.factors, colnames(df.gp.all))]
df.downreg.my.tb_od <- df.gp.all[match(my.downreg.factors, colnames(df.gp.all))]


