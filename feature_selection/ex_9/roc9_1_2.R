#setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_9/")
setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/ex_9/")

library(pROC)
library(ggplot2)


set.seed(12)

# To direct to the correct folder
date <- "2018-08-04/"
ex_dir <- "ex_9/"

# Features selected in Kaforou 2013
sel.gene.kaforou.tb_od <- read.csv("../../data/kaforou_2013/gene_tb_od_kaforou_2013.csv", header=TRUE, row.names = 1)

# Complete datasets
#df.gene.all <- read.csv("../../data/ex_9/gene_train_body.csv", header=TRUE, row.names = 1)
df.prot.all <- read.csv("../../data/ex_9/prot_train_body.csv", header=TRUE, row.names = 1)
df.gp.all <- cbind(df.prot.all, df.gene.all)
df.meta <- read.csv("../../data/ex_9/gp_train_meta.csv", header=TRUE, row.names = 1)

# Selected features for tb vs od
sel.gene.tb_od <- read.csv("../../data/ex_9/feat_sel/gene_tb_od_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.prot.tb_od <- read.csv("../../data/ex_9/feat_sel/prot_tb_od_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.gp.tb_od <- rbind(sel.prot.tb_od, sel.gene.tb_od)

# Reconstitute probe ids so kaforou stuff can be searched by id
all.probe.ids <- c()

for (i in (1+ncol(df.prot.all)):length(df.gp.all)){
  id.parts <- strsplit(colnames(df.gp.all)[i], split="_")
  recon.probe.id <- paste(id.parts[[1]][1], "_",id.parts[[1]][2], sep = "")
  all.probe.ids <- c(all.probe.ids, recon.probe.id)
}

#Get upreg and downreg factors for kaforou 2013

kaforou.upreg.factors <- c()
kaforou.downreg.factors <- c()

for (i in 1:nrow(sel.gene.kaforou.tb_od)){
  if (sel.gene.kaforou.tb_od$regulation_direction[i] == "Up"){
    kaforou.upreg.factors <- c(kaforou.upreg.factors, as.character(sel.gene.kaforou.tb_od$probe_id[i]))
    #print(paste("UP:", sel.gene.kaforou.tb_od$probe_id[i]))
  } else {
    kaforou.downreg.factors <- c(kaforou.downreg.factors, as.character(sel.gene.kaforou.tb_od$probe_id[i]))
    #print(paste("DOWN:", sel.gene.kaforou.tb_od$probe_id[i]))
  }
}

df.upreg.kaforou.tb_od <- df.gp.all[,match(kaforou.upreg.factors, all.probe.ids)]
df.downreg.kaforou.tb_od <- df.gp.all[,match(kaforou.downreg.factors, all.probe.ids)]

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

# Get DRS values
drs.kaforou <- c()
drs.my <- c()

for (i in 1:nrow(df.gp.all)){
  drs.kaforou.i <- sum(df.upreg.kaforou.tb_od[i,]) - sum(df.downreg.kaforou.tb_od[i,])
  drs.my.i <- sum(df.upreg.my.tb_od[i,]) - sum(df.downreg.my.tb_od[i,])
  drs.kaforou <- c(drs.kaforou, drs.kaforou.i)
  drs.my <- c(drs.my, drs.my.i)
}

# Get DRS values of HIV- tb and od patients

ind.hiv_neg.tb <- c()
ind.hiv_neg.od <- c()
ind.hiv_neg.tb_od <- c()


for (j in 1:nrow(df.prot.all)){
  if (df.meta$group[j] == 1){
    ind.hiv_neg.tb <- c(ind.hiv_neg.tb, j)
  }
}

for (j in 1:nrow(df.prot.all)){
  if (df.meta$group[j] == 6){
    ind.hiv_neg.od <- c(ind.hiv_neg.od, j)
  }
}

for (j in 1:nrow(df.prot.all)){
  if ((df.meta$group[j] == 1) || (df.meta$group[j] == 6)){
    ind.hiv_neg.tb_od <- c(ind.hiv_neg.tb_od, j)
  }
}

df.meta.hiv_neg.tb_od <- df.meta[ind.hiv_neg.tb_od,]

drs.kaforou.hiv_neg.tb <- drs.kaforou[ind.hiv_neg.tb]
drs.kaforou.hiv_neg.od <- drs.kaforou[ind.hiv_neg.od]
drs.kaforou.hiv_neg.tb_od <- drs.kaforou[ind.hiv_neg.tb_od]


drs.my.hiv_neg.tb <- drs.my[ind.hiv_neg.tb]
drs.my.hiv_neg.od <- drs.my[ind.hiv_neg.od]
drs.my.hiv_neg.tb_od <- drs.my[ind.hiv_neg.tb_od]

#Get mean and sd of drs
mean.drs.kaforou.hiv_neg.tb <- mean(drs.kaforou.hiv_neg.tb)
mean.drs.kaforou.hiv_neg.od <- mean(drs.kaforou.hiv_neg.od)
sd.drs.kaforou.hiv_neg.tb <- sd(drs.kaforou.hiv_neg.tb)
sd.drs.kaforou.hiv_neg.od <- sd(drs.kaforou.hiv_neg.od)


mean.drs.my.hiv_neg.tb <- mean(drs.my.hiv_neg.tb)
mean.drs.my.hiv_neg.od <- mean(drs.my.hiv_neg.od)
sd.drs.my.hiv_neg.tb <- sd(drs.my.hiv_neg.tb)
sd.drs.my.hiv_neg.od <- sd(drs.my.hiv_neg.od)

#Get threshold values
threshold.kaforou <- ((mean.drs.kaforou.hiv_neg.tb/sd.drs.kaforou.hiv_neg.tb)+(mean.drs.kaforou.hiv_neg.od/sd.drs.kaforou.hiv_neg.od))/((1/sd.drs.kaforou.hiv_neg.tb)+(1/sd.drs.kaforou.hiv_neg.od))

threshold.my <- ((mean.drs.my.hiv_neg.tb/sd.drs.my.hiv_neg.tb)+(mean.drs.my.hiv_neg.od/sd.drs.my.hiv_neg.od))/((1/sd.drs.my.hiv_neg.tb)+(1/sd.drs.my.hiv_neg.od))

# Get ROC based on DRS

kaforou.roc <- roc(df.meta.hiv_neg.tb_od$group, drs.kaforou.hiv_neg.tb_od)
kaforou.roc.smooth <- smooth(kaforou.roc)
png(paste("../../img/", ex_dir, date,"kaforou_2013_tb_od_roc.png", sep=""),
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) # smaller font size
plot(kaforou.roc.smooth, main="ROC curve for Kaforou 2013 selected factors, TB vs OD")
dev.off()

my.roc <- roc(df.meta.hiv_neg.tb_od$group, drs.my.hiv_neg.tb_od)
my.roc.smooth <- smooth(my.roc)

png(paste("../../img/", ex_dir, date,"tb_od_roc.png", sep=""),
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8) # smaller font size

plot(my.roc.smooth, main="ROC curve for selected factors, TB vs OD")
dev.off()

#Get details of ROCS
auc(my.roc)
