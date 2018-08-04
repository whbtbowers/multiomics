setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_9/")

tic()

library(ROC)
library(pROC)

set.seed(12)

# To direct to the correct folder
date <- "2018-08-01/"
ex_dir <- "ex_9/"

# Genes selected in Kaforou 2013
sel.gene.kaforou.tb_od <- read.csv("../../data/kaforou_2013/gene_tb_od_kaforou_2013.csv", header=TRUE, row.names = 1)

# Complete datasets
df.gene.all <- read.csv("../../data/ex_9/gene_train_body.csv", header=TRUE, row.names = 1)
df.prot.all <- read.csv("../../data/ex_9/prot_train_body.csv", header=TRUE, row.names = 1)
df.meta <- read.csv("../../data/ex_9/gp_train_meta.csv", header=TRUE, row.names = 1)

# Selected features
sel.gene.tb_od <- read.csv("../../data/ex_9/feat_sel/gene_tb_od_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)

sel.prot.tb_od <- read.csv("../../data/ex_9/feat_sel/prot_tb_od_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)

# Reconstitute probe ids so kaforou stuff can be searched by id
all.probe.ids <- c()

for (i in 1:length(head(df.gene.all))){
  id.parts <- strsplit(colnames(df.gene.all)[i], split="_")
  recon.probe.id <- paste(id.parts[[1]][1], "_",id.parts[[1]][2], sep = "")
  all.probe.ids <- c(all.probe.ids, recon.probe.id)
}

# Create binary response vector of sig genes from kaforou data
# All probes which appear in kaforou selection labelled 1

bin.vec.gene.kaforou <- rep(0, ncol(df.gene.all))
bin.vec.gene.kaforou[match(sel.gene.kaforou.tb_od$probe_id, all.probe.ids)] <- 1

# No protein data in kaforou data, so just adding vecor of all 0s.

bin.vec.prot.kaforou <- rep(0, ncol(df.prot.all))

bin.vec.kaforou <- c(bin.vec.prot.kaforou, bin.vec.gene.kaforou)

# Create binary predictor vector from my data.
# Significant proteins

bin.vec.prot <- rep(0, ncol(df.prot.all))
bin.vec.prot[match(sel.prot.tb_od$features, colnames(df.prot.all))] <- 1

bin.vec.gene <- rep(0, ncol(df.gene.all))
bin.vec.gene[match(sel.gene.tb_od$features, colnames(df.gene.all))] <- 1

bin.vec <- c(bin.vec.prot, bin.vec.gene)

# Draw ROC curve
roc.tb_od <- roc(bin.vec, bin.vec.kaforou)

# Principal components 1 & 2 by location

png(paste("../../img/", ex_dir, date,"tb_od_roc.png", sep=""),
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size


plot.roc(x=roc.tb_od, legacy.axes = TRUE, main = "TB vs OD")
dev.off()



# Get runtime
toc()