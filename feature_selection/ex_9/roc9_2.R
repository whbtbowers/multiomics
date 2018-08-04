setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_9/")

library(ROC)
library(pROC)
library(tictoc)

tic()
set.seed(12)

# To direct to the correct folder
date <- "2018-08-02/"
ex_dir <- "ex_9/"

# Features selected in Kaforou 2013
sel.gene.kaforou.tb_od <- read.csv("../../data/kaforou_2013/gene_tb_od_kaforou_2013.csv", header=TRUE, row.names = 1)
sel.gene.kaforou.tb_ltbi <- read.csv("../../data/kaforou_2013/gene_tb_ltbi_kaforou_2013.csv", header=TRUE, row.names = 1)
sel.gene.kaforou.tb_nontb <- read.csv("../../data/kaforou_2013/gene_tb_nontb_kaforou_2013.csv", header=TRUE, row.names = 1)

# Complete datasets
#df.gene.all <- read.csv("../../data/ex_9/gene_train_body.csv", header=TRUE, row.names = 1)
df.prot.all <- read.csv("../../data/ex_9/prot_train_body.csv", header=TRUE, row.names = 1)
df.meta <- read.csv("../../data/ex_9/gp_train_meta.csv", header=TRUE, row.names = 1)

# Selected features for tb vs od
sel.gene.tb_od <- read.csv("../../data/ex_9/feat_sel/gene_tb_od_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.prot.tb_od <- read.csv("../../data/ex_9/feat_sel/prot_tb_od_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)

# Selected features for tb vs ltbi
sel.gene.tb_ltbi <- read.csv("../../data/ex_9/feat_sel/gene_tb_ltbi_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.prot.tb_ltbi <- read.csv("../../data/ex_9/feat_sel/prot_tb_ltbi_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)

# Selected features for TB vs non-TB
sel.gene.tb_nontb <- read.csv("../../data/ex_9/feat_sel/gene_tb_nontb_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.prot.tb_nontb <- read.csv("../../data/ex_9/feat_sel/prot_tb_nontb_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)

# Reconstitute probe ids so kaforou stuff can be searched by id
all.probe.ids <- c()

for (i in 1:length(df.gene.all)){
  id.parts <- strsplit(colnames(df.gene.all)[i], split="_")
  recon.probe.id <- paste(id.parts[[1]][1], "_",id.parts[[1]][2], sep = "")
  recon.gene.sym <- id.parts[[1]][3]
  all.probe.ids <- c(all.probe.ids, recon.probe.id)
}

comp.sets <- list(
  list(sel.gene.tb_od, sel.prot.tb_od, sel.gene.kaforou.tb_od, "TB vs OD", "tb_od"),
  list(sel.gene.tb_ltbi, sel.prot.tb_ltbi, sel.gene.kaforou.tb_ltbi, "TB vs LTBI", "tb_ltbi"),
  list(sel.gene.tb_nontb, sel.prot.tb_nontb, sel.gene.kaforou.tb_nontb, "TB vs non-TB", "tb_nontb")
)

for (comp.set in comp.sets){
  sel.gene <- comp.set[1][[1]]
  sel.prot <- comp.set[2][[1]]
  sel.gene.kaforou <- comp.set[3][[1]]
  comp.verbose <- comp.set[4][[1]]
  comp.abbrv <- comp.set[5][[1]]
  
  # Examine which genes are shared between mine and kaforou 2013
  #sel.gene.tb_od$features
  sel.genes <- c()
  
  for (i in 1:length(sel.gene$features)){
    id.parts <- strsplit(toString(sel.gene$features[i]), split="_")
    recon.gene.sym <- id.parts[[1]][3]
    sel.genes <- c(sel.genes, recon.gene.sym)
  }
  
  counter <- 0
  
  for (i in 1:length(sel.genes)){
    if (sel.genes[i] %in% as.character(sel.gene.kaforou$gene_symbol)){
      counter <- counter + 1
      print(sel.genes[i])
    }
  }
  
  print(paste(counter, " overlapping genes in both feature selections for ", comp.verbose, sep=""))
  
  # Create binary response vector of sig genes from kaforou data
  # All probes which appear in kaforou selection labelled 1
  
  bin.vec.kaforou <- rep(0, ncol(df.gene.all))
  bin.vec.kaforou[match(sel.gene.kaforou$probe_id, all.probe.ids)] <- 1
  
  # My significant genes
  
  bin.vec <- rep(0, ncol(df.gene.all))
  bin.vec[match(sel.gene$features, colnames(df.gene.all))] <- 1
  
  # Draw ROC curve
  roc.tb_od <- roc(bin.vec, bin.vec.kaforou)
  
  png(paste("../../img/", ex_dir, date, comp.abbrv, "_roc.png", sep=""),
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  
  
  plot.roc(x=roc.tb_od, legacy.axes = TRUE, main = "TB vs OD")
  dev.off()
  
}
toc()
