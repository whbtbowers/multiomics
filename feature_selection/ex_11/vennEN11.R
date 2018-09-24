setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_11/")

library(VennDiagram)

sel.gene.tb_od <- read.csv("../../data/ex_11/feat_sel_2/gene_tb_od_hivneg_BH_LFC_EN_sig_factors.csv", header=TRUE, row.names = 1)
sel.gene.tb_ltbi <- read.csv("../../data/ex_11/feat_sel_2/gene_tb_ltbi_hivneg_BH_EN_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.gene.tb_nontb <- read.csv("../../data/ex_11/feat_sel_2/gene_tb_nontb_hivneg_BH_LFC_EN_sig_factors.csv", header=TRUE, row.names = 1)

date <- "2018-08-20/"
ex_dir <- "ex_11/"

  
venn.diagram(x = list(
  "TB vs OD" = sel.gene.tb_od$features,
  "TB vs LTBI" = sel.gene.tb_ltbi$features,
  "TB vs non-TB" = sel.gene.tb_nontb$features
  
),
filename = paste("../../img/", ex_dir, date, "_EN_sel_gene_venn.png", sep=""),
imagetype = "png",
col = "transparent",
fill = c("cornflowerblue", "green", "yellow"),
main = "Number of selected gene probes for each comparison"
)



