setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_12/")

library(VennDiagram)

sel.gp.tb_od <- read.csv("../../data/ex_12/feat_sel_1_2/gp_tb_od_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.gp.tb_ltbi <- read.csv("../../data/ex_12/feat_sel_1_2/gp_tb_ltbi_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)
sel.gp.tb_nontb <- read.csv("../../data/ex_12/feat_sel_1_2/gp_tb_nontb_BH_LFC_lasso_sig_factors.csv", header=TRUE, row.names = 1)

date <- "2018-08-22/"
ex_dir <- "ex_12/"

  
venn.diagram(x = list(
  "TB vs OD" = sel.gp.tb_od$features,
  "TB vs LTBI" = sel.gp.tb_ltbi$features,
  "TB vs non-TB" = sel.gp.tb_nontb$features

),
filename = paste("../../img/", ex_dir, date, "lasso_sel_gp_venn.png", sep=""),
imagetype = "png",
col = "transparent",
fill = c("cornflowerblue", "green", "yellow"),
main = "Number of selected features"
)



