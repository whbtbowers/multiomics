setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_11/")

library(tictoc)
library(stringr)

tic()

#sets <- c('gene', 'prot')
sets <- c('gene')
comps <- c('tb_od', 'tb_ltbi', 'tb_nontb')

for (comp in comps){
  
  sig.gp.emn.factors <- read.csv(paste("../../data/ex_11/feat_sel_1_2/gp_", comp, "_hivneg_BH_LFC_lasso_sig_factors.csv", sep=""), header=TRUE, row.names = 1)
  #Get rownames
  sig.gp.emn.factors.names <- sig.gp.emn.factors$features
  
  #Count number of gene and protein factors in post-lasso
  sig.gp.emn.factors.names.genecount <- 0
  
  for (i in 1:length(sig.gp.emn.factors.names)){
    if (str_detect(toString(sig.gp.emn.factors.names[i]), "ILMN_")){
      sig.gp.emn.factors.names.genecount <- sig.gp.emn.factors.names.genecount + 1
    }
  }
  
  sig.gp.emn.factors.names.protcount <- length(sig.gp.emn.factors.names) - sig.gp.emn.factors.names.genecount
  
  # 
  
  print(paste("Number of factors for ", comp, " gene and protein data after lasso using lambda+1se as threshold: ", nrow(sig.gp.emn.factors), ". Number of genes: ", sig.gp.emn.factors.names.genecount, ". Number of proteins: ", sig.gp.emn.factors.names.protcount, sep=""))

  for (set in sets){
      
      sig.factors <- read.csv(paste("../../data/ex_11/feat_sel_1_2/", set, "_", comp, "_BH_LFC_sig_factors.csv", sep=""), header=TRUE, row.names = 1)
      sig.emn.factors <- read.csv(paste("../../data/ex_11/feat_sel_1_2/", set, "_", comp, "_BH_LFC_lasso_sig_factors.csv", sep=""), header=TRUE, row.names = 1)
      
      print(paste("Number of factors for ", set, " and ", comp, " before lasso: ", nrow(sig.factors), ". After phase1 lasso using lambda+1se as threshold: ", nrow(sig.emn.factors), sep=""))
  }
}

toc()