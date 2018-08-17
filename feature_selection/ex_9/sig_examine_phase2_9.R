setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_9/")

library(tictoc)
library(stringr)

tic()

sets <- c('gene', 'prot')
comps <- c('tb_od', 'tb_ltbi', 'tb_nontb')

for (comp in comps){
  
  sig.phase2.factors <- read.csv(paste("../../data/ex_9/feat_sel_2/", comp, "_hivneg_BH_LFC_phase2_sig_factors.csv", sep=""), header=TRUE, row.names = 1)
  
  sig.phase2.emn.factors <- read.csv(paste("../../data/ex_9/feat_sel_2/", comp, "_hivneg_BH_LFC_lasso_phase2_sig_factors.csv", sep=""), header=TRUE, row.names = 1)
  
  #Get rownames
  sig.phase2.factors.names <- rownames(sig.phase2.factors)
  sig.phase2.emn.factors.names <- sig.phase2.emn.factors$features
  
  #Count number of gene and protein factors in pre-lasso
  sig.phase2.factors.names.genecount <- 0
  
  for (i in 1:length(sig.phase2.factors.names)){
    if (str_detect(toString(sig.phase2.factors.names[i]), "ILMN_")){
      sig.phase2.factors.names.genecount <- sig.phase2.factors.names.genecount + 1
    }
  }
  
  sig.phase2.factors.names.protcount <- length(sig.phase2.factors.names) - sig.phase2.factors.names.genecount
  
  #Count number of gene and protein factors in post-lasso
  sig.phase2.emn.factors.names.genecount <- 0
  
  for (i in 1:length(sig.phase2.emn.factors.names)){
    if (str_detect(toString(sig.phase2.emn.factors.names[i]), "ILMN_")){
      sig.phase2.emn.factors.names.genecount <- sig.phase2.emn.factors.names.genecount + 1
    }
  }
  
  sig.phase2.emn.factors.names.protcount <- length(sig.phase2.emn.factors.names) - sig.phase2.emn.factors.names.genecount
  
  # 
  
  print(paste("Number of factors for phase 2 ", comp, " with BH and LFC: ", nrow(sig.phase2.factors), ". Number of genes: ", sig.phase2.factors.names.genecount, ". Number of proteins: ", sig.phase2.factors.names.protcount, ". After lasso using lambda+1se as threshold: ", nrow(sig.phase2.emn.factors), ". Number of genes: ", sig.phase2.emn.factors.names.genecount, ". Number of proteins: ", sig.phase2.emn.factors.names.protcount, sep=""))

  for (set in sets){
      
      sig.emn.factors <- read.csv(paste("../../data/ex_9/feat_sel/", set, "_", comp, "_BH_LFC_lasso_sig_factors.csv", sep=""), header=TRUE, row.names = 1)
      
      print(paste("Number of factors for ", set, " and ", comp, " after phase1 lasso using lambda+1se as threshold: ", nrow(sig.emn.factors), sep=""))
  }
}



toc()


sig.factors <- read.csv(paste("../../data/ex_9/feat_sel/", set, comp, "_BH_LFC_sig_factors.csv", sep=""), header=TRUE, row.names = 1)