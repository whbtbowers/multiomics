setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_10/")

tic()


#sets <- c('gene', 'prot')
sets <- c('gene')
comps <- c('_tb_od', '_tb_ltbi', '_tb_nontb')

for (set in sets){
  for (comp in comps){
    
    
      sig.factors <- read.csv(paste("../../data/ex_10/feat_sel_1_2/", set, comp, "_BH_LFC_sig_factors.csv", sep=""), header=TRUE, row.names = 1)
      
      sig.emn.factors <- read.csv(paste("../../data/ex_10/feat_sel_1_2/", set, comp, "_BH_LFC_lasso_sig_factors.csv", sep=""), header=TRUE, row.names = 1)
      
      print(paste("Number of factors for ", set, " and ", comp, " with BH and LFC: ", nrow(sig.factors), ". After lasso using lambda+1se as threshold: ", nrow(sig.emn.factors), sep=""))
  }
}

toc()
