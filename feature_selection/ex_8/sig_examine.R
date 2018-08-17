setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_9/")
tic()

step1s <- c('_BH', '_LFC')
sets <- c('prot', 'gene')
comps <- c('_tb_ltbi', '_tb_od')

for (set in sets){
  for (comp in comps){
    for (step1 in step1s){
      sig.factors <- read.csv(paste("../../data/ex_9/feat_sel/", set, comp, step1, "_sig_factors.csv", sep=""), header=TRUE, row.names = 1)
      
      sig.emn.factors <- read.csv(paste("../../data/ex_9/feat_sel/", set, comp, step1, "_EMN_sig_factors.csv", sep=""), header=TRUE, row.names = 1)
      
      print(paste("Number of factors for ", set, " and ", comp, " with ", step1, ": ", nrow(sig.factors), ". After elastic net using lambda+1se as threshold: ", nrow(sig.emn.factors), sep=""))
    }
  }
}

toc()