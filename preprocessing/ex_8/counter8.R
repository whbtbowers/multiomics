setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_7/")

#df.all.meta.val <- read.csv("../../data/ex_7/gene_prot_validation_meta.csv", header=TRUE, row.names = 1)
#df.all.meta.val$group <- as.character(df.all.meta.val$group)


datasets <- list(
  list(df.gene.val, df.all.meta.val,"gene validation")
  ,
  list(df.gene.tt, df.all.meta.tt, "gene train/test")
  ,
  list(df.prot.val, df.all.meta.val, "protein validation")
  ,
  list(df.prot.tt, df.all.meta.tt, "protein train/test")
)

print(paste("Number of expressed genes: ", ncol(df.gene.val), sep=""))
print(paste("Number of expressed proteins: ", ncol(df.prot.val), sep=""))

for (h in 1:length(datasets)){
  data <- datasets[[h]][1][[1]]
  meta <- datasets[[h]][2][[1]]
  label <- datasets[[h]][3][[1]]
  
  print(paste("Number of expressed components in ", label, " dataset: ", ncol(data), sep=""))
  
  print(paste("Total patients: ", nrow(meta), sep=""))
  for (i in c("HIV+", "HIV-")){
    print(paste(i, ": ", length(meta$hiv.status[meta$hiv.status == i]), sep=""))
  }
  
  for (j in 1:6){
    print(paste(j, ": ", length(meta$group[meta$group == j])))
  }
}

