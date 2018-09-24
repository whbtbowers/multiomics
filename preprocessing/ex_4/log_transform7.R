setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/ex_5/")

#library(limma)

df.prot.data <- read.csv("../../data/ex_4/prot_data_body_filtimp_k20.csv", row.names=1, header=TRUE)
df.gene.data <- read.csv("../../data/ex_4/gene_data_body.csv", row.names=1, header=TRUE)
# Log transform dataset

df.gene.log_trans <- data.frame(log(as.matrix(inp.prot)))


datasets <- list(
  list(df.gene.data, "gene", "gene")
  ,
  list(df.prot.data, "protein", "prot")
)

for (i in 1:length(datasets)){
  data <- datasets[[i]][[1]]
  verbose <- datasets[[i]][[2]]
  abbrv <- datasets[[i]][[3]]
  
  mostneg <- min(data)
  
  data.corr <- data + abs(mostneg)
  
  data.corr.log <- log(data.corr)
  
  write.csv(data.corr.log, file=paste("../../data/ex_5/", abbrv, "_data_body_logtransformed.csv", sep=""), row.names=TRUE)
  
}
