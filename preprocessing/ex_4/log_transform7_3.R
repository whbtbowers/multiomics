setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_7/")

#Initial datasets
df.gene.body <- read.csv("../../data/ex_7/gene_data_body.csv", header=TRUE, row.names = 1)
df.prot.body <- read.csv("../../data/ex_7/prot_data_body.csv", header=TRUE, row.names = 1)

# Normalised datasets
#df.gene.body.wn <- read.csv("../../data/ex_7/gene_data_body_wn.csv", header=TRUE, row.names = 1)
#df.prot.body.wn <- read.csv("../../data/ex_7/prot_data_body_wn.csv", header=TRUE, row.names = 1)

library(preprocessCore)

datasets <- list(
  list(df.gene.body, "gene", "gene")
  ,
  list(df.prot.body, "protein", "prot")
  #,
  #list(df.gene.body.wn, "normalised gene", "gene_wn")
  #,
  #list(df.prot.body.wn, "normalised protein", "prot_wn")
)

for (i in 1:length(datasets)){
  data <- datasets[[i]][[1]]
  verbose <- datasets[[i]][[2]]
  abbrv <- datasets[[i]][[3]]
  
  mostneg <- min(data)

  data.corr <- data + abs(mostneg) + (1e-05)

  data.corr.logged <- log2(data.corr)
  
  write.csv(data.corr.logged, file=paste("../../data/ex_7/", abbrv, "_data_body_lt2.csv", sep=""), row.names=TRUE)
  
  # Quantile normalisation
  
  data.corr.logged.qn <- as.data.frame(normalize.quantiles(as.matrix(data.corr.logged)), 
                                       row.names = rownames(data), 
                                       col.names = colnames(data)
  )
  colnames(data.corr.logged.qn) <- colnames(data)
  
  write.csv(data.corr.logged.qn, file=paste("../../data/ex_7/", abbrv, "_data_body_lt2_qn.csv", sep=""), row.names=TRUE)
  
}
