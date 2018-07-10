setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/")

library(illuminaio)
library(stringr)

inp.probe_labels <- readBGX(file.path("../data/HumanHT-12_V4_0_R2_15002873_B.bgx"))

inp.gene_data <- read.csv("../data/gene_ex2/gene_data_ex.csv",row.names=1, header=TRUE)
inp.gene_data.body <- inp.gene_data[,-c(1:7)]


df.probe_labels.probes <- inp.probe_labels$probes
df.probe_labels.controls <-inp.probe_labels$controls

start_time <- Sys.time()

ref_id <- df.probe_labels.probes$Array_Address_Id

labels.gene_data <- colnames(inp.gene_data.body)

#Get gene names
ind.gene_labels <- c()
gene_labels <- c()

for (i in 1:length(labels.gene_data)){
  ind <- match(as.numeric(substring(labels.gene_data[i],2)), ref_id)
  ind.gene_labels <- c(ind.gene_labels, ind)
  #gene_labels <- c(gene_labels, df.probe_labels.probes$ILMN_Gene[ind])
}

gene_labels <- df.probe_labels.probes$ILMN_Gene[ind.gene_labels]

colnames(inp.gene_data)[8:ncol(inp.gene_data)] <- gene_labels

write.csv(inp.gene_data, "../data/gene_ex2/labd_gene_data_ex.csv", row.names=TRUE)

inp.gene_data.body.scaled <- scale(inp.gene_data[,-c(1:7)])

inp.gene_data.scaled <- cbind(inp.gene_data[,c(1:7)],inp.gene_data.body.scaled)

write.csv(inp.gene_data.scaled, "../data/gene_ex2/labd_gene_data_ex_scaled.csv", row.names=TRUE)
