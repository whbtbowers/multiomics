setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/ex_4/")

inp.prot <- read.csv("../../data/ex_4/gene_data_body.csv", row.names=1, header=TRUE)

# Log transform dataset

df.gene.log_trans <- log(inp.prot)

