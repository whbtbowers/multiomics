setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/mix_ex")

df.prot_data.hiv_pos <- read.csv("../../data/gene_ex/gene_data_hiv_pos.csv", header=TRUE, row.names = 1)

df.prot_data.hiv_pos.body <- df.prot_data.hiv_pos[,-c(1,2)]

