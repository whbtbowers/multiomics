df.gene_data <- read.csv("../data/gene_data_filtimp_k10.csv"=1)

subset.gene_data <- df.gene_data[1:100, 1:100]

write.csv(subset.gene_data, file="../data/100sq_subset_gene_data_filtimp_k10.csv")