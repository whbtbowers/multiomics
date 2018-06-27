library(heatmap3)

subset.gene_data <- read.csv("../data/100sq_subset_gene_data_filtimp_k10.csv", header=TRUE, row.names=1)



heatmap3(subset.gene_data,
         Rowv=FALSE,
         Colv=FALSE,
         margins=c(10,5),
         xlab="Genes",
         ylab="Patients"
         )