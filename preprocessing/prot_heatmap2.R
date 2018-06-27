library(heatmap3)

df.prot_data <- read.csv("../data/protein_data_filtimp_k10.csv")

heatmap3(df.prot_data[,-c(1, 2, 3)],
         Rowv=FALSE,
         Colv=FALSE,
         margins=c(10,5),
         xlab="Proteins",
         ylab="Patients"
         )