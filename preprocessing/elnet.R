inp.data <- read.csv("../data/protein_data_filtimp_k10.csv")
#print(inp.data)

library("glmnet")
fit<-glmnet(inp.data)
