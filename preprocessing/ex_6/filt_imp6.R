setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_6/")

inp.data <- read.csv("../../data/ex_6/prot_data_body.csv", header=TRUE, row.names = 1) # Original

sprintf("%d empty cells out of %d in dataset giving frequency of %f", sum(is.na(inp.data)), (nrow(inp.data)*ncol(inp.data)), ((sum(is.na(inp.data)))/(nrow(inp.data)*ncol(inp.data))))

row.del.count <- 0
col.del.count <- 0

sprintf("Initial data consist of %d columns and %d rows", ncol(inp.data), nrow(inp.data))

# Remove column <20% filled
for (i in rev(1:ncol(inp.data))){
  n.na.col <- sum(is.na(inp.data[,i]))
  f.na.col <- n.na.col/nrow(inp.data)
  if (f.na.col >= 0.8){
    print(i)
    inp.data <- inp.data[,-i]
    col.del.count <- col.del.count + 1
  }
}

sprintf("%d columns removed", col.del.count)
sprintf("Data consist of %d columns and %d rows", ncol(inp.data), nrow(inp.data))

# Remove rows <20% filled
for (i in rev(1:nrow(inp.data))){
  n.na.row <- sum(is.na(inp.data[i,]))
  f.na.row <- n.na.row/ncol(inp.data)
  if (f.na.col >= 0.8){
    inp.data <- inp.data[-i,]
    row.del.count <- row.del.count + 1
  }
  return
}

sprintf("%d rows removed", row.del.count)
sprintf("Data consist of %d columns and %d rows", ncol(inp.data), nrow(inp.data))

sprintf("%d empty cells out of %d in filtered dataset giving frequency of %f", sum(is.na(inp.data)), (nrow(inp.data)*ncol(inp.data)), ((sum(is.na(inp.data)))/(nrow(inp.data)*ncol(inp.data))))



library(bnstruct)


#Impute remaining empty cells
print("Initiating imputation")

#X.t <- data.frame(t(inp.data))

#X.imp <- data.frame(knn.impute(as.matrix(X.t)), k=20) #Use k=20 like Wang
X.imp <- t(data.frame(knn.impute(as.matrix(t(inp.data))), k=20)) #Use k=20 like Wang

print("Writing to CSV")
write.csv(X.imp[-nrow(X.imp),], file="../../data/ex_6/prot_data_body_logtransformed_filtimp_k20.csv", row.names=TRUE)
