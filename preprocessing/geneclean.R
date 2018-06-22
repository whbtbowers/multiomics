inp.data <- read.table("../data/expr_probes_adult_TB.txt", sep='\t', header=TRUE)
inp.data <- data.frame(inp.data)

sprintf("%d empty cells out of %d in dataset giving frequency of %f", sum(is.na(inp.data)), (nrow(inp.data)*ncol(inp.data)), ((sum(is.na(inp.data)))/(nrow(inp.data)*ncol(inp.data))))

row.del.count <- 0
col.del.count <- 0

sprintf("Initial data consist of %d columns and %d rows", ncol(inp.data), nrow(inp.data))

# Remove column <10% filled
for (i in rev(1:ncol(inp.data))){
  n.na.col <- sum(is.na(inp.data[,i]))
  f.na.col <- n.na.col/nrow(inp.data)
  if (f.na.col >= 0.9){
    inp.data <- inp.data[,-i]
    col.del.count <- col.del.count + 1
  }
}

sprintf("%d columns removed", col.del.count)
sprintf("Data consist of %d columns and %d rows", ncol(inp.data), nrow(inp.data))

# Remove rows <10% filled
for (i in rev(1:nrow(inp.data))){
  n.na.row <- sum(is.na(inp.data[i,]))
  f.na.row <- n.na.row/ncol(inp.data)
  if (f.na.col >= 0.9){
    inp.data <- inp.data[-i,]
    row.del.count <- row.del.count + 1
  }
  return
}

sprintf("%d rows removed", row.del.count)
sprintf("Data consist of %d columns and %d rows", ncol(inp.data), nrow(inp.data))

sprintf("%d empty cells out of %d in filtered dataset giving frequency of %f", sum(is.na(inp.data)), (nrow(inp.data)*ncol(inp.data)), ((sum(is.na(inp.data)))/(nrow(inp.data)*ncol(inp.data))))
