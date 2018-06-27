library('ggplot2')

inp.data = read.csv("../data/protein_data.csv")

# print(head(inp.data))

print(dim(inp.data))

g = 5

sprintf("Data contains %d columns and %d rows", ncol(inp.data), nrow(inp.data))

# To report back number of columns and rows deleted
cols.del = 0

# Check for sparsely pop'd rows (<10%)
for (i in 1:ncol(inp.data)){
  col.na = sum(is.na(inp.data[i]))
  freq.na = (col.na/nrow(inp.data[i]))
  if (freq.na >= 0.9){
    print(i)
    cols.del = cols.del + 1
    col.ind = c(col.ind, inp.data[i])
    inp.data <- inp.data[-i,]
  }

}

print(cols.del)
