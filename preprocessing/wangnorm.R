setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/ex_4/")

df.gene <- read.csv("../../data/ex_4/gene_data_body.csv", row.names=1, header=TRUE)

df.gene.wnorm <- data.frame()

for (i in 1:ncol(df.gene)){
  f <- df.gene[, i]
  E.f <- mean(f)
  var.f <- var(f)
 
  # Wang normalise
  
  f.wnorm <- data.frame((f-E.f)/sqrt(var.f))
  
  if (ncol(df.gene.wnorm) == 0){
    df.gene.wnorm <- f.wnorm
  } else {
    df.gene.wnorm <- cbind(df.gene.wnorm, f.wnorm)
  }
  
}

# Add row and column names columns
colnames(df.gene.wnorm) <- colnames(df.gene)
rownames(df.gene.wnorm) <- rownames(df.gene)

# Write to CSV
write.csv(df.gene.wnorm, "../../data/ex_4/gene_data_body_wnormed.csv", row.names=TRUE)
