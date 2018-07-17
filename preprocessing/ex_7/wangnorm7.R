setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_6")

df.prot.data <- read.csv("../../data/ex_6/prot_data_body_filtimp_k20.csv", header=TRUE, row.names = 1) # Corrected and log transformed protein data
df.gene.data <- read.csv("../../data/ex_6/gene_data_body.csv", header=TRUE, row.names = 1) # Corrected and log transformed protein data


datasets <- list(
  list(df.gene.data, "gene")
  ,
  list(df.prot.data, "prot")
)

for (i in 1:length(datasets)){
  data <- datasets[[i]][[1]]
  abbrv <- datasets[[i]][[2]]
  
  df.wnorm <- data.frame()
  
  for (i in 1:ncol(data)){
    f <- data[, i]
    E.f <- mean(f)
    var.f <- var(f)
    
    # Wang normalise
    
    f.wnorm <- data.frame((f-E.f)/sqrt(var.f))
    
    if (ncol(df.wnorm) == 0){
      df.wnorm <- f.wnorm
    } else {
      df.wnorm <- cbind(df.wnorm, f.wnorm)
    }
    
  }
  
  # Add row and column names columns
  colnames(df.wnorm) <- colnames(data)
  rownames(df.wnorm) <- rownames(data)
  
  # Write to CSV
  write.csv(df.wnorm, paste("../../data/ex_6/", abbrv, "_data_body_wnormed.csv", sep=""), row.names=TRUE)
  
}
