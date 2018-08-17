setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_6")

library(VennDiagram)

df.gene.meta <- read.csv("../../data/ex_6/gene_data_meta.csv", header=TRUE, row.names = 1, na.strings = c("", "NA")) # Gene meta
df.prot.meta <- read.csv("../../data/ex_6/prot_data_meta.csv", header=TRUE, row.names = 1, na.strings = c("", "NA")) # Protein meta

# Remove patients with missing protein or array IDs

for (i in rev(1:nrow(df.gene.meta))){
  if((is.na(df.gene.meta$prot.id[i])) || (is.na(df.gene.meta$array.id[i]))){
    df.gene.meta <- df.gene.meta[-i,]
  }
}

for (i in rev(1:nrow(df.prot.meta))){
  if((is.na(df.prot.meta$prot.id[i])) || (is.na(df.prot.meta$array.id[i]))){
    df.prot.meta <- df.prot.meta[-i,]
  }
}

# Set date
date <- "2018-07-16"

testcols <- list(
  list("prot_id", 1),
  list("array_id", 2)
)

for (i in 1:length(testcols)){
  abrv <- testcols[[i]][[1]][1]
  colnum <- testcols[[i]][[2]][1]
  
  # Create venn diagram
  
  venn.diagram(x = list(
    "Gene data" = df.gene.meta[,colnum],
    "Protein data" = df.prot.meta[,colnum]
  ),
  filename = paste("../../img/ex_6/", date, "/venn_", abrv,".png", sep=""),
  imagetype= "png"
  )
}


