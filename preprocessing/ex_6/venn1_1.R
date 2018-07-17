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

# Get indices of patients based on TB status
ind.prot.tb <- which(df.prot.meta$tb.status %in% "TB")
ind.gene.tb <- which(df.gene.meta$tb.status %in% "TB")

ind.prot.ltbi <- which(df.prot.meta$tb.status %in% "LTBI")
ind.gene.ltbi <- which(df.gene.meta$tb.status %in% "LTBI")

ind.prot.od <- which(df.prot.meta$tb.status %in% "TB")
ind.gene.od <- which(df.gene.meta$tb.status %in% "TB")

#Get patients based on TB status

df.prot.tb <- df.prot.meta[ind.prot.tb,]
df.gene.tb <- df.gene.meta[ind.gene.tb,]

df.prot.ltbi <- df.prot.meta[ind.prot.ltbi,]
df.gene.ltbi <- df.gene.meta[ind.gene.ltbi,]

df.prot.od <- df.prot.meta[ind.prot.od,]
df.gene.od <- df.gene.meta[ind.gene.od,]

# Set date
date <- "2018-07-16"

datasets <- list(
  list(df.prot.meta, df.gene.meta, ""),
  list(df.prot.tb, df.gene.tb, "_tb"),
  list(df.prot.ltbi, df.gene.ltbi, "_ltbi"),
  list(df.prot.od, df.gene.od, "_od")
)

for (h in 1:length(datasets)){
  
  inf.prot.meta <- datasets[[h]][[1]]
  inf.gene.meta <- datasets[[h]][[2]]
  inf.abrv <- datasets[[h]][[3]]
  
  testcols <- list(
    list("prot_id", 1),
    list("array_id", 2)
  )
  
  for (i in 1:length(testcols)){
    abrv <- testcols[[i]][[1]][1]
    colnum <- testcols[[i]][[2]][1]
    
    # Create venn diagram
    
    venn.diagram(x = list(
      "Gene data" = inf.gene.meta[,colnum],
      "Protein data" = inf.prot.meta[,colnum]
      ),
    filename = paste("../../img/ex_6/", date, "/venn_", abrv, inf.abrv, ".png", sep=""),
    imagetype= "png"
    )
  }
}


