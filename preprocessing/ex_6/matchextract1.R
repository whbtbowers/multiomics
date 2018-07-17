setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_6")

# Gene data
df.gene.meta <- read.csv("../../data/ex_6/gene_data_meta.csv", header=TRUE, row.names = 1, na.strings = c("", "NA")) # Gene meta
df.gene.body <- read.csv("../../data/ex_6/gene_data_body.csv", header=TRUE, row.names = 1, na.strings = c("", "NA")) # Gene body

# Protein data
df.prot.meta <- read.csv("../../data/ex_6/prot_data_meta.csv", header=TRUE, row.names = 1, na.strings = c("", "NA")) # Protein meta
df.prot.body <-read.csv("../../data/ex_6/prot_data_body_filtimp_k20.csv", header=TRUE, row.names = 1, na.strings = c("", "NA")) # Gene meta

# Remove patients with missing protein or array IDs

for (i in rev(1:nrow(df.gene.meta))){
  if((is.na(df.gene.meta$prot.id[i])) || (is.na(df.gene.meta$array.id[i]))){
    df.gene.meta <- df.gene.meta[-i,]
    df.gene.body <- df.gene.body[-i,]
  }
}

for (i in rev(1:nrow(df.prot.meta))){
  if((is.na(df.prot.meta$prot.id[i])) || (is.na(df.prot.meta$array.id[i]))){
    df.prot.meta <- df.prot.meta[-i,]
    df.prot.body <- df.prot.body[-i,]
  }
}

# Vectors containing indices from protein and gene data
ind.gene.meta <- c()
ind.prot.meta <- c()

for (i in 1:nrow(df.gene.meta)){
  if (df.gene.meta$array.id[i] %in% df.prot.meta$array.id){
    ind.gene.meta <- c(ind.gene.meta, i)
    ind.prot.meta <- c(ind.prot.meta, match(df.gene.meta$array.id[i], df.prot.meta$array.id))
  }
} 

# Stick with gene meta, tack on protein data seeing as order already set.

# Extract patients from both datasets using indices
df.all.meta <- df.gene.meta[ind.gene.meta,]
df.sel.prot.body <- df.prot.body[ind.prot.meta,]
df.sel.gene.body <- df.gene.body[ind.gene.meta,]

df.all.body <- cbind(df.sel.prot.body, df.sel.gene.body)

# Write to CSV
write.csv(df.all.meta, paste("../../data/ex_7/gene_prot_data_meta.csv", sep=""), row.names=TRUE)
write.csv(df.all.body, paste("../../data/ex_7/gene_prot_data_body.csv", sep=""), row.names=TRUE)
write.csv(df.sel.prot.body, paste("../../data/ex_7/prot_data_body.csv", sep=""), row.names=TRUE)
write.csv(df.sel.gene.body, paste("../../data/ex_7/gene_data_body.csv", sep=""), row.names=TRUE)
