library(stringr)

df.gene_data <- read.csv("../data/gene_data.csv")
#df.gene_data <- read.csv("../data/gene_data_sub100.csv", header=TRUE, row.names = 1)
df.gene_inf <- read.csv("../data/gene_inf.csv", header=TRUE)

#df.gene_data.subset <- df.gene_data[1:100,1:100]
#write.csv(df.gene_data.subset, file="../data/gene_data_sub100.csv")

#df.gene_inf$inf.status <- toString(df.gene_inf$inf.status)

# Call indices lists for gene info

ind.gene.hiv_pos <- c()
ind.gene.hiv_neg <- c()

for (i in 1:nrow(df.gene_inf)){
  if (str_detect(df.gene_inf$inf.status[i], "HIV")){
    if (str_detect(df.gene_inf$inf.status[i], "HIV-")){
      ind.gene.hiv_neg <- c(ind.gene.hiv_neg, i)
    }else{
      ind.gene.hiv_pos <- c(ind.gene.hiv_pos, i)
    }
  }
}

# List of all indices
ind.gene.hiv_all <- c(ind.gene.hiv_pos, ind.gene.hiv_neg)

# Subset gene inf
df.gene_inf.hiv_all <- df.gene_inf[ind.gene.hiv_all,]
df.gene_inf.hiv_pos <- df.gene_inf[ind.gene.hiv_pos,]
df.gene_inf.hiv_neg <- df.gene_inf[ind.gene.hiv_neg,]

#Get indices of big dataset
ind.gene_data.hiv_all <- c()
ind.gene_data.hiv_pos <- c()
ind.gene_data.hiv_neg <- c()

# Select gene data indices from big dataset based on array ID

print("Finding indices of all patients in big dataset")

for (i in 1:nrow(df.gene_inf.hiv_all)){
  for (j in 1:nrow(df.gene_data)){
    if (str_detect(toString(df.gene_inf.hiv_all$array.id[i]), toString(df.gene_data$array.id[j]))){
      ind.gene_data.hiv_all <- c(ind.gene_data.hiv_all, j)
    }
  }
}

ind.gene_data.hiv_all <- unique(ind.gene_data.hiv_all)

# Select HIV+ gene indices data from big dataset based on array ID

print("Finding indices of HIV+ patients in big dataset")

for (i in 1:nrow(df.gene_inf.hiv_pos)){
  for (j in 1:nrow(df.gene_data)){
    if (str_detect(toString(df.gene_inf.hiv_pos$array.id[i]), toString(df.gene_data$array.id[j]))){
      ind.gene_data.hiv_pos <- c(ind.gene_data.hiv_pos, j)
    }
  }
}

ind.gene_data.hiv_pos <- unique(ind.gene_data.hiv_pos)

# Select HIV- gene indices data from big dataset based on array ID

print("Finding indices of HIV- patients in big dataset")

for (i in 1:nrow(df.gene_inf.hiv_neg)){
  for (j in 1:nrow(df.gene_data)){
    if (str_detect(toString(df.gene_inf.hiv_neg$array.id[i]), toString(df.gene_data$array.id[j]))){
      ind.gene_data.hiv_neg <- c(ind.gene_data.hiv_neg, j)
    }
  }
}

ind.gene_data.hiv_neg <- unique(ind.gene_data.hiv_neg)

# Use indices to select rows from big dataset

print("Subsetting big dataset")

df.gene_data.hiv_all <- df.gene_data[ind.gene_data.hiv_all,]
df.gene_data.hiv_pos <- df.gene_data[ind.gene_data.hiv_pos,]
df.gene_data.hiv_neg <- df.gene_data[ind.gene_data.hiv_neg,]

# Write CSVs

print("Writing CSVs")

write.csv(df.gene_data.hiv_all, file="../data/gene_ex/gene_data_hiv_all.csv")
write.csv(df.gene_data.hiv_pos, file="../data/gene_ex/gene_data_hiv_pos.csv")
write.csv(df.gene_data.hiv_neg, file="../data/gene_ex/gene_data_hiv_neg.csv")