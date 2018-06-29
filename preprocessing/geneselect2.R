setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/")
library(stringr)

df.gene_data <- read.csv("../data/gene_data.csv")
#df.gene_data <- read.csv("../data/gene_data_sub100.csv", header=TRUE, row.names = 1)
df.gene_inf <- read.csv("../data/gene_inf.csv", header=TRUE)

for (i in rev(1:nrow(df.gene_inf))){
  if (df.gene_inf$array.id[i] == "no"){
    df.gene_inf <- df.gene_inf[-i,]
  } else if (is.na(df.gene_inf$p1.sample.id[i])){
    df.gene_inf <- df.gene_inf[-i,]
  }
}

#df.gene_data.subset <- df.gene_data[1:100,1:100]
#write.csv(df.gene_data.subset, file="../data/gene_data_sub100.csv")

#df.gene_inf$inf.status <- toString(df.gene_inf$inf.status)

# Call indices lists for gene info

ind.gene_inf.hiv_neg <- c()
ind.gene_inf.hiv_pos <- c()

for (i in 1:nrow(df.gene_inf)){
  if (str_detect(df.gene_inf$inf.status[i],"HIV-")){
    ind.gene_inf.hiv_neg <- c(ind.gene_inf.hiv_neg, i)
  } else if(str_detect(df.gene_inf$inf.status[i],"HIV+")){
    ind.gene_inf.hiv_pos <- c(ind.gene_inf.hiv_pos, i)
  }
}

# Subset gene inf
df.gene_inf.hiv_pos <- df.gene_inf[ind.gene_inf.hiv_pos,]
df.gene_inf.hiv_neg <- df.gene_inf[ind.gene_inf.hiv_neg,]

#Get indices of patients in gene_inf
inf_ind.gene_inf.hiv_neg <- c()
inf_ind.gene_inf.hiv_pos <- c()

#Get indices of patients in big dataset
ind.gene_data.hiv_neg <- c()
ind.gene_data.hiv_pos <- c()

# Grab HIV- patients from big dataset based on array id
for (i in 1:nrow(df.gene_inf.hiv_neg)){
  for (j in 1:nrow(df.gene_data)){
    if (str_detect(toString(df.gene_inf.hiv_neg$array.id[i]), toString(df.gene_data$array.id[j]))){
      inf_ind.gene_inf.hiv_neg <- c(inf_ind.gene_inf.hiv_neg, i)
      ind.gene_data.hiv_neg <- c(ind.gene_data.hiv_neg, j) 
    }
  }
}

# Grab HIV+ patients from big dataset based on array id
for (i in 1:nrow(df.gene_inf.hiv_pos)){
  for (j in 1:nrow(df.gene_data)){
    if (str_detect(toString(df.gene_inf.hiv_pos$array.id[i]), toString(df.gene_data$array.id[j]))){
      inf_ind.gene_inf.hiv_pos <- c(inf_ind.gene_inf.hiv_pos, i)
      ind.gene_data.hiv_pos <- c(ind.gene_data.hiv_pos, j) 
    }
  }
}

#Create disease status column
df.inf_status.hiv_neg <- data.frame(df.gene_inf.hiv_neg$inf.status[unique(inf_ind.gene_inf.hiv_neg)])
df.inf_status.hiv_pos <- data.frame(df.gene_inf.hiv_pos$inf.status[unique(inf_ind.gene_inf.hiv_pos)])

#combine into single HIV- data frame
df.gene_data.hiv_neg <- df.gene_data[unique(ind.gene_data.hiv_neg),]
array_id.hiv_neg <-df.gene_data.hiv_neg$array.id
df.hiv_neg <- data.frame(array_id.hiv_neg, df.inf_status.hiv_neg, df.gene_data.hiv_neg)
df.hiv_neg <- df.hiv_neg[,-3]
colnames(df.hiv_neg)[1] <- "array_id"
colnames(df.hiv_neg)[2] <- "inf.status"

#combine into single HIV+ data frame
df.gene_data.hiv_pos <- df.gene_data[unique(ind.gene_data.hiv_pos),]
array_id.hiv_pos <-df.gene_data.hiv_pos$array.id
df.hiv_pos <- data.frame(array_id.hiv_pos, df.inf_status.hiv_pos, df.gene_data.hiv_pos)
df.hiv_pos <- df.hiv_pos[,-3]
colnames(df.hiv_pos)[1] <- "array.id"
colnames(df.hiv_pos)[2] <- "inf.status"


# Write CSVs

print("Writing CSVs")
write.csv(df.hiv_pos, file="../data/gene_ex/gene_data_hiv_pos.csv")
write.csv(df.hiv_neg, file="../data/gene_ex/gene_data_hiv_neg.csv")