setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_12/")

set.seed(12)

print("Loading metadata")
df.all.meta <- read.csv("../../data/ex_12/gp_data_meta.csv", header=TRUE, row.names = 1)
df.all.meta$group <- as.character(df.all.meta$group)

print("Loading gene data - this will take a while")
df.gene.body <- read.csv("../../data/ex_12/gene_data_body_l2t_qn.csv", header=TRUE, row.names = 1) 

print("Loading protein data")
df.prot.body <- read.csv("../../data/ex_12/prot_data_body_l2t_qn.csv", header=TRUE, row.names = 1)

# Total number of patients in validation set
tot.num_20 <- 0

# Row numbers of patients in validation set
tot.ind.validation <- c()

#Select random indices of patients by disese category

print("Randomly sampling patients from each disease category")

for (i in 1:6){
  
  group_len <- length(df.all.meta$group[df.all.meta$group == i])
  num_20 <- round(0.2*group_len)
  tot.num_20 <- tot.num_20 + num_20
  
  #Get row numbers of patients in disease group
  ind.group <- c()
  
  for (j in 1:nrow(df.all.meta)){
    if (df.all.meta$group[j] == i){
      ind.group <- c(ind.group, j)
    }
  }
  
  #Get validation split of group
  ind.validation <- sample(ind.group, num_20, replace=FALSE)
  tot.ind.validation <- c(tot.ind.validation, ind.validation)
  
}

#Take validation and test/train splits of all data

print("Generating validation and training sets")

df.gene.body.validation <- df.gene.body[tot.ind.validation,]
df.gene.body.test_train <- df.gene.body[-tot.ind.validation,]

df.prot.body.validation <- df.prot.body[tot.ind.validation,]
df.prot.body.test_train <- df.prot.body[-tot.ind.validation,]

df.meta.validation <- df.all.meta[tot.ind.validation,]
df.meta.test_train <-df.all.meta[-tot.ind.validation,]

# Create CSVs

print("Writing CSVs")

write.csv(df.gene.body.validation, paste("../../data/ex_12/gene_validation_body.csv", sep=""),row.names=TRUE)
write.csv(df.gene.body.test_train, paste("../../data/ex_12/gene_train_body.csv", sep=""),row.names=TRUE)

write.csv(df.prot.body.validation, "../../data/ex_12/prot_validation_body.csv",row.names=TRUE)
write.csv(df.prot.body.test_train, "../../data/ex_12/prot_train_body.csv",row.names=TRUE)

write.csv(df.meta.validation, "../../data/ex_12/gp_validation_meta.csv", row.names=TRUE)
write.csv(df.meta.test_train, "../../data/ex_12/gp_train_meta.csv", row.names=TRUE)

