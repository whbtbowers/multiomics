setwd("/project/home17/whb17/Documents/project3/project_files/feature_selection/g_ex/")

df.gene_data <- read.csv("../../data/gene_ex2/labd_24sel_gene_data_ex_scaled.csv", header=TRUE, row.names = 1) # Scaled, labelled, 24 selected by elastic net
df.gene_data$group <- as.character(df.gene_data$group)

# Total number of patients in valudation set
tot.num_20 <- 0

# Row numbers of intividuals in validation set
tot.ind.validation <- c()

for (i in 1:6){
  
  group_no <- i
  group_len <- length(df.gene_data$group[df.gene_data$group==group_no])
  num_20 <- round(0.2*group_len)
  tot.num_20 <- tot.num_20 + num_20
  
  #Get row numbers of patients in disease group
  ind.group <- c()
  
  for (i in 1:nrow(df.gene_data)){
    if (df.gene_data$group[i] == group_no){
      ind.group <- c(ind.group, i)
    }
  }
  
  #Get validation split of group
  ind.validation <- sample(ind.group, num_20, replace=FALSE)
  tot.ind.validation <- c(tot.ind.validation, ind.validation)
}
#Sample group 1
group_no <- 1
group_len <- length(df.gene_data$group[df.gene_data$group==group_no])
num_20 <- round(0.2*group_len)

# Get row number of patients in disease group 
ind.group <- c()

for (i in 1:nrow(df.gene_data)){
  if (df.gene_data$group[i] == group_no){
    ind.group <- c(ind.group, i)
  }
}

#Take validation split of all data
df.validation <- df.gene_data[tot.ind.validation,]
df.test_train <- df.gene_data[-tot.ind.validation,]

write.csv(df.validation,"../../data/gene_ex2/gene_validation.csv",row.names=TRUE)
write.csv(df.test_train,"../../data/gene_ex2/gene_test_train.csv",row.names=TRUE)

