setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_7/")

df.all.meta <- read.csv("../../data/ex_7/gene_prot_data_meta.csv", header=TRUE, row.names = 1)
df.all.meta$group <- as.character(df.all.meta$group)

df.gene.body <- read.csv("../../data/ex_7/gene_data_body_lt2_qn.csv", header=TRUE, row.names = 1) 

df.prot.body <- read.csv("../../data/ex_7/prot_data_body_lt2_qn.csv", header=TRUE, row.names = 1)

datasets <- list(
  list(df.gene.body, "gene", "gene")
  ,
  list(df.prot.body,  "protein", "prot")
)

for (i in 1:length(datasets)){
  data <- datasets[[i]][[1]]
  verbose <- datasets[[i]][[2]]
  abbrv <- datasets[[i]][[3]]
  
  # Select HIV- patients
  #ind.hiv_neg <- c()
  
  #for (j in 1:nrow(data)){
    #if ((df.all.meta$group[j] == 1) || (df.all.meta$group[j] == 3) || (df.all.meta$group[j] == 6)){
      #ind.hiv_neg <- c(ind.hiv_neg, j)
    #}
  #}
  
  #data.hiv_neg <- data[ind.hiv_neg,]
  #meta.hiv_neg <- df.all.meta[ind.hiv_neg,]
  
  # Total number of patients in validation set
  tot.num_20 <- 0
  
  # Row numbers of patients in validation set
  tot.ind.validation <- c()
  
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
  
  df.body.validation <- data[tot.ind.validation,]
  df.body.test_train <- data[-tot.ind.validation,]
  
  df.meta.validation <- df.all.meta[tot.ind.validation,]
  df.meta.test_train <-df.all.meta[-tot.ind.validation,]
  
  write.csv(df.body.validation, paste("../../data/ex_7/", abbrv, "_validation_body.csv", sep=""),row.names=TRUE)
  write.csv(df.body.test_train, paste("../../data/ex_7/", abbrv, "_testtrain_body.csv", sep=""),row.names=TRUE)
  
  write.csv(df.meta.validation, paste("../../data/ex_7/gene_prot_validation_meta.csv", sep=""),row.names=TRUE)
  write.csv(df.meta.test_train, paste("../../data/ex_7/gene_prot_testtrain_meta.csv", sep=""),row.names=TRUE)
}
