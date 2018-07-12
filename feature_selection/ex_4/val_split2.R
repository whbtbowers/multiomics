setwd("/project/home17/whb17/Documents/project3/project_files/feature_selection/ex_4/")

df.gene.body <- read.csv("../../data/ex_4/sel_gene_data_body.csv", header=TRUE, row.names = 1) # Normalised like Wang, labelled, 32 selected by elastic net
df.gene.meta <- read.csv("../../data/ex_4/gene_data_meta.csv", header=TRUE, row.names = 1)
df.gene.meta$group <- as.character(df.gene.meta$group)

df.prot.body <- read.csv("../../data/ex_4/prot_data_body.csv", header=TRUE, row.names = 1) # Normalised like Wang, labelled, 32 selected by elastic net
df.prot.meta <- read.csv("../../data/ex_4/prot_data_meta.csv", header=TRUE, row.names = 1)
df.prot.meta$group <- as.character(df.prot.meta$group)

datasets <- list(
  list(df.gene.body, df.gene.meta, "gene", "gene")
  ,
  list(df.prot.body, df.prot.meta, "protein", "prot")
)

for (i in 1:length(datasets)){
  data <- datasets[[i]][[1]]
  meta <- datasets[[i]][[2]]
  verbose <- datasets[[i]][[3]]
  abbrv <- datasets[[i]][[4]]
  
  # Select HIV- patients
  ind.hiv_neg <- c()
  
  for (j in 1:nrow(data)){
    if ((meta$group[j] == 1) || (meta$group[j] == 3) || (meta$group[j] == 6)){
      ind.hiv_neg <- c(ind.hiv_neg, j)
    }
  }
  
  data.hiv_neg <- data[ind.hiv_neg,]
  meta.hiv_neg <- meta[ind.hiv_neg,]
  
  # Total number of patients in validation set
  tot.num_30 <- 0
  
  # Row numbers of patients in validation set
  tot.ind.validation <- c()
  
  for (i in 1:6){
    
    group_no <- 1
    group_len <- length(df.gene.meta$group[df.gene.meta$group == group_no])
    num_30 <- round(0.3*group_len)
    tot.num_30 <- tot.num_30 + num_30
    
    #Get row numbers of patients in disease group
    ind.group <- c()
    
    for (j in 1:nrow(df.gene.meta)){
      if (df.gene.meta$group[j] == group_no){
        ind.group <- c(ind.group, j)
      }
    }
    
    #Get validation split of group
    ind.validation <- sample(ind.group, num_30, replace=FALSE)
    tot.ind.validation <- c(tot.ind.validation, ind.validation)
    
  }
  
  #Take validation and test/train splits of all data
  
  df.body.validation <- df.gene.body[tot.ind.validation,]
  df.body.test_train <- df.gene.body[-tot.ind.validation,]
  
  df.meta.validation <- df.gene.meta[tot.ind.validation,]
  df.meta.test_train <- df.gene.meta[-tot.ind.validation,]
  
  write.csv(df.body.validation, paste("../../data/ex_4/", abbrv, "_validation_body.csv"),row.names=TRUE)
  write.csv(df.meta.test_train, paste("../../data/ex_4/", abbrv, "_testtrain_body.csv"),row.names=TRUE)
  
  write.csv(df.meta.validation, paste("../../data/ex_4/", abbrv, "_validation_meta.csv"),row.names=TRUE)
  write.csv(df.meta.test_train, paste("../../data/ex_4/", abbrv, "_testtrain_meta.csv"),row.names=TRUE)
}
