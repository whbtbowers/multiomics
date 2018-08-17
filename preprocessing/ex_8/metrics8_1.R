setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_8/")

library(ggplot2)
library(stats)
library(ggfortify)
library(heatmap3)

# Meta data and body data for genes and proteins
df.all.meta <- read.csv("../../data/ex_8/gp_data_meta.csv", header=TRUE, row.names = 1)  # Meta data
df.all.meta$group <- as.character(df.all.meta$group)

# Just genes
df.gene.body <- read.csv("../../data/ex_8/gene_data_body.csv", header=TRUE, row.names = 1)  # Initial
df.gene.body.lt2 <- read.csv("../../data/ex_8/gene_data_body_lt2.csv", header=TRUE, row.names = 1)  # Log transformed
df.gene.body.lt2.qn <- read.csv("../../data/ex_8/gene_data_body_lt2_qn.csv", header=TRUE, row.names = 1)  # Log2 transformed THEN quantile normalised

# Just proteins
df.prot.body <- read.csv("../../data/ex_8/prot_data_body.csv", header=TRUE, row.names = 1)  # Initial
df.prot.body.lt2 <- read.csv("../../data/ex_8/prot_data_body_lt2.csv", header=TRUE, row.names = 1)  # Log2 transformed
df.prot.body.lt2.qn <- df.prot.body.lt.wn <- read.csv("../../data/ex_8/prot_data_body_lt2_qn.csv", header=TRUE, row.names = 1)  # Log2 transformed THEN quantile normalised

# Set date and labels
date <- "2018-07-30"
ex_dir <- "ex_8/"

datasets <- list(
  #list(df.gene.body, "Initial gene", "gene")
  #,
  #list(df.prot.body, "Initial protein", "prot")
  #,
  #list(df.gene.body.lt2, "Log2 transformed gene", "gene_lt2")
  #,
  #list(df.prot.body.lt2, "Log2 transformed protein", "prot_lt2")
  #,
  #list(df.gene.body.lt2.qn, "Log2 transformed then quantile normalised gene", "gene_lt2_qn")
  #,
  list(df.prot.body.lt2.qn, "Log2 transformed then quantile normalised protein", "prot_lt2_qn")
  )

for (i in 1:length(datasets)){
  data <- datasets[[i]][[1]]
  verbose <- datasets[[i]][[2]]
  abbrv <- datasets[[i]][[3]]
  
  # Select HIV- 
  ind.hiv_neg <- c()
  
  for (j in 1:nrow(data)){
    if ((df.all.meta$group[j] == 1) || (df.all.meta$group[j] == 3) || (df.all.meta$group[j] == 6)){
      ind.hiv_neg <- c(ind.hiv_neg, j)
    }
  }
  
  data.hiv_neg <- data[ind.hiv_neg,]
  meta.hiv_neg <- df.all.meta[ind.hiv_neg,]
  
  # Principal components 1 & 2 by location
  
  png(paste("../../img/", ex_dir, date, "/",  abbrv, "_pca_pc1_pc2_bysite.png", sep=""),
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  
  
  autoplot(prcomp(data.hiv_neg),
           data =meta.hiv_neg,
           colour = "site",
           main=paste( verbose, "data", sep=" ")
  )
  
  dev.off()
  
  # Principal components 3 & 4 by by location
  
  png(paste("../../img/", ex_dir, date, "/",  abbrv, "_pca_pc3_pc4_bysite.png", sep=""),
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  
  autoplot(prcomp(data.hiv_neg),
           x=3, 
           y=4,
           data =meta.hiv_neg,
           colour = "site",
           main=paste( verbose, "data", sep=" ")
  )
  
  dev.off()
  
  # Principal components 1 & 2 by infection status
  
  png(paste("../../img/", ex_dir, date, "/",  abbrv, "_pca_pc1_pc2_byinf.png", sep=""),
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  
  
  autoplot(prcomp(data.hiv_neg),
           data =meta.hiv_neg,
           colour = "group",
           main=paste( verbose, "data", sep=" ")
  )
  
  dev.off()
  
  # Principal components 3 & 4 by infection status
  
  png(paste("../../img/", ex_dir, date, "/",  abbrv, "_pca_pc3_pc4_byinf.png", sep=""),
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  
  autoplot(prcomp(data.hiv_neg),
           x=3, 
           y=4,
           data =meta.hiv_neg,
           colour = "group",
           main=paste( verbose, "data", sep=" ")
  )
  
  dev.off()
  
}


