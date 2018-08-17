setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_8/")

library(ggplot2)
library(stats)
library(ggfortify)

# Meta data and body data for genes and proteins
df.meta <- read.csv("../../data/ex_9/gp_data_meta.csv", header=TRUE, row.names = 1)  # Meta data
df.meta$group <- as.character(df.meta$group)

# Just genes
df.gene.body <- read.csv("../../data/ex_8/gene_data_body.csv", header=TRUE, row.names = 1)  # Initial
df.gene.body.lt2 <- read.csv("../../data/ex_9/gene_data_body_l2t.csv", header=TRUE, row.names = 1)  # Log transformed
df.gene.body.lt2.qn <- read.csv("../../data/ex_9/gene_data_body_l2t_qn.csv", header=TRUE, row.names = 1)  # Log2 transformed THEN quantile normalised

# Just proteins
df.prot.body <- read.csv("../../data/ex_8/prot_data_body.csv", header=TRUE, row.names = 1)  # Initial
df.prot.body.lt2 <- read.csv("../../data/ex_9/prot_data_body_l2t.csv", header=TRUE, row.names = 1)  # Log2 transformed
df.prot.body.lt2.qn <- read.csv("../../data/ex_9/prot_data_body_l2t_qn.csv", header=TRUE, row.names = 1)  # Log2 transformed THEN quantile normalised

# Set date and labels
date <- "2018-07-31/"
ex_dir <- "ex_9/"

datasets <- list(
  list(df.gene.body, "Initial gene", "gene")
  #,
  #list(df.prot.body, "Initial protein", "prot")
  #,
  #list(df.gene.body.lt2, "Log2 transformed gene", "gene_l2t")
  #,
  #list(df.prot.body.lt2, "Log2 transformed protein", "prot_l2t")
  #,
  #list(df.gene.body.lt2.qn, "Log2 transformed then quantile normalised gene", "gene_l2t_qn")
  #,
  #list(df.prot.body.lt2.qn, "Log2 transformed then quantile normalised protein", "prot_l2t_qn")
  )

for (i in 1:length(datasets)){
  data <- datasets[[i]][[1]]
  verbose <- datasets[[i]][[2]]
  abbrv <- datasets[[i]][[3]]
  
  # Select HIV- 
  ind.hiv_neg <- c()
  
  for (j in 1:nrow(data)){
    if ((df.meta$group[j] == 1) || (df.meta$group[j] == 3) || (df.meta$group[j] == 6)){
      ind.hiv_neg <- c(ind.hiv_neg, j)
    }
  }
  
  data.hiv_neg <- data[ind.hiv_neg,]
  meta.hiv_neg <- df.meta[ind.hiv_neg,]
  
  #Select probe and patient labels
  lab.probe <- colnames(data.hiv_neg)
  lab.patient <- rownames(data.hiv_neg)
 
  
  # Principal components 1 & 2 by location
  
  png(paste("../../img/", ex_dir, date,   abbrv, "_pca_pc1_pc2_bysite.png", sep=""),
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
  
  png(paste("../../img/", ex_dir, date,   abbrv, "_pca_pc3_pc4_bysite.png", sep=""),
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
  
  png(paste("../../img/", ex_dir, date,   abbrv, "_pca_pc1_pc2_byinf.png", sep=""),
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
  
  png(paste("../../img/", ex_dir, date,   abbrv, "_pca_pc3_pc4_byinf.png", sep=""),
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
  
  # Transform data so boxplot gives distribution of each patient
  
  data.t <- as.data.frame(t(data))
  colnames(data.t) <- lab.patient
  rownames(data.t) <- lab.probe
  
  
  #Boxplot
  png(paste("../../img/", ex_dir, date, abbrv, "_bplot.png", sep=""),
      width = 10000,        
      height = 2500,
      res = 300,          
      pointsize = 12)        
  
  plot.new()
  par(xaxt="n", mar=c(10,5,3,1))
  boxplot(data.t, col="red")
  lablist<-as.vector(colnames(data.t))
  axis(1, at=seq(1, ncol(data.t), by=1), labels = FALSE)
  text(seq(1, ncol(data.t), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
  title(paste(verbose, "data", sep=" "))
  dev.off()
  
}


