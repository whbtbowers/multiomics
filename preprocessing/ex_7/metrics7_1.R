setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_7/")

library(ggplot2)
library(stats)
library(ggfortify)
library(heatmap3)

# Meta data and body data for genes and proteins
df.all.meta <- read.csv("../../data/ex_7/gene_prot_data_meta.csv", header=TRUE, row.names = 1)  # Meta data
df.all.body <- read.csv("../../data/ex_7/gene_prot_data_body.csv", header=TRUE, row.names = 1)  # Initial
df.all.body.lt <- read.csv("../../data/ex_7/gene_prot_data_body_lt.csv", header=TRUE, row.names = 1)  # Log transformed
df.all.body.wn <- read.csv("../../data/ex_7/gene_prot_data_body_wn.csv", header=TRUE, row.names = 1)  # Normalised
df.all.body.lt.wn <- read.csv("../../data/ex_7/gene_prot_lt_data_body_wn.csv", header=TRUE, row.names = 1)  # Log transformed THEN normalised
df.all.body.wn.lt <- read.csv("../../data/ex_7/gene_prot_wn_data_body_lt.csv", header=TRUE, row.names = 1)  # Normalised THEN log transformed

df.all.meta$group <- as.character(df.all.meta$group)

# Just genes
df.gene.body <- read.csv("../../data/ex_7/gene_data_body.csv", header=TRUE, row.names = 1)  # Initial
df.gene.body.lt <- read.csv("../../data/ex_7/gene_data_body_lt.csv", header=TRUE, row.names = 1)  # Log transformed
df.gene.body.wn <- read.csv("../../data/ex_7/gene_data_body_wn.csv", header=TRUE, row.names = 1)  # Normalised
df.gene.body.lt.wn <- read.csv("../../data/ex_7/gene_lt_data_body_wn.csv", header=TRUE, row.names = 1)  # Log transformed THEN normalised
df.gene.body.wn.lt <- read.csv("../../data/ex_7/gene_wn_data_body_lt.csv", header=TRUE, row.names = 1)  # Normalised THEN log transformed

# Just proteins
df.prot.body <- read.csv("../../data/ex_7/prot_data_body.csv", header=TRUE, row.names = 1)  # Initial
df.prot.body.lt <- read.csv("../../data/ex_7/prot_data_body_lt.csv", header=TRUE, row.names = 1)  # Log transformed
df.prot.body.wn <- read.csv("../../data/ex_7/prot_data_body_wn.csv", header=TRUE, row.names = 1)  # Normalised
df.prot.body.lt.wn <- read.csv("../../data/ex_7/prot_lt_data_body_wn.csv", header=TRUE, row.names = 1)  # Log transformed THEN normalised
df.prot.body.wn.lt <- read.csv("../../data/ex_7/prot_wn_data_body_lt.csv", header=TRUE, row.names = 1)  # Normalised THEN log transformed

# Set date and labels
date <- "2018-07-17"

datasets <- list(
  #list(df.gene.body, "Initial gene", "gene")
  #,
  #list(df.prot.body, "Initial protein", "prot")
  #,
  #list(df.gene.body.lt, "Log transformed gene", "gene_lt")
  #,
  #list(df.prot.body.lt, "Log transformed protein", "prot_lt")
  #,
  #list(df.gene.body.wn, "Normalised gene", "gene_wn")
  #,
  #list(df.prot.body.wn, "Normalised protein", "prot_wn")
  #,
  #list(df.gene.body.lt.wn, "Log transformed then normalised gene", "gene_lt_wn")
  #,
  #list(df.prot.body.lt.wn, "Log transformed then normalised protein", "prot_lt_wn")
  #,
  #list(df.gene.body.wn.lt, "Normalised then log transformed gene", "gene_wn_lt")
  #,
  list(df.prot.body.wn.lt, "Normalised then log transformed protein", "prot_wn_lt")
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

    
  #Boxplot
  png(paste("../../img/ex_7/", date, "/",   abbrv, "_bplot.png", sep=""),
      width = 10000,        
      height = 2500,
      res = 300,          
      pointsize = 12)        
  
  plot.new()
  par(xaxt="n", mar=c(10,5,3,1))
  boxplot(data.hiv_neg, col="gold")
  lablist<-as.vector(colnames(data.hiv_neg))
  axis(1, at=seq(1, ncol(data.hiv_neg), by=1), labels = FALSE)
  text(seq(1, ncol(data.hiv_neg), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
  title(paste( verbose, "data", sep=" "))
  dev.off()
  
  # Heatmap
  png(paste("../../img/ex_7/", date, "/",  abbrv, "_hmap.png", sep=""),
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  
  heatmap3(data.hiv_neg,
           #Rowv=NA,
           #Colv=NA,
           margins=c(10,5),
           main= paste(verbose, "data", sep=" "),
           xlab= "prots",
           ylab= "Patients",
  )
  dev.off()
  
  # Principal components 1 & 2 by location
  
  png(paste("../../img/ex_7/", date, "/",  abbrv, "_pca_pc1_pc2_bysite.png", sep=""),
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
  
  png(paste("../../img/ex_7/", date, "/",  abbrv, "_pca_pc3_pc4_bysite.png", sep=""),
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
  
  png(paste("../../img/ex_7/", date, "/",  abbrv, "_pca_pc1_pc2_byinf.png", sep=""),
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
  
  png(paste("../../img/ex_7/", date, "/",  abbrv, "_pca_pc3_pc4_byinf.png", sep=""),
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


