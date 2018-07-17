setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_6/")

library(ggplot2)
library(stats)
library(ggfortify)
library(heatmap3)


#df.prot.data <- read.csv("../../data/ex_6/prot_data_body_filtimp_k20.csv", header=TRUE, row.names = 1) # Corrected and log transformed protein data
df.prot.data <- read.csv("../../data/ex_6/prot_data_body_wnormed.csv", header=TRUE, row.names = 1) # Wang normalised protein data
df.prot.meta <- read.csv("../../data/ex_6/prot_data_meta.csv", header=TRUE, row.names = 1) # Protein meta data

df.prot.meta$group <- as.character(df.prot.meta$group)

#df.gene.data <- read.csv("../../data/ex_6/gene_data_body.csv", header=TRUE, row.names = 1) # Corrected and log transformed protein data
df.gene.data <- read.csv("../../data/ex_6/gene_data_body_wnormed.csv", header=TRUE, row.names = 1) # Wang normalised gene data
df.gene.meta <- read.csv("../../data/ex_6/gene_data_meta.csv", header=TRUE, row.names = 1) # Gene meta data

df.gene.meta$group <- as.character(df.gene.meta$group)

# Set date and labels
label_verbose <- "Normalised"
label_abbrv <- "wnormed_"
date <- "2018-07-13"

datasets <- list(
  list(df.gene.data, df.gene.meta, "gene", "gene")
  ,
  list(df.prot.data, df.prot.meta, "protein", "prot")
)

for (i in 1:length(datasets)){
  data <- datasets[[i]][[1]]
  meta <- datasets[[i]][[2]]
  verbose <- datasets[[i]][[3]]
  abbrv <- datasets[[i]][[4]]
  
  # Select HIV- 
  ind.hiv_neg <- c()
  
  for (j in 1:nrow(data)){
    if ((meta$group[j] == 1) || (meta$group[j] == 3) || (meta$group[j] == 6)){
      ind.hiv_neg <- c(ind.hiv_neg, j)
    }
  }
  
  
  data.hiv_neg <- data[ind.hiv_neg,]
  meta.hiv_neg <- meta[ind.hiv_neg,]

    
  #Boxplot
  png(paste("../../img/ex_6/", date, "/", label_abbrv,  abbrv, "_bplot.png", sep=""),
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
  title(paste(label_verbose, verbose, "data", sep=" "))
  dev.off()
  
  # Heatmap
  png(paste("../../img/ex_6/", date, "/", label_abbrv, abbrv, "_hmap.png", sep=""),
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  
  heatmap3(data.hiv_neg,
           #Rowv=NA,
           #Colv=NA,
           margins=c(10,5),
           main= paste(label_verbose, verbose, "data", sep=" "),
           xlab= "prots",
           ylab= "Patients",
  )
  dev.off()
  
  # Principal components 1 & 2 by infection status
  
  png(paste("../../img/ex_6/", date, "/", label_abbrv, abbrv, "_pca_pc1_pc2_byinf.png", sep=""),
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  
  
  autoplot(prcomp(data.hiv_neg),
           data =meta.hiv_neg,
           colour = "group",
           main=paste(label_verbose, verbose, "data", sep=" ")
  )
  
  dev.off()
  
  # Principal components 3 & 4 by infection status
  
  png(paste("../../img/ex_6/", date, "/", label_abbrv, abbrv, "_pca_pc3_pc4_byinf.png", sep=""),
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  
  autoplot(prcomp(data.hiv_neg),
           x=3, 
           y=4,
           data =meta.hiv_neg,
           colour = "group",
           main=paste(label_verbose, verbose, "data", sep=" ")
  )
  
  dev.off()
  

}


