setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_7/")

library(ggplot2)
library(stats)
library(ggfortify)
library(heatmap3)

df.prot.meta <- read.csv("../../data/ex_7/gene_prot_data_meta.csv", header=TRUE, row.names = 1)
df.prot.body <- read.csv("../../data/ex_7/prot_data_body.csv", header=TRUE, row.names = 1)
df.prot.body.lt2 <- read.csv("../../data/ex_7/prot_data_body_lt2.csv", header=TRUE, row.names = 1)
df.prot.body.lt2.qn <- read.csv("../../data/ex_7/prot_data_body_lt2_qn.csv", header=TRUE, row.names = 1)

# To direct to the correct folder
date <- "2018-07-25/"
ex_dir <- "ex_7/"

#Check for presence of GBP6
for i in 

# Sets to run through automated graphing

datasets <- list(
  list(df.prot.body, "initial", "prot")
  ,
  list(df.prot.body.lt2, "log2 transformed", "prot_lt2")
  ,
  list(df.prot.body.lt2.qn, "log2 transformed and quantile normalised", "prot_lt2_qn")
)

for (d in 1:length(datasets)){
  ds.data <- datasets[[d]][[1]]
  ds.ttl_adj <- datasets[[d]][[2]]
  ds.abrv <- datasets[[d]][[3]]
  
  #Select HIV- patients
  
  ind.hiv_neg <- c()
  
  for (i in 1:nrow(df.prot.body)){
    if((df.prot.meta$group[i] == 1) || (df.prot.meta$group[i] == 3) || (df.prot.meta$group[i] == 6)){
      ind.hiv_neg <- c(ind.hiv_neg, i)
    }
  }
  
  ds.data.hiv_neg <- ds.data[ind.hiv_neg,]
  df.prot.meta.hiv_neg <- df.prot.meta[ind.hiv_neg,]
  
  # Create melt table across seperate data and metadata
  Label <- c()
  variable <- c()
  value <- c()
  
  for (i in 1:nrow(ds.data.hiv_neg)){
    for (j in 1:ncol(ds.data.hiv_neg)){
      Label <- c(Label, as.character(df.prot.meta.hiv_neg$tb.status[i]))
      variable <- c(variable, colnames(ds.data)[j])
      value <- c(value, ds.data[i,j])
    }
  }
  
  prot.hiv_neg.melt <- data.frame(Label, variable, value)
  
  #Boxplot

  plot.melt <- ggplot(data = prot.hiv_neg.melt, 
                      aes(x=variable, y=value)) + 
    geom_boxplot(aes(fill=Label)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title=paste("Distribution of", ds.ttl_adj, "protein expression values", sep=" "),
         x ="Proteins", 
         y = "Expression"
    )
  
  ggsave(paste("../../img/", ex_dir, date, "/", ds.abrv, "_bplot.png", sep=""),
         plot = plot.melt,
         device = "png",
         width = 10,        
         height = 5,
         dpi = 300,
         limitsize = FALSE
  )
  
}