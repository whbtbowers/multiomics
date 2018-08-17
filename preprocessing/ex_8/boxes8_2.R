setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_8/")

library(ggplot2)
library(stats)
library(ggfortify)
library(heatmap3)

df.gene.meta <- read.csv("../../data/ex_8/gp_data_meta.csv", header=TRUE, row.names = 1)
df.gene.body <- read.csv("../../data/ex_8/gene_data_body.csv", header=TRUE, row.names = 1)
df.gene.body.lt2 <- read.csv("../../data/ex_8/gene_data_body_lt2.csv", header=TRUE, row.names = 1)
df.gene.body.lt2.qn <- read.csv("../../data/ex_8/gene_data_body_lt2_qn.csv", header=TRUE, row.names = 1)

# To direct to the correct folder
date <- "2018-07-26/"
ex_dir <- "ex_8/"

#Check for presence of gene of choice

choicegene <- "GBP6"

obs_genes <- c()

for (i in 1:ncol(df.gene.body)){
  obs_genes <- c(obs_genes, strsplit(colnames(df.gene.body[i]), "_")[[1]][3])
}

ind.choicegene <- match(choicegene, obs_genes)

# Sets to run through automated graphing

datasets <- list(
  list(df.gene.body, "initial", "gene")
  ,
  list(df.gene.body.lt2, "log2 transformed", "gene_lt2")
  ,
  list(df.gene.body.lt2.qn, "log2 transformed and quantile normalised", "lt2_qn")
)

for (d in 1:length(datasets)){
  ds.data <- datasets[[d]][[1]]
  ds.ttl_adj <- datasets[[d]][[2]]
  ds.abrv <- datasets[[d]][[3]]
  
  # Ensure only that first GBP6 probe's data is used
  ds.data.choicegene <- ds.data[,ind.choicegene]
  
  #Select HIV- patients
  
  ind.hiv_neg <- c()
  
  for (i in 1:length(ds.data.choicegene)){
    if((df.gene.meta$group[i] == 1) || (df.gene.meta$group[i] == 3) || (df.gene.meta$group[i] == 6)){
      ind.hiv_neg <- c(ind.hiv_neg, i)
    }
  }
  
  ds.data.choicegene.hiv_neg <- ds.data.choicegene[ind.hiv_neg]
  df.gene.meta.hiv_neg <- df.gene.meta[ind.hiv_neg,]
  
  # Create melt table across seperate data and metadata
  Label <- c()
  variable <- c()
  value <- c()
  
  for (i in 1:length(ds.data.choicegene.hiv_neg)){
    Label <- c(Label, as.character(df.gene.meta.hiv_neg$tb.status[i]))
    variable <- c(variable, ds.ttl_adj)
    value <- c(value, ds.data.choicegene[i])
  } 
  
  gene.hiv_neg.melt <- data.frame(Label, variable, value)

  #Boxplot
  
  plot.melt <- ggplot(data = gene.hiv_neg.melt, 
                      aes(x=variable, y=value)) + 
    geom_boxplot(aes(fill=Label)) + 
    theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
    labs(title=paste("Distribution of", choicegene, "expression values", sep=" "),
         x ="Transformation", 
         y = "Expression"
    )
  
  ggsave(paste("../../img/", ex_dir, date, choicegene, "_",  ds.abrv, "_trans_bplot.png", sep=""),
         plot = plot.melt,
         device = "png",
         width = 10,        
         height = 5,
         dpi = 300,
         limitsize = FALSE
  )
}
