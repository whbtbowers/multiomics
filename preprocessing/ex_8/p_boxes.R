setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_8/")

library(limma)
library(heatmap3)
library(SNFtool)
library(glmnet)

set.seed(12)

df.gene.body <- read.csv("../../data/ex_8/gene_data_body_lt2_qn.csv", header=TRUE, row.names = 1)  # Gene test/train set

# To direct to the correct folder
date <- "2018-07-30/"
ex_dir <- "ex_8/"

# So boxplots samples
t.df.gene.body <- t(df.gene.body)

#Boxplot
png(paste("../../img/", ex_dir, date, "gene_l2t_qn_patient_bplot.png", sep=""),
    width = 10000,        
    height = 2500,
    res = 300,          
    pointsize = 12)        

plot.new()
par(xaxt="n", mar=c(10,5,3,1))
boxplot(t.df.gene.body, col="gold")
lablist<-as.vector(colnames(t.df.gene.body))
axis(1, at=seq(1, ncol(t.df.gene.body), by=1), labels = FALSE)
text(seq(1, ncol(t.df.gene.body), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
dev.off()
