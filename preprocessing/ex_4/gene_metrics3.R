setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/ex_4/")

df.gene.data <- read.csv("../../data/ex_4/gene_data_body_wnormed.csv", header=TRUE, row.names = 1) # Wang normalised
#df.gene.data <- read.csv("../../data/ex_4/gene_data_body.csv", header=TRUE, row.names = 1) # Original
df.gene.meta <- read.csv("../../data/ex_4/gene_data_meta.csv", header=TRUE, row.names = 1) # Meta data

df.gene.meta$group <- as.character(df.gene.meta$group)

#Select HIV- patients

ind.hiv_neg <- c()

for (i in 1:nrow(df.gene.data)){
  if((df.gene.meta$group[i] == 1) || (df.gene.meta$group[i] == 3) || (df.gene.meta$group[i] == 6)){
    ind.hiv_neg <- c(ind.hiv_neg, i)
  }
}

df.gene.data <- df.gene.data[ind.hiv_neg,]
df.gene.meta <- df.gene.meta[ind.hiv_neg,]
  
library(ggplot2)
library(stats)
library(ggfortify)
library(heatmap3)

#Boxplot unscaled data
png("img/2018-07-10/gene_wnorm_bplot.png",
    width = 10000,        
    height = 2500,
    res = 300,          
    pointsize = 12)        

plot.new()
par(xaxt="n", mar=c(10,5,3,1))
boxplot(df.gene.data, col="gold")
lablist<-as.vector(colnames(df.gene.data))
axis(1, at=seq(1, ncol(df.gene.data), by=1), labels = FALSE)
text(seq(1, ncol(df.gene.data), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
title("Initial gene data")
dev.off()

# Heatmap
png("img/2018-07-10/gene_wnorm_hmap.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap3(df.gene.data,
         Rowv=NA,
         Colv=NA,
         margins=c(10,5),
         main= "Gene data",
         xlab= "genes",
         ylab= "Patients",
)
dev.off()

# Principal components 1 & 2 by location
path.pca1_2 <- paste("img/2018-07-10/gene_wnorm_pca_pc1_pc2_bysite.png", sep="")
png(path.pca1_2,
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.gene.data), data = df.gene.meta, colour = "site", main="Gene data")
dev.off()

# Principal components 3 & 4 by location

path.pca3_4 <- paste("img/2018-07-10/gene_wnorm_pca_pc3_pc4_bysite.png", sep="")
png(path.pca3_4,
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.gene.data), x=3, y=4, data = df.gene.meta, colour = "site", main="Gene data")
dev.off()

# Principal components 1 & 2 by infection status
path.pca1_2 <- paste("img/2018-07-10/gene_wnorm_pca_pc1_pc2_byinf.png", sep="")
png(path.pca1_2,
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.gene.data), data = df.gene.meta, colour = "group", main="Gene data")
dev.off()

# Principal components 3 & 4 by infection status

path.pca3_4 <- paste("img/2018-07-10/gene_wnorm_pca_pc3_pc4_byinf.png", sep="")
png(path.pca3_4,
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.gene.data), x=3, y=4, data = df.gene.meta, colour = "group", main="Gene data")
dev.off()
