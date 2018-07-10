#setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/g_ex/")
setwd("/project/home17/whb17/Documents/project3/project_files/feature_selection/g_ex/")

#df.prot_data <- read.csv("../../data/gene_ex2/gene_data_ex.csv", header=TRUE, row.names = 1) #Initial, uncorrected, unlabelled
df.prot_data <- read.csv("../../data/gene_ex2/labd_24sel_gene_data_ex_scaled.csv", header=TRUE, row.names = 1) # Scaled, labelled, 24 selected by elastic net
df.prot_data$group <- as.character(df.prot_data$group)

#Select HIV- patients

ind.hiv_neg <- c()

for (i in 1:nrow(df.gene_data)){
  if((df.gene_data$group[i] == 1) || (df.gene_data$group[i] == 3) || (df.gene_data$group[i] == 6)){
    ind.hiv_neg <- c(ind.hiv_neg, i)
  }
}

df.prot_data <- df.gene_data[ind.hiv_neg,]

library(ggplot2)
library(stats)
library(ggfortify)
library(heatmap3)

  
# Give only main body data
df.prot_data.body <- df.prot_data[,-c(1:7)]

#####################################################
#                                                   #
#                  Unadjusted data                  #
#                                                   #
#####################################################

#Boxplot unscaled data
path.bplot <- paste("img/2018-07-06/24sel_gene_bplot.png", sep="")
png(path.bplot,
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

plot.new()
par(xaxt="n", mar=c(10,5,3,1))
boxplot(df.prot_data.body, col="gold")
lablist<-as.vector(colnames(df.prot_data.body))
axis(1, at=seq(1, ncol(df.prot_data.body), by=1), labels = FALSE)
text(seq(1, ncol(df.prot_data.body), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
title("Gene data")
dev.off()

# Heatmap for HIV+ and HIV-
path.hmap <- paste("img/2018-07-06/24sel_gene_hmap.png", sep="")
png(path.hmap,
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap3(df.prot_data.body,
         Rowv=NA,
         Colv=NA,
         margins=c(10,5),
         main= "Gene data",
         xlab= "genes",
         ylab= "Patients",
)
dev.off()

# Principal components 1 & 2 by location
path.pca1_2 <- paste("img/2018-07-06/24sel_gene_pca_pc1_pc2_bysite.png", sep="")
png(path.pca1_2,
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.body), data = df.prot_data, colour = "site", main="Gene data")
dev.off()

# Principal components 3 & 4 by location

path.pca3_4 <- paste("img/2018-07-06/24sel_gene_pca_pc3_pc4_bysite.png", sep="")
png(path.pca3_4,
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.body), x=3, y=4, data = df.prot_data, colour = "site", main="Gene data")
dev.off()

# Principal components 1 & 2 by infection status
path.pca1_2 <- paste("img/2018-07-06/24sel_gene_pca_pc1_pc2_byinf.png", sep="")
png(path.pca1_2,
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.body), data = df.prot_data, colour = "group", main="Gene data")
dev.off()

# Principal components 3 & 4 by infection status

path.pca3_4 <- paste("img/2018-07-06/24sel_gene_pca_pc3_pc4_byinf.png", sep="")
png(path.pca3_4,
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.body), x=3, y=4, data = df.prot_data, colour = "group", main="Gene data")
dev.off()

#####################################################
#                                                   #
#                    Scaled data                    #
#                                                   #
#####################################################

df.prot_data.body.scaled <- scale(df.prot_data.body)

#Boxplot scaled data
path.bplot <- paste("img/2018-07-06/24sel_gene_bplot_scaled.png", sep="")
png(path.bplot,
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

plot.new()
par(xaxt="n", mar=c(10,5,3,1))
boxplot(df.prot_data.body.scaled, col="gold")
lablist<-as.vector(colnames(df.prot_data.body))
axis(1, at=seq(1, ncol(df.prot_data.body), by=1), labels = FALSE)
text(seq(1, ncol(df.prot_data.body), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
title("Gene data")
dev.off()

# Heatmap for scaled data
path.hmap <- paste("img/2018-07-06/24sel_gene_hmap_scaled.png", sep="")
png(path.hmap,
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap3(df.prot_data.body.scaled,
         Rowv=NA,
         Colv=NA,
         margins=c(10,5),
         main= "Gene data",
         xlab= "genes",
         ylab= "Patients",
)
dev.off()

# Principal components 1 & 2 for scaled data

path.pca1_2 <- paste("img/2018-07-06/24sel_gene_pca_pc1_pc2_scaled_bysite.png", sep="")
png(path.pca1_2,
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.body.scaled), data = df.prot_data, colour = "site", main="Gene data")
dev.off()

# Principal components 3 & 4 for scaled data

path.pca3_4 <- paste("img/2018-07-06/24sel_gene_pca_pc3_pc4_scaled_bysite.png", sep="")
png(path.pca3_4,
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.body.scaled), x=3, y=4, data = df.prot_data, colour = "site", main="Gene data")
dev.off()

