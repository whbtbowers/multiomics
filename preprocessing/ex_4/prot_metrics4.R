setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/ex_4/")

df.prot.data <- read.csv("../../data/ex_4/prot_data_body_wnormed.csv", header=TRUE, row.names = 1) # Wang normalised
#df.prot.data <- read.csv("../../data/ex_4/prot_data_body_filtimp_k20.csv", header=TRUE, row.names = 1) # Original
df.prot.meta <- read.csv("../../data/ex_4/prot_data_meta.csv", header=TRUE, row.names = 1) # Meta data

df.prot.meta$group <- as.character(df.prot.meta$group)

# Set label and date
label <- "sel_prot"
date <- "2018-07-11"



#Select HIV- patients

ind.hiv_neg <- c()

for (i in 1:nrow(df.prot.data)){
  if((df.prot.meta$group[i] == 1) || (df.prot.meta$group[i] == 3) || (df.prot.meta$group[i] == 6)){
    ind.hiv_neg <- c(ind.hiv_neg, i)
  }
}

df.prot.data <- df.prot.data[ind.hiv_neg,]
df.prot.meta <- df.prot.meta[ind.hiv_neg,]
  
library(ggplot2)
library(stats)
library(ggfortify)
library(heatmap3)

#Boxplot unscaled data
png(paste("img/", date, "/", label, "_bplot.png", sep=""),
    width = 10000,        
    height = 2500,
    res = 300,          
    pointsize = 12)        

plot.new()
par(xaxt="n", mar=c(10,5,3,1))
boxplot(df.prot.data, col="gold")
lablist<-as.vector(colnames(df.prot.data))
axis(1, at=seq(1, ncol(df.prot.data), by=1), labels = FALSE)
text(seq(1, ncol(df.prot.data), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
title("Initial Normalised protein data")
dev.off()

# Heatmap
png(paste("img/", date, "/", label, "_hmap.png", sep=""),
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap3(df.prot.data,
         Rowv=NA,
         Colv=NA,
         margins=c(10,5),
         main= "Normalised protein data",
         xlab= "prots",
         ylab= "Patients",
)
dev.off()

# Principal components 1 & 2 by location
#path.pca1_2 <- paste("img/2018-07-11/prot_wnorm_pca_pc1_pc2_bysite.png", sep="")
#png(path.pca1_2,
#    width = 5*300,        # 5 x 300 pixels
#    height = 5*300,
#    res = 300,            # 300 pixels per inch
#    pointsize = 8)        # smaller font size

#autoplot(prcomp(df.prot.data), data = df.prot.meta, colour = "site", main="Normalised protein data")
#dev.off()

# Principal components 3 & 4 by location

#path.pca3_4 <- paste("img/2018-07-11/prot_wnorm_pca_pc3_pc4_bysite.png", sep="")
#png(path.pca3_4,
#    width = 5*300,        # 5 x 300 pixels
#    height = 5*300,
#   res = 300,            # 300 pixels per inch
#    pointsize = 8)        # smaller font size

#autoplot(prcomp(df.prot.data), x=3, y=4, data = df.prot.meta, colour = "site", main="Normalised protein data")
#dev.off()

# Principal components 1 & 2 by infection status
png(paste("img/", date, "/", label, "_pca_pc1_pc2_byinf", sep=""),
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot.data), data = df.prot.meta, colour = "group", main="Normalised protein data")
dev.off()

# Principal components 3 & 4 by infection status

path.pca3_4 <- paste("img/2018-07-11/prot_wnorm_pca_pc3_pc4_byinf.png", sep="")
png(path.pca3_4,
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot.data), x=3, y=4, data = df.prot.meta, colour = "group", main="Normalised protein data")
dev.off()

