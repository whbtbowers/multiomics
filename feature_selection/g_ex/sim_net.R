setwd("/project/home17/whb17/Documents/project3/project_files/feature_selection/g_ex/")

library(heatmap3)
library(SNFtool)

df.gene.data<- read.csv("../../data/gene_ex2/gene_test_train.csv", header=TRUE, row.names = 1)
df.gene.data.body <- df.gene.data[,-c(1:7)]

#Select HIV- patients

ind.hiv_neg <- c()

for (i in 1:nrow(df.gene.data)){
  if((df.gene.data$group[i] == 1) || (df.gene.data$group[i] == 3) || (df.gene.data$group[i] == 6)){
    ind.hiv_neg <- c(ind.hiv_neg, i)
  }
}

df.gene.data.body.hiv_neg <- df.gene.data.body[ind.hiv_neg,]
df.gene.data.hiv_neg <- df.gene.data[ind.hiv_neg,]

diff.gene.hiv_neg <- as.matrix(dist(df.gene.data.body.hiv_neg))

AffMat.gene.hiv_neg <- affinityMatrix(diff.gene.hiv_neg, K=10)

#Create and save heatmap

png("img/2018-07-10/gene_affmat_mat_hmap1.png",
    width = 3000,        # 5 x 300 pixels
    height = 3000,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

#heatmap3(AffMat.gene.hiv_neg,
#        margins=c(10,5),
#        main= "HIV- patients, affinity matrix",
#        showColDendro=FALSE,
#        showRowDendro=FALSE
#        )

heatmap(AffMat.gene.hiv_neg,
        margins=c(10,5),
        main= "HIV- patients, affinity matrix"
        )

dev.off()


