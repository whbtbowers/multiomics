setwd("/project/home17/whb17/Documents/project3/project_files/feature_selection/g_ex/")

library(limma)

df.gene_data <- read.csv("../../data/gene_ex2/labd_gene_data_ex.csv", header=TRUE, row.names = 1) #Load in data
df.gene_data$group <- as.character(df.gene_data$group) # Ensure group is treated as categorical rather than continuous
df.gene_data.body <- df.gene_data[,-c(1:7)]

#Select HIV- patients

ind.hiv_neg <- c()

for (i in 1:nrow(df.gene_data)){
  if((df.gene_data$group[i] == 1) || (df.gene_data$group[i] == 3) || (df.gene_data$group[i] == 6)){
    ind.hiv_neg <- c(ind.hiv_neg, i)
  }
}

df.gene_data.body.hiv_neg <- df.gene_data.body[ind.hiv_neg,]
df.gene_data.hiv_neg <- df.gene_data[ind.hiv_neg,]

tb_status.factor <- factor(df.gene_data.hiv_neg$tb.status)
site.factor <- factor(df.gene_data.hiv_neg$site)

mat.design <- model.matrix(~site.factor+tb_status.factor)
colnames(mat.design) <- c("Intercept", "OD", "TB", "ML")

fit <- lmFit(t(df.gene_data.body.hiv_neg), mat.design)
fit <- eBayes(fit, trend=TRUE, robust=TRUE)
results <- decideTests(fit)
summary(results)
topTable(fit, coef="ML")


png("img/2018-07-04/gene_diff_ex.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8        # smaller font size
    )
  
plotMD(fit,coef="ML",status=results[,4],values=c(1,-1),hl.col=c("red","blue"))
dev.off()