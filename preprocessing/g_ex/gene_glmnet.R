setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/g_ex/")

df.gene_data.hiv_pos <- read.csv("../../data/gene_ex/gene_data_hiv_pos.csv", header=TRUE, row.names = 1)
df.gene_data.hiv_neg <- read.csv("../../data/gene_ex/gene_data_hiv_neg.csv", header=TRUE, row.names = 1)



## Convert to category number
df.gene_data.hiv_pos$inf.status <- as.character(df.gene_data.hiv_pos$inf.status)
df.gene_data.hiv_pos$inf.status[df.gene_data.hiv_pos$inf.status == "TB+/HIV+"] <- 1
df.gene_data.hiv_pos$inf.status[df.gene_data.hiv_pos$inf.status == "LTBI/HIV+"] <- 2
df.gene_data.hiv_pos$inf.status[df.gene_data.hiv_pos$inf.status == "HIV+/Inf Not TB"] <- 3

df.gene_data.hiv_pos.body <- data.frame(df.gene_data.hiv_pos[,-c(1,2)])

library(glmnet)

df.prot_glm <- glmnet(x=as.matrix(as.numeric(df.gene_data.hiv_pos.body)), as.matrix(df.gene_data.hiv_pos$inf.status), family="mgaussian")

png("../preprocessing/prot_glmnet.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

plot(df.prot_glm, label=TRUE)

dev.off()

coef(df.prot_glm,s=0.1)