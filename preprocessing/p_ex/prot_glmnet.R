setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/p_ex/")

df.prot_data.hiv_pos <- read.csv("../../data/protein_extraction_retry/protein_data_hiv_pos_filtimp_k10.csv", header=TRUE, row.names = 1)
df.prot_data.hiv_neg <- read.csv("../../data/protein_extraction_retry/protein_data_hiv_neg_filtimp_k10.csv", header=TRUE, row.names = 1)

df.prot_data.hiv_pos.body <- df.prot_data.hiv_pos[,-c(1,2)]
df.prot_data.hiv_neg.body <- df.prot_data.hiv_neg[,-c(1,2)]

library(glmnet)

# Get matrices of GLM data
df.prot_glm.hiv_pos <- glmnet(data.matrix(df.prot_data.hiv_pos.body), data.matrix(df.prot_data.hiv_pos$inf.status), family="mgaussian")
df.prot_glm.hiv_neg <- glmnet(data.matrix(df.prot_data.hiv_neg.body), data.matrix(df.prot_data.hiv_neg$inf.status), family="mgaussian")


png("img/2018-07-02/prot_glmnet_hiv_pos.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

plot(df.prot_glm.hiv_pos, label=TRUE)
dev.off()

coef(df.prot_glm.hiv_pos,s=0.1)

png("img/2018-07-02/prot_glmnet_hiv_pos.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

plot(df.prot_glm.hiv_pos, label=TRUE)
dev.off()

coef(df.prot_glm.hiv_pos,s=0.1)