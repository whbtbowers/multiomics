setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/")

df.prot_data <- read.csv("../data/protein_data_hivplus_filtimp_k10.csv", header=TRUE, row.names = 1)

df.prot_data.body <- df.prot_data[,-c(1,2)]

library(glmnet)

df.prot_glm <- glmnet(data.matrix(df.prot_data.body), 
                      data.matrix(df.prot_data$inf.status), 
                      family="gaussian"
                      )

png("../preprocessing/prot_glmnet.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8         # smaller font size
    )

plot(df.prot_glm, label=TRUE)

dev.off()

coef(df.prot_glm,s=0.1)
