setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/p_ex/")

library(glmnet)

df.prot_data.hiv_pos <- read.csv("../../data/protein_extraction_retry/protein_data_hiv_pos_filtimp_k10.csv", header=TRUE, row.names = 1)


##Convert infection status to category number
df.prot_data.hiv_pos$inf.status <- as.character(df.prot_data.hiv_pos$inf.status)
df.prot_data.hiv_pos$inf.status[df.prot_data.hiv_pos$inf.status == "TB+/HIV+"] <- 1
df.prot_data.hiv_pos$inf.status[df.prot_data.hiv_pos$inf.status == "LTBI/HIV+"] <- 2
df.prot_data.hiv_pos$inf.status[df.prot_data.hiv_pos$inf.status == "HIV+/Inf Not TB"] <- 3

# Extract just data
df.prot_data.hiv_pos.body <- df.prot_data.hiv_pos[,-c(1,2)]

#df.prot_data.hiv_pos.body <- scale(df.prot_data.hiv_pos.body)

prot_data.hiv_pos.body.glm <- glmnet(as.matrix(df.prot_data.hiv_pos.body),
                                     as.matrix(df.prot_data.hiv_pos$inf.status),
                                     family="mgaussian"
                                     )

path.glm <- paste("img/2018-07-02/prot_glmnet_hiv_pos.png", sep="")

png(path.glm,
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8         # smaller font size
)

plot(prot_data.hiv_pos.body.glm)

dev.off()

coef(prot_data.hiv_pos.body.glm, s=0.5)

