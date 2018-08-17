#setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/ex_9/")
setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_9/")

#library("lumi")
library(preprocessCore)
library(ggfortify)
library(ggplot2)
library(BiocInstaller)
library(SNFtool)

set.seed(12)

#df.gene.body <- read.csv("../../data/ex_8/gene_data_body.csv", header=TRUE, row.names = 1)
#df.gene.body <- read.csv("../../data/ex_9/prot_data_body_l2t_qn.csv", header=TRUE, row.names = 1)
df.gene.body <- read.csv("../../data/ex_8/prot_data_body.csv", header=TRUE, row.names = 1)
df.meta <- read.csv("../../data/ex_8/gp_data_meta.csv", header=TRUE, row.names = 1)
df.meta$group <- as.character(df.meta$group)

# Get row and column names
probe_names <- colnames(df.gene.body)
patient_labels <- rownames(df.gene.body)

# To direct to the correct folder
date <- "2018-08-17/"
ex_dir <- "ex_10/"

# So boxplots samples
#t.df.gene.body <- t(df.gene.body)

#Keep same orientation but be lazy
t.df.gene.body <- t(df.gene.body)


# Correct for values less than or equal to 0
if (min(t.df.gene.body) <= 0){
  corr <- abs(min(t.df.gene.body)) + 0.005
  t.df.gene.body <- t.df.gene.body + corr
}


# log2 transform data only
t.df.gene.body.log2 <- as.data.frame(log2(t.df.gene.body))

# Retranspose for saving
#df.gene.body.log2 <- as.data.frame(t(t.df.gene.body.log2))
#colnames(df.gene.body.log2) <- probe_names
#rownames(df.gene.body.log2) <- patient_labels

colnames(t.df.gene.body.log2) <- probe_names
rownames(t.df.gene.body.log2) <- patient_labels

write.csv(t.df.gene.body.log2, "../../data/ex_9/prot_data_body_l2t_BO.csv")

# quantile normalise data only
t.df.gene.body.qn <- as.data.frame(normalize.quantiles(as.matrix(t.df.gene.body)))


# quantile normalise data after log2 transformation
t.df.gene.body.log2.qn <- as.data.frame(normalize.quantiles(as.matrix(t.df.gene.body.log2)))

#df.gene.body.log2.qn <- as.data.frame(t(t.df.gene.body.log2.qn))
#colnames(df.gene.body.log2.qn) <- probe_names
#rownames(df.gene.body.log2.qn) <- patient_labels

colnames(t.df.gene.body.log2.qn) <- probe_names
rownames(t.df.gene.body.log2.qn) <- patient_labels

write.csv(t.df.gene.body.log2.qn, "../../data/ex_9/prot_data_body_l2t_qn_BO.csv")
