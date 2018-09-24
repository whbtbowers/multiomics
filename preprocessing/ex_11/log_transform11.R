#setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/ex_11/")
setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_11/")

#library("lumi")
library(preprocessCore)
library(ggfortify)
library(ggplot2)
library(BiocInstaller)
library(SNFtool)

set.seed(12)

df.gene.body <- read.csv("../../data/ex_11/gene_data_body.csv", header=TRUE, row.names = 1)
df.prot.body <- read.csv("../../data/ex_11/prot_data_body.csv", header=TRUE, row.names = 1)
df.meta <- read.csv("../../data/ex_11/gp_data_meta.csv", header=TRUE, row.names = 1)
df.meta$group <- as.character(df.meta$group)

# Get row and column names
patient_labels <- rownames(df.gene.body)
probe_names <- colnames(df.gene.body)
prot_names <- colnames(df.prot.body)


# To direct to the correct folder
date <- "2018-08-20/"
ex_dir <- "ex_11/"

# So boxplots samples
#t.df.gene.body <- t(df.gene.body)

#Keep same orientation but be lazy
t.df.gene.body <- t(df.gene.body)
t.df.prot.body <- t(df.prot.body)

# Correct for values less than or equal to 0
if (min(t.df.gene.body) <= 0){
  corr <- abs(min(t.df.gene.body)) + 0.005
  t.df.gene.body <- t.df.gene.body + corr
}

# Correct for values less than or equal to 0
if (min(t.df.prot.body) <= 0){
  corr <- abs(min(t.df.prot.body)) + 0.005
  t.df.prot.body <- t.df.prot.body + corr
}


# log2 transform both datasets
t.df.gene.body.log2 <- as.data.frame(log2(t.df.gene.body))
t.df.prot.body.log2 <- as.data.frame(log2(t.df.prot.body))

# Retranspose for saving
df.gene.body.log2 <- as.data.frame(t(t.df.gene.body.log2))
colnames(df.gene.body.log2) <- probe_names
rownames(df.gene.body.log2) <- patient_labels

#colnames(t.df.gene.body.log2) <- probe_names
#rownames(t.df.gene.body.log2) <- patient_labels

df.prot.body.log2 <- as.data.frame(t(t.df.prot.body.log2))
colnames(df.prot.body.log2) <- prot_names
rownames(df.prot.body.log2) <- patient_labels

write.csv(df.gene.body.log2, "../../data/ex_11/gene_data_body_l2t.csv")
write.csv(df.prot.body.log2, "../../data/ex_11/prot_data_body_l2t.csv")

# Quantile normalise gene data after log2 transformation
t.df.gene.body.log2.qn <- as.data.frame(normalize.quantiles(as.matrix(t.df.gene.body.log2)))

# Standard normalise protein data after log2 transformation
t.df.prot.body.log2.sn <- as.data.frame(standardNormalization(as.matrix(t.df.prot.body.log2)))

# Retransform for saving
df.gene.body.log2.qn <- as.data.frame(t(t.df.gene.body.log2.qn))
colnames(df.gene.body.log2.qn) <- probe_names
rownames(df.gene.body.log2.qn) <- patient_labels


df.prot.body.log2.sn <- as.data.frame(t(t.df.prot.body.log2.sn))
colnames(df.prot.body.log2.sn) <- prot_names
rownames(df.prot.body.log2.sn) <- patient_labels

write.csv(df.gene.body.log2.qn, "../../data/ex_11/gene_data_body_l2t_qn.csv")
write.csv(df.prot.body.log2.sn, "../../data/ex_11/prot_data_body_l2t_sn.csv")
