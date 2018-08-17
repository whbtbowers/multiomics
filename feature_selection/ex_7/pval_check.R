setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_7/")

library(limma)
library(glmnet)

set.seed(12)

# Protein data
#df.prot.tt <- read.csv("../../data/ex_7/prot_data_body_lt2.csv", header=TRUE, row.names = 1)  # Protein l2t data

df.prot.tt <- read.csv("../../data/ex_7/gene_testtrain_body.csv", header=TRUE, row.names = 1)  # Protein l2t/qn test/train set


# Complete meta
df.all.meta.tt <- read.csv("../../data/ex_7/gene_prot_testtrain_meta.csv", header=TRUE, row.names = 1)
df.all.meta.tt$group <- as.character(df.all.meta.tt$group)

# To direct to the correct folder
date <- "2018-07-20/"
ex_dir <- "ex_7/"

# Parameters
alpha = 0.5
K = 20

#Select HIV- TB vs LTBI patients
ind.hiv_neg.tb_ltbi <- c()

for (i in 1:nrow(df.all.meta.tt)){
  if((df.all.meta.tt$group[i] == 1) || (df.all.meta.tt$group[i] == 3)){
    ind.hiv_neg.tb_ltbi <- c(ind.hiv_neg.tb_ltbi, i)
  }
}

set.data.hiv_neg.tb_ltbi <- df.prot.tt[ind.hiv_neg.tb_ltbi,]
set.meta.hiv_neg.tb_ltbi <- df.all.meta.tt[ind.hiv_neg.tb_ltbi,]

tb_status.factor <- factor(set.meta.hiv_neg.tb_ltbi$tb.status)
site.factor <- factor(set.meta.hiv_neg.tb_ltbi$site)
sex.factor <- factor(set.meta.hiv_neg.tb_ltbi$sex)

mat.design <- model.matrix(~1 + tb_status.factor + sex.factor + site.factor)
#mat.design <- model.matrix(~0 + tb_status.factor)

fit <- lmFit(t(set.data.hiv_neg.tb_ltbi), mat.design)
fit <- eBayes(fit, trend=TRUE, robust=TRUE)
results <- decideTests(fit)
print(summary(results))
tab.res <- topTable(fit, coef = "tb_status.factorTB", n=ncol(df.prot.tt))

# Select by BH corrected P-values
ind.dif_ex <- c()


for (i in 1:length(tab.res$adj.P.Val)){
  if (tab.res$adj.P.Val[i] > 0.05){
    ind.dif_ex <- c(ind.dif_ex, i)
  }
}

sig_P = tab.res$adj.P.Val[ind.dif_ex]
sig_factor = rownames(tab.res)[ind.dif_ex]
df.sig_factor = data.frame(sig_factor)
write.csv(sig_factor, "../../data/ex_7/sig.txt")
