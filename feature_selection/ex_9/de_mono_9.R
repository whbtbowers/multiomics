setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_9/")

library(limma)
library(heatmap3)
library(SNFtool)
library(glmnet)

set.seed(12)

#df.gene.body <- read.csv("../../data/ex_9/gene_train_body.csv", header=TRUE, row.names = 1)  # Protein test/train set

#df.prot.body <- read.csv("../../data/ex_9/prot_train_body.csv", header=TRUE, row.names = 1)  # Protein train set

df.prot.body <- read.csv("../../data/ex_9/prot_data_body_l2t_qn_BO.csv", header=TRUE, row.names = 1)

#df.prot.body <- df.gene.body

#df.meta <- read.csv("../../data/ex_9/gp_train_meta.csv", header=TRUE, row.names = 1)
df.meta <- read.csv("../../data/ex_9/gp_data_meta.csv", header=TRUE, row.names = 1)
df.meta$group <- as.character(df.meta$group)

# To direct to the correct folder
date <- "2018-07-31/"
ex_dir <- "ex_9/"

# Parameters
alpha = 0.5
K = 20

#Select HIV- TB and OD patients

ind.hiv_neg.tb_od <- c()

for (i in 1:nrow(df.meta)){
  if((df.meta$group[i] == 1) || (df.meta$group[i] == 6)){
    ind.hiv_neg.tb_od <- c(ind.hiv_neg.tb_od, i)
  }
}

df.prot.hiv_neg.tb_od <- df.prot.body[ind.hiv_neg.tb_od,]
df.meta.hiv_neg.tb_od <- df.meta[ind.hiv_neg.tb_od,]

#############################
## Limma-based DE analysis ##
#############################

# Make factors for analysis to consider

fac.sex <- factor(df.meta.hiv_neg.tb_od$sex)
fac.site <- factor(df.meta.hiv_neg.tb_od$site) # Correct for site? Maybe use as intercept
fac.tb <- factor(df.meta.hiv_neg.tb_od$tb.status)

#design <- model.matrix(~fac.site + fac.sex + fac.tb)
design <- model.matrix(~0 + fac.site + fac.sex + fac.tb)

fit <- lmFit(t(df.prot.hiv_neg.tb_od), design)
fit <- eBayes(fit, trend=TRUE, robust=TRUE)
results <- decideTests(fit)
summary(results)
tab.res <- topTable(fit, coef=4, n=22)

png(paste("../../img/", ex_dir, date, "prot_tb_od_meandiff_BO.png", sep=""),
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8        # smaller font size
)

plotMD(fit, coef=4, status=results[,4], values=c(1,-1), hl.col=c("red","blue"))
dev.off()

# Check BH corrected values

ind.dif_ex <- c()


for (i in 1:length(tab.res$adj.P.Val)){
  if (tab.res$adj.P.Val[i] < 0.05){
    ind.dif_ex <- c(ind.dif_ex, i)
  }
}

sig_P = tab.res$adj.P.Val[ind.dif_ex]
sig_factor = rownames(tab.res)[ind.dif_ex]

# Get columns from significant factors

###########################
## Elastic net selection ##
###########################

# Fit to elastic net

fit.glm <- glmnet(as.matrix(df.prot.hiv_neg.tb_od), 
                  df.meta.hiv_neg.tb_od$group, 
                  family="gaussian", 
                  alpha=alpha
)

png(paste("../../img/", ex_dir, date, "prot_tb_od_glmnet_coeff_BO.png", sep=""),
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8        # smaller font size
)

plot(fit.glm, label=TRUE)
dev.off()

# Cross-validated analysis of coefficients

cvfit <- cv.glmnet(data.matrix(df.prot.hiv_neg.tb_od), 
                   data.matrix(as.numeric(df.meta.hiv_neg.tb_od$group)), 
                   family="gaussian",
                   alpha=alpha
)

png(paste("../../img/", ex_dir, date, "prot_tb_od_cv_glmnet_coeff.png", sep=""),
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8        # smaller font size
)

plot(cvfit, main="Cross-validated coefficient plot for HIV- TB vs OD protein data")
dev.off()

# Out of curiosity, comparing rige, EMN, and lasso
foldid=sample(1:K,size=length(data.matrix(as.numeric(df.meta.hiv_neg.tb_od$group))),replace=TRUE)
cv1=cv.glmnet(data.matrix(df.prot.hiv_neg.tb_od),data.matrix(as.numeric(df.meta.hiv_neg.tb_od$group)),foldid=foldid,alpha=1)
cv.5=cv.glmnet(data.matrix(df.prot.hiv_neg.tb_od),data.matrix(as.numeric(df.meta.hiv_neg.tb_od$group)),foldid=foldid,alpha=.5)
cv0=cv.glmnet(data.matrix(df.prot.hiv_neg.tb_od),data.matrix(as.numeric(df.meta.hiv_neg.tb_od$group)),foldid=foldid,alpha=0)

png(paste("../../img/", ex_dir, date, "prot_tb_od_alpha_comp.png", sep=""),
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8        # smaller font size
)

par(mfrow=c(2,2))
plot(cv1);plot(cv.5);plot(cv0)
plot(log(cv1$lambda),cv1$cvm,pch=19,col="red",xlab="log(Lambda)",ylab=cv1$name)
points(log(cv.5$lambda),cv.5$cvm,pch=19,col="grey")
points(log(cv0$lambda),cv0$cvm,pch=19,col="blue")
legend("topleft",legend=c("alpha= 1","alpha= .5","alpha 0"),pch=19,col=c("red","grey","blue"))
dev.off()

