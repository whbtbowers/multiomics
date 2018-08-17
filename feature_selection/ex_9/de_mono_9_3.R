setwd("/home/whb17/Documents/project3/project_files/feature_selection/ex_9/")

library(limma)
library(heatmap3)
library(SNFtool)
library(glmnet)

set.seed(12)

#df.gene.body <- read.csv("../../data/ex_9/gene_train_body.csv", header=TRUE, row.names = 1)  # Protein test/train set

df.prot.body <- read.csv("../../data/ex_9/prot_train_body.csv", header=TRUE, row.names = 1)  # Protein training set

# Barplot to quickly check distributions
#par(xaxt="n", mar=c(10,5,3,1))
#boxplot(df.prot.body, col="red")
#lablist<-as.vector(colnames(df.prot.body))
#axis(1, at=seq(1, ncol(df.prot.body), by=1), labels = FALSE)
#text(seq(1, ncol(df.prot.body), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
#df.prot.body <- df.gene.body

df.meta <- read.csv("../../data/ex_9/gp_train_meta.csv", header=TRUE, row.names = 1)
#df.meta <- read.csv("../../data/ex_9/gp_data_meta.csv", header=TRUE, row.names = 1)
df.meta$group <- as.character(df.meta$group)

# To direct to the correct folder
date <- "2018-08-08/"
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
tab.res <- topTable(fit, coef=4, n=ncol(df.prot.body))

png(paste("../../img/", ex_dir, date, "prot_tb_od_meandiff.png", sep=""),
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8        # smaller font size
)

plotMD(fit, coef=4, status=results[,4], values=c(1,-1), hl.col=c("red","blue"))
dev.off()

# Check BH corrected values

ind.dif_ex <- c()


for (k in 1:length(tab.res$adj.P.Val)){
  if (tab.res$adj.P.Val[k] < 0.05){
    ind.dif_ex <- c(ind.dif_ex, k)
  }
}

sig_P = tab.res$adj.P.Val[ind.dif_ex]
sig_factor = rownames(tab.res)[ind.dif_ex]

sel.df.prot.hiv_neg.tb_od <- df.prot.hiv_neg.tb_od[match(sig_factor, colnames(df.prot.hiv_neg.tb_od))]

# Get columns from significant factors

###########################
## Elastic net selection ##
###########################

# Cross-validated analysis of coefficients to find optmal lamda cut-off

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

#Convert to dataframe for easier sorting
coef.cvfit.lambda1se <- coef(cvfit, s=cvfit$lambda.1se)

features <- coef.cvfit.lambda1se@Dimnames[[1]][which(coef.cvfit.lambda1se != 0 ) ]

df.coef.cvfit.lambda1se <- data.frame(
  features, 
  coefs <- coef.cvfit.lambda1se[which(coef.cvfit.lambda1se != 0 ) ] 
)
colnames(df.coef.cvfit.lambda1se) <- c("features", "coefficients")

#Remove intercept
if(df.coef.cvfit.lambda1se$features[1] == "(Intercept)"){
  df.coef.cvfit.lambda1se <- df.coef.cvfit.lambda1se[-1,]
  features <- features[-1]
}

# Straightforward lasso to extract beta coefficients


# Fit to elastic net

fit.glm <- glmnet(as.matrix(df.prot.hiv_neg.tb_od), 
                  df.meta.hiv_neg.tb_od$group,
                  family="gaussian", 
                  alpha=alpha
)

png(paste("../../img/", ex_dir, date, "prot_tb_od_glmnet_coeff.png", sep=""),
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8        # smaller font size
)

plot(fit.glm, label=TRUE)
dev.off()


