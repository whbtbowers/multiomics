setwd("/project/home17/whb17/Documents/project3/project_files/feature_selection/g_ex/")

library(limma)

#df.gene_data <- read.csv("../../data/gene_ex2/labd_gene_data_ex.csv", header=TRUE, row.names = 1) #Load in inscaled data
df.gene_data <- read.csv("../../data/gene_ex2/labd_gene_data_ex_scaled.csv", header=TRUE, row.names = 1) #Load in scaled data

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

# Analyse for differential expression
tb_status.factor <- factor(df.gene_data.hiv_neg$tb.status)
site.factor <- factor(df.gene_data.hiv_neg$site)

mat.design <- model.matrix(~site.factor+tb_status.factor)
colnames(mat.design) <- c("Intercept", "OD", "TB", "ML")

fit <- lmFit(t(df.gene_data.body.hiv_neg), mat.design)
fit <- eBayes(fit, trend=TRUE, robust=TRUE)
results <- decideTests(fit)
summary(results)
tab.res <- topTable(fit, coef="ML", n=ncol(df.gene_data.body.hiv_neg))

# Toggle if graph composed
comp_img <- "Y"
#comp_img <- "N"

if (comp_img == "Y"){
  png("img/2018-07-05/gene_diff_ex_scaled.png",
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8        # smaller font size
      )
    
  plotMD(fit,coef="ML",status=results[,4],values=c(1,-1),hl.col=c("red","blue"))
  dev.off()
}


bh.P_value <- p.adjust(tab.res$P.Value, method="BH", n=length(tab.res$P.Value))

# Get significantly differentially expressed genes

ind.dif_eq <- c()

for (i in 1:length(bh.P_value)){
  if (bh.P_value[i] < 0.5){
    ind.dif_eq <- c(ind.dif_eq, i)
  }
}

library(glmnet)

fit.prot_glm <- glmnet(as.matrix(df.gene_data.body.hiv_neg), df.gene_data.hiv_neg$group, family="gaussian", alpha=0.5)

if (comp_img == "Y"){
  png("img/2018-07-09/gene_coefficients.png",
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8        # smaller font size
  )
  plot(fit.prot_glm, label=TRUE)
  dev.off()
}

coeff.prot_glm <- coef(fit.prot_glm, s=0.4)
inds.prot_glm <- which(coeff.prot_glm!=0)

choice_genes <- row.names(coeff.prot_glm)[inds.prot_glm]
`%ni%` <- Negate(`%in%`)
choice_genes <-choice_genes[choice_genes %ni% '(Intercept)']

# Get indices from big protein data for chosen genes
ind.sel.prot_data <- c()

for (i in 1:length(choice_genes)){
  ind.sel.prot_data <- c(ind.sel.prot_data, match(choice_genes[i], colnames(df.gene_data.body.hiv_neg)))
}

# Replot with fewer genes
fit.sel.prot_glm <- glmnet(as.matrix(df.gene_data.body.hiv_neg[,ind.sel.prot_data]), df.gene_data.hiv_neg$group, family="gaussian", alpha=0.5)

if (comp_img == "Y"){
  png("img/2018-07-09/selected_gene_coefficients.png",
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8        # smaller font size
  )
  plot(fit.sel.prot_glm, label=TRUE)
  dev.off()
}

df.de_sel.gene_data <- cbind(df.gene_data[,c(1:7)], df.gene_data.body[,ind.sel.prot_data])

write.csv(df.de_sel.gene_data, "../../data/gene_ex2/labd_24sel_gene_data_ex_scaled.csv", row.names=TRUE)
