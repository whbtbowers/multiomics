setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/ex_9/")
#setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_9/")

library("lumi")
library(preprocessCore)
library(ggfortify)
library(ggplot2)

set.seed(12)
df.gene.body <- read.csv("../../data/ex_8/gene_data_body.csv", header=TRUE, row.names = 1)
#df.prot.body <- read.csv("../../data/ex_8/prot_data_body.csv", header=TRUE, row.names = 1)
df.meta <- read.csv("../../data/ex_8/gp_data_meta.csv", header=TRUE, row.names = 1)
df.meta$group <- as.character(df.meta$group)

#t.df.gene.body <- read.csv("../../data/ex_8/gene_data_body_t.csv", header=TRUE, row.names = 1)




#transform all files
#write.csv(t(df.meta), "../../data/ex_8/gp_data_meta_t.csv")
#write.csv(t(df.gene.body), "../../data/ex_8/gene_data_body_t.csv")

#t.df.gene.body <- lumiR.batch("../../data/ex_8/gene_data_body_t.csv", sep=",")

# Get row and column names
probe_names <- colnames(df.gene.body)
patient_labels <- rownames(df.gene.body)

# To direct to the correct folder
date <- "2018-07-30/"
ex_dir <- "ex_9/"

# So boxplots samples
t.df.gene.body <- t(df.gene.body)

#Boxplot initial data
png(paste("../../img/", ex_dir, date, "gene_patient_bplot.png", sep=""),
    width = 10000,        
    height = 2500,
    res = 300,          
    pointsize = 12)        

plot.new()
par(xaxt="n", mar=c(10,5,3,1))
boxplot(t.df.gene.body, col="red")
lablist<-as.vector(colnames(t.df.gene.body))
axis(1, at=seq(1, ncol(t.df.gene.body), by=1), labels = FALSE)
text(seq(1, ncol(t.df.gene.body), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
title("Initial gene data")
dev.off()

# PCA plot
# Principal components 1 & 2 by infection status

#png(paste("../../img/ex_6/", date, "/", label_abbrv, abbrv, comp_abbrv, "_pca_pc1_pc2_byinf.png", sep=""),
#    width = 5*300,        # 5 x 300 pixels
#    height = 5*300,
#    res = 300,            # 300 pixels per inch
#    pointsize = 8)        # smaller font size


#autoplot(prcomp(comp_data),
#         data = comp_meta,
#         colour = "group",
#         main=paste(label_verbose, verbose, comp_verbose, "data", sep=" ")
#)

#dev.off()

# Principal components 3 & 4 by infection status

#png(paste("../../img/ex_6/", date, "/", label_abbrv, abbrv, comp_abbrv, "_pca_pc3_pc4_byinf.png", sep=""),
#    width = 5*300,        # 5 x 300 pixels
#    height = 5*300,
#    res = 300,            # 300 pixels per inch
#   pointsize = 8)        # smaller font size
#
#autoplot(prcomp(comp_data),
#         x=3, 
#         y=4,
#         data = comp_meta,
#         colour = "group",
#         main=paste(label_verbose, verbose, comp_verbose, "data", sep=" ")
#)

#dev.off()

#}


# log2 transform data
t.df.gene.body.log2 <- as.data.frame(log2(t.df.gene.body))
#t.df.gene.body.log2 <- lumiT(as.matrix(t.df.gene.body), method='log2') # Must be eSet for lumiT

#Boxplot logged data
png(paste("../../img/", ex_dir, date, "l2t_gene_patient_bplot.png", sep=""),
    width = 10000,        
    height = 2500,
    res = 300,          
    pointsize = 12)        

plot.new()
par(xaxt="n", mar=c(10,5,3,1))
boxplot(t.df.gene.body.log2, col="red")
lablist<-as.vector(colnames(t.df.gene.body.log2))
axis(1, at=seq(1, ncol(t.df.gene.body.log2), by=1), labels = FALSE)
text(seq(1, ncol(t.df.gene.body.log2), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
title("Log2 transformed gene data")
dev.off()

# Retranspose for PCA
df.gene.body.log2 <- as.data.frame(t(t.df.gene.body.log2))
colnames(df.gene.body.log2) <- probe_names
rownames(df.gene.body.log2) <- patient_labels

write.csv(df.gene.body.log2, "../../data/ex_9/gene_data_body_l2t.csv")

# PCA plot
# Principal components 1 & 2 by infection status

png(paste("../../img/", ex_dir, date, "gene_l2t_pca_pc1_pc2_byinf.png", sep=""),
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size


autoplot(prcomp(df.gene.body.log2),
         data = df.meta,
         colour = "group",
         main="Log2 transformed gene data"
)

dev.off()

# Principal components 3 & 4 by infection status

png(paste("../../img/", ex_dir, date, "/", "gene_l2t_pca_pc3_pc4_byinf.png", sep=""),
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
   pointsize = 8)        # smaller font size

autoplot(prcomp(df.gene.body.log2),
         x=3, 
         y=4,
         data = df.meta,
         colour = "group",
         main="Log2 transformed gene data"
)

dev.off()

  
# quantile normalise data
t.df.gene.body.log2.qn <- as.data.frame(normalize.quantiles(t.df.gene.body.log2))
#t.df.gene.body.log2.qn <- lumiN(as.matrix(t.df.gene.body.log2), method='quantile')

#Boxplot logged and quantile normalised data
png(paste("../../img/", ex_dir, date, "l2t_qn_gene_patient_bplot.png", sep=""),
    width = 10000,        
    height = 2500,
    res = 300,          
    pointsize = 12)        

plot.new()
par(xaxt="n", mar=c(10,5,3,1))
boxplot(t.df.gene.body.log2.qn, col="red")
lablist<-as.vector(colnames(t.df.gene.body))
axis(1, at=seq(1, ncol(t.df.gene.body.log2.qn), by=1), labels = FALSE)
text(seq(1, ncol(t.df.gene.body.log2.qn), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
title("Log2 transformed and quantile normalised gene data")
dev.off()

df.gene.body.log2.qn <- as.data.frame(t(t.df.gene.body.log2.qn))
colnames(df.gene.body.log2.qn) <- probe_names
rownames(df.gene.body.log2.qn) <- patient_labels

write.csv(df.gene.body.log2.qn, "../../data/ex_9/gene_data_body_l2t_qn.csv")
