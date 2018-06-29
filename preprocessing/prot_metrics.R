library(ggplot2)
library(stats)
library(ggfortify)

df.prot_data.all <- read.csv("../data/protein_extraction_retry/protein_data_hiv_all_filtimp_k10.csv", header=TRUE, row.names = 1)
df.prot_data.hiv_pos <- read.csv("../data/protein_extraction_retry/protein_data_hiv_pos_filtimp_k10.csv", header=TRUE, row.names = 1)
df.prot_data.hiv_neg <- read.csv("../data/protein_extraction_retry/protein_data_hiv_neg_filtimp_k10.csv", header=TRUE, row.names = 1)

df.prot_data.all.body <- df.prot_data.all[,-c(1,2)]
df.prot_data.hiv_pos.body <- df.prot_data.hiv_pos[,-c(1,2)]
df.prot_data.hiv_neg.body <- df.prot_data.hiv_neg[,-c(1,2)]

#####################################################
#                                                   #
#                Both HIV categories                #
#                                                   #
#####################################################

# Boxplot for HIV+ and HIV-

png("../data/protein_extraction_retry/img/prot_bplot_all.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

plot.new()
par(xaxt="n", mar=c(10,5,1,1))
boxplot(df.prot_data.all.body, col="gold")
lablist<-as.vector(colnames(df.prot_data.all.body))
axis(1, at=seq(1, ncol(df.prot_data.all.body), by=1), labels = FALSE)
text(seq(1, ncol(df.prot_data.all.body), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
dev.off()

# Heatmap for HIV+ and HIV-

png("../data/protein_extraction_retry/img/prot_hmap_all.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap3(df.prot_data.all.body,
         Rowv=NA,
         Colv=NA,
         margins=c(10,5),
         xlab="Proteins",
         ylab="Patients",
)
dev.off()

# Principal components 1 & 2 for HIV+ and HIV-

png("../data/protein_extraction_retry/img/prot_pca_pc1_pc2_hiv_all.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.all.body), data = df.prot_data.all, colour = "inf.status")
dev.off()

# Principal components 3 & 4 for HIV+ and HIV-

png("../data/protein_extraction_retry/img/prot_pca_pc3_pc4_hiv_all.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.all.body), x=3, y=4, data = df.prot_data.all, colour = "inf.status")
dev.off()

#####################################################
#                                                   #
#                        HIV+                       #
#                                                   #
#####################################################

# Boxplot for HIV+

png("../data/protein_extraction_retry/img/prot_bplot_hiv_pos.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

plot.new()
par(xaxt="n", mar=c(10,5,1,1))
boxplot(df.prot_data.hiv_pos.body, col="gold")
lablist<-as.vector(colnames(df.prot_data.hiv_pos.body))
axis(1, at=seq(1, ncol(df.prot_data.hiv_pos.body), by=1), labels = FALSE)
text(seq(1, ncol(df.prot_data.hiv_pos.body), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
dev.off()

# Heatmap for HIV+

png("../data/protein_extraction_retry/img/prot_hmap_hiv_pos.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap3(df.prot_data.hiv_pos.body,
         Rowv=NA,
         Colv=NA,
         margins=c(10,5),
         xlab="Proteins",
         ylab="Patients",
)
dev.off()

# Principal components 1 & 2 for HIV+

png("../data/protein_extraction_retry/img/prot_pca_pc1_pc2_hiv_pos.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.hiv_pos.body), data = df.prot_data.hiv_pos, colour = "inf.status")
dev.off()

# Principal components 3 & 4 for HIV+

png("../data/protein_extraction_retry/img/prot_pca_pc3_pc4_hiv_pos.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.hiv_pos.body), x=3, y=4, data = df.prot_data.hiv_pos, colour = "inf.status")
dev.off()

#####################################################
#                                                   #
#                        HIV-                       #
#                                                   #
#####################################################

# Boxplot for HIV-

png("../data/protein_extraction_retry/img/prot_bplot_hiv_neg.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

plot.new()
par(xaxt="n", mar=c(10,5,1,1))
boxplot(df.prot_data.hiv_neg.body, col="gold")
lablist<-as.vector(colnames(df.prot_data.hiv_neg.body))
axis(1, at=seq(1, ncol(df.prot_data.hiv_neg.body), by=1), labels = FALSE)
text(seq(1, ncol(df.prot_data.hiv_neg.body), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
dev.off()

# Heatmap for HIV-

png("../data/protein_extraction_retry/img/prot_hmap_hiv_neg.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap3(df.prot_data.hiv_neg.body,
         Rowv=NA,
         Colv=NA,
         margins=c(10,5),
         xlab="Proteins",
         ylab="Patients",
)
dev.off()

# Principal components 1 & 2 for HIV-

png("../data/protein_extraction_retry/img/prot_pca_pc1_pc2_hiv_neg.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.hiv_neg.body), data = df.prot_data.hiv_neg, colour = "inf.status")
dev.off()

# Principal components 3 & 4 for HIV-

png("../data/protein_extraction_retry/img/prot_pca_pc3_pc4_hiv_neg.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.hiv_neg.body), x=3, y=4, data = df.prot_data.hiv_neg, colour = "inf.status")
dev.off()

## Scale all data ##
df.all.body.scaled <- scale(df.prot_data.all.body)
df.hiv_pos.body.scaled <- scale(df.prot_data.hiv_pos.body)
df.hiv_neg.body.scaled <- scale(df.prot_data.hiv_neg.body)

#####################################################
#                                                   #
#            Both HIV categories scaled             #
#                                                   #
#####################################################

# Boxplot for scaled HIV+ and HIV-

png("../data/protein_extraction_retry/img/prot_bplot_all_scaled.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

plot.new()
par(xaxt="n", mar=c(10,5,1,1))
boxplot(df.all.body.scaled, col="gold")
lablist<-as.vector(colnames(df.prot_data.all.body))
axis(1, at=seq(1, ncol(df.prot_data.all.body[-8]), by=1), labels = FALSE)
text(seq(1, ncol(df.prot_data.all.body[-8]), by=1), par("usr")[3] - 0.2, labels = lablist[-8], srt = 90, pos = 2, xpd = TRUE)
dev.off()

# Heatmap for scaled HIV+ and HIV-

png("../data/protein_extraction_retry/img/prot_hmap_all_scaled.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size


heatmap3(df.all.body.scaled,
         Rowv=NA,
         Colv=NA,
         margins=c(10,5),
         xlab="Proteins",
         ylab="Patients",
)
dev.off()

# First 2 principal components for HIV+ and HIV-

png("../data/protein_extraction_retry/img/prot_pca_pc1_pc2_hiv_all_scaled.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.all.body.scaled), data = df.prot_data.all, colour = "inf.status")
dev.off()

# Principal components 3 & 4 for HIV+ and HIV-

png("../data/protein_extraction_retry/img/prot_pca_pc3_pc4_hiv_all_scaled.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.all.body.scaled), x=3, y=4, data = df.prot_data.all, colour = "inf.status")
dev.off()

#####################################################
#                                                   #
#                    HIV+ scaled                    #
#                                                   #
#####################################################

# Boxplot for HIV+ scaled

png("../data/protein_extraction_retry/img/prot_bplot_hiv_pos_scaled.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

plot.new()
par(xaxt="n", mar=c(10,5,1,1))
boxplot(df.hiv_pos.body.scaled, col="gold")
lablist<-as.vector(colnames(df.hiv_pos.body.scaled))
axis(1, at=seq(1, ncol(df.prot_data.hiv_pos.body), by=1), labels = FALSE)
text(seq(1, ncol(df.prot_data.hiv_pos.body), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
dev.off()

# Heatmap for HIV+ scaled

png("../data/protein_extraction_retry/img/prot_hmap_hiv_pos_scaled.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap3(df.hiv_pos.body.scaled,
         Rowv=NA,
         Colv=NA,
         margins=c(10,5),
         xlab="Proteins",
         ylab="Patients",
)
dev.off()
 
# Principal components 1 & 2 for HIV+ scaled

png("../data/protein_extraction_retry/img/prot_pca_pc1_pc2_hiv_pos_scaled.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.hiv_pos.body.scaled), data = df.prot_data.hiv_pos, colour = "inf.status")
dev.off()

# Principal components 3 & 4 for HIV+ scaled

png("../data/protein_extraction_retry/img/prot_pca_pc3_pc4_hiv_pos_scaled.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.hiv_pos.body.scaled), x=3, y=4, data = df.prot_data.hiv_pos, colour = "inf.status")
dev.off()

#####################################################
#                                                   #
#                    HIV- scaled                    #
#                                                   #
#####################################################

# Boxplot for HIV- scaled

png("../data/protein_extraction_retry/img/prot_bplot_hiv_neg_scaled.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

plot.new()
par(xaxt="n", mar=c(10,5,1,1))
boxplot(df.hiv_neg.body.scaled, col="gold")
lablist<-as.vector(colnames(df.prot_data.hiv_pos.body))
axis(1, at=seq(1, ncol(df.prot_data.hiv_pos.body), by=1), labels = FALSE)
text(seq(1, ncol(df.prot_data.hiv_pos.body), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
dev.off()

# Heatmap for HIV- scaled

png("../data/protein_extraction_retry/img/prot_hmap_hiv_neg_scaled.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap3(df.hiv_neg.body.scaled,
         Rowv=NA,
         Colv=NA,
         margins=c(10,5),
         xlab="Proteins",
         ylab="Patients",
)
dev.off()

# Principal components 1 & 2 for HIV- scaled

png("../data/protein_extraction_retry/img/prot_pca_pc1_pc2_hiv_neg_scaled.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.hiv_neg.body.scaled), data = df.prot_data.hiv_neg, colour = "inf.status")
dev.off()

# Principal components 3 & 4 for HIV- scaled

png("../data/protein_extraction_retry/img/prot_pca_pc3_pc4_hiv_neg_scaled.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.hiv_neg.body.scaled), x=3, y=4, data = df.prot_data.hiv_neg, colour = "inf.status")
dev.off()
