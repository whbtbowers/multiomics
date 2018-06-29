library(ggplot2)

df.prot_data.all <- read.csv("../data/protein_extraction_retry/protein_data_hiv_all_filtimp_k10.csv", header=TRUE, row.names = 1)

df.prot_data.all.body <- df.prot_data.all[,-c(1,2)]

# HIV+ and HIV-

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

df.body.scaled <- scale(df.prot_data.all.body)

# HIV+ and HIV-, scaled

png("../data/protein_extraction_retry/img/prot_bplot_scaled.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

plot.new()
par(xaxt="n", mar=c(10,5,1,1))
boxplot(df.body.scaled, col="gold")
lablist<-as.vector(colnames(df.prot_data.all.body))
axis(1, at=seq(1, ncol(df.prot_data.all.body[-8]), by=1), labels = FALSE)
text(seq(1, ncol(df.prot_data.all.body[-8]), by=1), par("usr")[3] - 0.2, labels = lablist[-8], srt = 90, pos = 2, xpd = TRUE)
dev.off()

