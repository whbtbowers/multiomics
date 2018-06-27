df.prot_data.all <- read.csv("../data/protein_extraction_retry/protein_data_hiv_all_filtimp_k10.csv", header=TRUE, row.names = 1)

df.prot_data.all.body <- df.prot_data.all[,-c(1,2)]

par(xaxt="n")
boxplot(df.prot_data.all.body, col="gold",)
lablist<-as.vector(colnames(df.prot_data.all.body))
axis(1, at=seq(1, ncol(df.prot_data.all.body), by=1), labels = FALSE)
text(seq(1, ncol(df.prot_data.all.body), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
#dev.off()
