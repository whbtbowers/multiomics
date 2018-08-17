setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_9/")

library("lumi")

df.prot.body <- read.csv("../../data/ex_8/prot_data_body.csv", header=TRUE, row.names = 1)


# To direct to the correct folder
date <- "2018-07-30/"
ex_dir <- "ex_9/"

# So boxplots samples
t.df.prot.body <- t(df.prot.body)

#Boxplot initial data
png(paste("../../img/", ex_dir, date, "prot_patient_bplot.png", sep=""),
    width = 10000,        
    height = 2500,
    res = 300,          
    pointsize = 12)        

plot.new()
par(xaxt="n", mar=c(10,5,3,1))
boxplot(t.df.prot.body, col="green")
lablist<-as.vector(colnames(t.df.prot.body))
axis(1, at=seq(1, ncol(t.df.prot.body), by=1), labels = FALSE)
text(seq(1, ncol(t.df.prot.body), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
title("Initial protein data")
dev.off()

# PCA plot
# Principal components 1 & 2 by infection status

png(paste("../../img/ex_6/", date, "/", label_abbrv, abbrv, comp_abbrv, "_pca_pc1_pc2_byinf.png", sep=""),
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size


autoplot(prcomp(comp_data),
         data = comp_meta,
         colour = "group",
         main=paste(label_verbose, verbose, comp_verbose, "data", sep=" ")
)

dev.off()

# Principal components 3 & 4 by infection status

png(paste("../../img/ex_6/", date, "/", label_abbrv, abbrv, comp_abbrv, "_pca_pc3_pc4_byinf.png", sep=""),
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(comp_data),
         x=3, 
         y=4,
         data = comp_meta,
         colour = "group",
         main=paste(label_verbose, verbose, comp_verbose, "data", sep=" ")
)

dev.off()

}

