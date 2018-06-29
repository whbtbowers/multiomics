setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/g_ex/")

df.prot_data.hiv_pos <- read.csv("../../data/gene_ex/gene_data_hiv_pos.csv", header=TRUE, row.names = 1)
df.prot_data.hiv_neg <- read.csv("../../data/gene_ex/gene_data_hiv_neg.csv", header=TRUE, row.names = 1)

list.datasets <- list(list(df.prot_data.hiv_pos, "HIV positive patients", "hiv_pos"),
                      list(df.prot_data.hiv_neg, "HIV negative patients", "hiv_neg")
                      )

# Create graphs to display metrics

for (set in list.datasets){
  data <- set[1]
  verbose <- set[2]
  abbrv <- set[3]
  
  # Give only main body data
  data.body <- data[[1]][,-c(1,2)]
  
  #####################################################
  #                                                   #
  #                  Unadjusted data                  #
  #                                                   #
  #####################################################
  
  #Boxplot unscaled data
  path.bplot <- paste("../../data/gene_ex/img/gene_bplot_", abbrv[1], ".png", sep="")
  png(path.bplot,
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  
  plot.new()
  par(xaxt="n", mar=c(10,5,3,1))
  boxplot(data.body, col="gold")
  lablist<-as.vector(colnames(data.body))
  axis(1, at=seq(1, ncol(data.body), by=1), labels = FALSE)
  text(seq(1, ncol(data.body), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
  title(verbose)
  dev.off()
  
  # Heatmap for HIV+ and HIV-
  path.hmap <- paste("../../data/gene_ex/img/gene_hmap_", abbrv[1], ".png", sep="")
  png(path.hmap,
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  
  heatmap3(data.body,
           Rowv=NA,
           Colv=NA,
           margins=c(10,5),
           main= verbose,
           xlab= "genes",
           ylab= "Patients",
  )
  dev.off()
  
  # Principal components 1 & 2
  path.pca1_2 <- paste("../../data/gene_ex/img/gene_pca_pc1_pc2_", abbrv[1], ".png", sep="")
  png(path.pca1_2,
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  
  autoplot(prcomp(data.body), data = data[[1]], colour = "inf.status", main=verbose[1])
  dev.off()
  
  # Principal components 3 & 4
  
  path.pca3_4 <- paste("../../data/gene_ex/img/gene_pca_pc3_pc4_", abbrv[1], ".png", sep="")
  png(path.pca3_4,
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  
  autoplot(prcomp(data.body), x=3, y=4, data = data[[1]], colour = "inf.status", main=verbose[1])
  dev.off()
  
  #####################################################
  #                                                   #
  #                    Scaled data                    #
  #                                                   #
  #####################################################
  
  data.body.scaled <- scale(data.body)
  
  #Boxplot scaled data
  path.bplot <- paste("../../data/gene_ex/img/gene_bplot_", abbrv[1], "_scaled.png", sep="")
  png(path.bplot,
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  
  plot.new()
  par(xaxt="n", mar=c(10,5,3,1))
  boxplot(data.body.scaled, col="gold")
  lablist<-as.vector(colnames(data.body))
  axis(1, at=seq(1, ncol(data.body), by=1), labels = FALSE)
  text(seq(1, ncol(data.body), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
  title(verbose)
  dev.off()
  
  # Heatmap for scaled data
  path.hmap <- paste("../../data/gene_ex/img/gene_hmap_", abbrv[1], "_scaled.png", sep="")
  png(path.hmap,
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  
  heatmap3(data.body.scaled,
           Rowv=NA,
           Colv=NA,
           margins=c(10,5),
           main= verbose,
           xlab= "genes",
           ylab= "Patients",
  )
  dev.off()
  
  # Principal components 1 & 2 for scaled data
  
  path.pca1_2 <- paste("../../data/gene_ex/img/gene_pca_pc1_pc2_", abbrv[1], "_scaled.png", sep="")
  png(path.pca1_2,
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  
  autoplot(prcomp(data.body.scaled), data = data[[1]], colour = "inf.status", main=verbose[1])
  dev.off()
  
  # Principal components 3 & 4 for scaled data
  
  path.pca3_4 <- paste("../../data/gene_ex/img/gene_pca_pc3_pc4_", abbrv[1], "_scaled.png", sep="")
  png(path.pca3_4,
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8)        # smaller font size
  
  autoplot(prcomp(data.body.scaled), x=3, y=4, data = data[[1]], colour = "inf.status", main=verbose[1])
  dev.off()
  
}