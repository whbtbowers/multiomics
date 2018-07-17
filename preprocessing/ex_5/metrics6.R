setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_6/")

library(ggplot2)
library(stats)
library(ggfortify)
library(heatmap3)


#df.prot.data <- read.csv("../../data/ex_5/prot_data_body_logtransformed_filtimp_k20.csv", header=TRUE, row.names = 1) # Corrected and log transformed protein data
df.prot.data <- read.csv("../../data/ex_5/prot_data_body_logtransformed_filtimp_k20_wnormed.csv", header=TRUE, row.names = 1) # Wang normalised protein data
df.prot.meta <- read.csv("../../data/ex_4/prot_data_meta.csv", header=TRUE, row.names = 1) # Protein meta data

df.prot.meta$group <- as.character(df.prot.meta$group)

#df.gene.data <- read.csv("../../data/ex_5/gene_data_body_logtransformed.csv", header=TRUE, row.names = 1) # Corrected and log transformed protein data
df.gene.data <- read.csv("../../data/ex_5/gene_data_body_logtransformed_wnormed.csv", header=TRUE, row.names = 1) # Wang normalised gene data
df.gene.meta <- read.csv("../../data/ex_4/gene_data_meta.csv", header=TRUE, row.names = 1) # Gene meta data

df.gene.meta$group <- as.character(df.gene.meta$group)

# Set date and labels
label_verbose <- "Initial"
label_abbrv <- ""
date <- "2018-07-13"

datasets <- list(
  list(df.gene.data, df.gene.meta, "gene", "gene")
  ,
  list(df.prot.data, df.prot.meta, "protein", "prot")
)

for (i in 1:length(datasets)){
  data <- datasets[[i]][[1]]
  meta <- datasets[[i]][[2]]
  verbose <- datasets[[i]][[3]]
  abbrv <- datasets[[i]][[4]]
  
  # Select HIV- TB & OD and TB & LTBI patients
  ind.hiv_neg.tb_od <- c()
  ind.hiv_neg.tb_ltbi <- c()
  
  for (j in 1:nrow(data)){
    if ((meta$group[j] == 1) || (meta$group[j] == 6)){
      ind.hiv_neg.tb_od <- c(ind.hiv_neg.tb_od, j)
    }
  }
  
  for (j in 1:nrow(data)){
    if ((meta$group[j] == 1) || (meta$group[j] == 3)){
      ind.hiv_neg.tb_ltbi <- c(ind.hiv_neg.tb_ltbi, j)
    }
  }
  
  data.hiv_neg.tb_od <- data[ind.hiv_neg.tb_od,]
  meta.hiv_neg.tb_od <- meta[ind.hiv_neg.tb_od,]
 
  data.hiv_neg.tb_ltbi <- data[ind.hiv_neg.tb_ltbi,]
  meta.hiv_neg.tb_ltbi <- meta[ind.hiv_neg.tb_ltbi,]
  
  #Organise comparisons
  comps <- list(
    list(data.hiv_neg.tb_od, meta.hiv_neg.tb_od, "TB vs OD", "_tb_od")
    ,
    list(data.hiv_neg.tb_ltbi, meta.hiv_neg.tb_ltbi, "TB vs LTBI", "_tb_ltbi")
  )
  
  for (k in 1:length(comps)){
    comp_data <- comps[[k]][[1]]
    comp_meta <- comps[[k]][[2]]
    comp_verbose <- comps[[k]][[3]]
    comp_abbrv <- comps[[k]][[4]]
    
    #Boxplot
    png(paste("../../img/ex_6/", date, "/", label_abbrv,  abbrv, comp_abbrv, "_bplot.png", sep=""),
        width = 10000,        
        height = 2500,
        res = 300,          
        pointsize = 12)        
    
    plot.new()
    par(xaxt="n", mar=c(10,5,3,1))
    boxplot(comp_data, col="gold")
    lablist<-as.vector(colnames(comp_data))
    axis(1, at=seq(1, ncol(comp_data), by=1), labels = FALSE)
    text(seq(1, ncol(comp_data), by=1), par("usr")[3] - 0.2, labels = lablist, srt = 90, pos = 2, xpd = TRUE)
    title(paste(label_verbose, verbose, comp_verbose, "data", sep=" "))
    dev.off()
    
    # Heatmap
    png(paste("../../img/ex_6/", date, "/", label_abbrv, abbrv, comp_abbrv, "_hmap.png", sep=""),
        width = 5*300,        # 5 x 300 pixels
        height = 5*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8)        # smaller font size
    
    heatmap3(comp_data,
             #Rowv=NA,
             #Colv=NA,
             margins=c(10,5),
             main= paste(label_verbose, verbose, comp_verbose, "data", sep=" "),
             xlab= "prots",
             ylab= "Patients",
    )
    dev.off()
    
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
}

