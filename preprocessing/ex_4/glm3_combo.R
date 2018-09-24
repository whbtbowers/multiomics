setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_5/")

library(glmnet)

# Gene data
#df.gene.data <- read.csv("../../data/ex_4/gene_data_body_wnormed.csv", header=TRUE, row.names = 1) # Wang normalised
#df.gene.data <- read.csv("../../data/ex_4/gene_data_body.csv", header=TRUE, row.names = 1) # Original
df.gene.data <- read.csv("../../data/ex_5/gene_data_body_logtransformed_wnormed.csv", header=TRUE, row.names = 1) # Wang normalised gene data
df.gene.meta <- read.csv("../../data/ex_4/gene_data_meta.csv", header=TRUE, row.names = 1) # Meta
#df.gene.meta$group <- as.character(df.gene.meta$group)

# Protein data
#df.prot.data <- read.csv("../../data/ex_4/prot_data_body_wnormed.csv", header=TRUE, row.names = 1) # Wang normalised
#df.prot.data <- read.csv("../../data/ex_4/prot_data_body.csv", header=TRUE, row.names = 1) # Original
df.prot.data <- read.csv("../../data/ex_5/prot_data_body_logtransformed_filtimp_k20_wnormed.csv", header=TRUE, row.names = 1) # Wang normalised protein data
df.prot.meta <- read.csv("../../data/ex_4/prot_data_meta.csv", header=TRUE, row.names = 1) # Meta
#df.prot.meta$group <- as.character(df.prot.meta$group)

datasets <- list(
  list(df.gene.data, df.gene.meta, "gene", "gene")
  ,
  list(df.prot.data, df.prot.meta, "protein", "prot")
)

#Set various metrics
date <- "2018-07-13"
alphas <- c(0, 0.5, 1)

for (alpha in alphas){
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
      comp.data <- comps[[k]][[1]]
      comp.meta <- comps[[k]][[2]]
      comp.verbose <- comps[[k]][[3]]
      comp.abbrv <- comps[[k]][[4]]
      
      #Plot individual coefficient graphs
      
      glm.data <- glmnet(data.matrix(comp.data), 
                          data.matrix(comp.meta$group), 
                          family="gaussian",
                          alpha=alpha
      )
      
      png(paste("../../img/ex_5/",  date, "/", abbrv, comp.abbrv, "_glmnet_alpha", alpha, ".png", sep=""),
          width = 5*300,        # 5 x 300 pixels
          height = 5*300,
          res = 300,            # 300 pixels per inch
          pointsize = 8         # smaller font size
      )
      plot(glm.data, label=TRUE)
      
      dev.off()
    }
    
    #Plot superimposed cross validated graph 
    
    png(paste("../../img/ex_5/",  date, "/", abbrv,"_superimposed_cv_glmnet_alpha", alpha, ".png", sep=""),
        width = 5*300,        # 5 x 300 pixels
        height = 5*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8         # smaller font size
    )
    
    for (k in 1:length(comps)){
      comp.data <- comps[[k]][[1]]
      comp.meta <- comps[[k]][[2]]
      comp.verbose <- comps[[k]][[3]]
      comp.abbrv <- comps[[k]][[4]]
      
      cvfit <- cv.glmnet(data.matrix(comp.data), 
                         data.matrix(comp.meta$group), 
                         family="gaussian",
                         alpha=alpha
      )
      plot(cvfit)
      par(new = TRUE)
    }
    
    dev.off()
    
    
  }
}