setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/ex_5/")

library(glmnet)

# Gene data
#df.gene.data <- read.csv("../../data/ex_4/gene_data_body_wnormed.csv", header=TRUE, row.names = 1) # Wang normalised
#df.gene.data <- read.csv("../../data/ex_4/gene_data_body.csv", header=TRUE, row.names = 1) # Original
df.gene.data <- read.csv("../../data/ex_5/gene_data_body_logtransformed_wnormed.csv", header=TRUE, row.names = 1) # Wang normalised gene data
df.gene.meta <- read.csv("../../data/ex_4/gene_data_meta.csv", header=TRUE, row.names = 1) # Meta
df.gene.meta$group <- as.character(df.gene.meta$group)

# Protein data
#df.prot.data <- read.csv("../../data/ex_4/prot_data_body_wnormed.csv", header=TRUE, row.names = 1) # Wang normalised
#df.prot.data <- read.csv("../../data/ex_4/prot_data_body.csv", header=TRUE, row.names = 1) # Original
df.prot.data <- read.csv("../../data/ex_5/prot_data_body_logtransformed_filtimp_k20_wnormed.csv", header=TRUE, row.names = 1) # Wang normalised protein data
df.prot.meta <- read.csv("../../data/ex_4/prot_data_meta.csv", header=TRUE, row.names = 1) # Meta
df.prot.meta$group <- as.character(df.prot.meta$group)

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
  
  # Select HIV- patients
  ind.hiv_neg <- c()
  
  for (j in 1:nrow(data)){
    if ((meta$group[j] == 1) || (meta$group[j] == 3) || (meta$group[j] == 6)){
      ind.hiv_neg <- c(ind.hiv_neg, j)
    }
  }
  
  data.hiv_neg <- data[ind.hiv_neg,]
  meta.hiv_neg <- meta[ind.hiv_neg,]
  
  glm.data.hiv_neg <- glmnet(data.matrix(data.hiv_neg), 
                      data.matrix(meta.hiv_neg$group), 
                      family="gaussian",
                      alpha=0.5
  )
  
  png(paste("../../img/ex5/2018-07-11/", abbrv ,"_glmnet.png", sep=""),
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8         # smaller font size
  )
  plot(glm.data.hiv_neg, label=TRUE)
  
  dev.off()
  
  #Get coefficients
  coeff.glm.data.hiv_neg <- coef(glm.data.hiv_neg, s=c(0.1,0.2,0.3,0.4,0.5))
  
  # To few proteins to eliminate. usic s=0.5, returns 32 genes
  
  if (verbose == "gene"){
    nonzero.coeff.glm.data.hiv_neg <- coeff.glm.data.hiv_neg[,5][coeff.glm.data.hiv_neg[,5]!=0][-1] # Get nonzero coefficients without intercept
    ind.nonzero.coeff.glm.data.hiv_neg <- match(names(nonzero.coeff.glm.data.hiv_neg), colnames(data))
    
    # Select genes above set lamda threshold
    sel.data.hiv_neg <- data.hiv_neg[,ind.nonzero.coeff.glm.data.hiv_neg]
    sel.data <- data[,ind.nonzero.coeff.glm.data.hiv_neg]
    
    #GLM plot selected genes
    glm.sel.data.hiv.neg <- glmnet(data.matrix(sel.data.hiv_neg), 
                                   data.matrix(meta.hiv_neg$group), 
                                   family="gaussian", 
                                   alpha=0.5
                                   )
    
    png(paste("../../img/ex5/2018-07-11/sel_", abbrv ,"_glmnet.png", sep=""),
        width = 5*300,        # 5 x 300 pixels
        height = 5*300,
        res = 300,            # 300 pixels per inch
        pointsize = 8         # smaller font size
    )
    plot(glm.sel.data.hiv.neg, label=TRUE)
    
    dev.off()
    
    write.csv(sel.data, file=paste("../../data/ex_4/sel_", abbrv, "_data_body.csv", sep=""), row.names=TRUE)
    
  }
}
