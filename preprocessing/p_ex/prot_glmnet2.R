setwd("/project/home17/whb17/Documents/project3/project_files/preprocessing/p_ex/")

library(glmnet)

df.prot_data.hiv_pos <- read.csv("../../data/protein_extraction_retry/protein_data_hiv_pos_filtimp_k10.csv", header=TRUE, row.names = 1)
df.prot_data.hiv_neg <- read.csv("../../data/protein_extraction_retry/protein_data_hiv_neg_filtimp_k10.csv", header=TRUE, row.names = 1)

list.datasets <- list(list(df.prot_data.hiv_neg, "HIV negative patients", "hiv_neg"),
                      list(df.prot_data.hiv_pos, "HIV positive patients", "hiv_pos")
)

for (set in list.datasets){
  data <- set[[1]]
  verbose <- set[2][[1]]
  abbrv <- set[3][[1]]
  
  ##Convert infection status to category number
  data$inf.status <- as.character(data$inf.status)
  data$inf.status[data$inf.status == "TB+/HIV+"] <- as.numeric(1)
  data$inf.status[data$inf.status == "LTBI/HIV+"] <- as.numeric(2)
  data$inf.status[data$inf.status == "HIV+/Inf Not TB"] <- as.numeric(3)
  
  data.body <- data[,-c(1,2)]
  
  data.body.glm <- glmnet(x=as.matrix(data.body), y=as.matrix(data$inf.status), family="mgaussian")
  
  path.png <- paste("img/2018-07-02/prot_glmnet_", abbrv, ".png", sep="")
  png(path.png,
      width = 5*300,        # 5 x 300 pixels
      height = 5*300,
      res = 300,            # 300 pixels per inch
      pointsize = 8         # smaller font size
  )
  
  plot(data.body.glm, label=TRUE)
  dev.off()
  
  print(verbose)
  
  coef(data.body.glm, s=0.1)
}
