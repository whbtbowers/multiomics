setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_7")

library(VennDiagram)

df.meta <- read.csv("../../data/ex_7/gene_prot_data_meta.csv", header=TRUE, row.names = 1, na.strings = c("", "NA")) # Complete meta

#for(i in 1:6){
#  print(paste("## GROUP ", i, " ##", sep=""))
#  print(length(df.meta$group[df.meta$group == as.character(i)]))
#  print(round(length(df.meta$group[df.meta$group == as.character(i)]) * 0.2))
#  print(round(length(df.meta$group[df.meta$group == as.character(i)]) * 0.8))
#}

# Set date
date <- "2018-07-16"

# Get indices of patients based on TB status
ind.tb <- which(df.meta$tb.status %in% "TB")
ind.ltbi <- which(df.meta$tb.status %in% "LTBI")
ind.od <- which(df.meta$tb.status %in% "TB")

#Get patients based on TB status
df.tb <- df.meta[ind.tb,]
df.ltbi <- df.meta[ind.ltbi,]
df.od <- df.meta[ind.od,]

inf_cats <- list(
  list(df.tb, "tb"),
  list(df.ltbi, "ltbi"),
  list(df.od, "od")
)


for (i in 1:length(inf_cats)){
  data <- inf_cats[[i]][[1]][1]
  abrv <- inf_cats[[i]][[2]][1]

  # Create venn diagram
  
  venn.diagram(x = list(
    "Gene data" = df.gene.meta[,colnum],
    "Protein data" = df.prot.meta[,colnum]
  ),
  filename = paste("../../img/ex_7/", date, "/venn_", abrv,".png", sep=""),
  imagetype= "png"
  )
}


