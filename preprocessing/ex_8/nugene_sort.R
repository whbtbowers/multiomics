setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_8/")

library(stringr)

#tot_gene_data <- read.delim("~/Documents/project3/project_files/data/tot_gene_data.txt", header=TRUE, row.names=1) # Initial text file
#write.csv(tot_gene_data, file=paste("../../data/ex_8/tot_gene_data.csv", sep=""), row.names=TRUE) # Write new csv
df.gene_init <- read.csv("../../data/ex_8/tot_gene_data.csv", header=TRUE, row.names = 1)

df.gene_init$CHROMOSOME <- NULL
df.gene_init$SYNONYMS <- NULL
df.gene_init$SEARCH_KEY <- NULL
df.gene_init$ILMN_GENE <- NULL
df.gene_init$SEARCH_KEY <- NULL

# subset to test
#df.gene_init <- df.gene_init[1:100, 1:100]

# Add gene symbol to ID and remove column
symbols <- df.gene_init$SYMBOL

for (i in 1:nrow(df.gene_init)){
  rownames(df.gene_init)[i] <- paste(rownames(df.gene_init)[i], "_", symbols[i], sep = "")
}

df.gene_init <- df.gene_init[,-1]

#Remove Pvalue columns

for (i in rev(1:ncol(df.gene_init))){
  if(str_detect(colnames(df.gene_init)[i],"Pval")){
    df.gene_init <- df.gene_init[,-i]
  }
}

# Leave colname as just array name

ar_lab <- c()

for (i in 1:ncol(df.gene_init)){
  ar_lab <- c(ar_lab, strsplit(colnames(df.gene_init)[i], ".AVG_")[[1]][1])
}

colnames(df.gene_init) <- ar_lab

t.df.gene_init <- t(df.gene_init)

write.csv(t.df.gene_init, "../../data/ex_8/tot_gene_data_phase2.csv")

#colnames(t.df.gene_init) <- t.df.gene_init[1,]



