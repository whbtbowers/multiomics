setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_7/")

library(illuminaio)
library(stringr)

inp.probe_labels <- readBGX(file.path("../../data/HumanHT-12_V4_0_R2_15002873_B.bgx"))

#inp.gene_data <- read.csv("../../data/ex_7/gene_data_body_lt_qn.csv",row.names=1, header=TRUE)
#inp.gene_labels <- read.csv("../../data/ex_7/sig.txt", row.names=1, header=TRUE)
inp.sig_genes <- read.csv("../../data/ex_7/gene_tb_od_sig_factors_alpha0.5.csv", row.names=1, header=TRUE)

df.probe_labels.probes <- inp.probe_labels$probes
df.probe_labels.controls <-inp.probe_labels$controls

write.csv(df.probe_labels.probes, "../../data/ex_7/probes.csv")

# For sig factors:
labs.sig_gene <- rownames(inp.sig_genes)


#start_time <- Sys.time()

# Find certain genes

#ind.found <- c()

#to_find <- c("ACTB", "HUPO")

#for(gene in to_find){
#  ind.found <- c(ind.found, match(gene, df.probe_labels.probes$ILMN_Gene))
#}



ref_id <- df.probe_labels.probes$Array_Address_Id

#labels.gene_data <- colnames(inp.gene_data.body)

#Get gene names
ind.gene_labels <- c()
gene_labels <- c()
combo_labs <- c()
synonyms <- c()

# If list of sig factors
for (i in 1:length(labs.sig_gene)){
  ind <- match(as.numeric(substring(labs.sig_gene[i],2)), ref_id)
  ind.gene_labels <- c(ind.gene_labels, ind)
  gene_labels <- c(gene_labels, df.probe_labels.probes$ILMN_Gene[ind])
  combo_labs <- c(combo_labs, paste(labs.sig_gene[i], "_", gene_labels[i], sep=""))
  synonyms <- c(synonyms, df.probe_labels.probes$Synonyms[ind])
}

for (i in 2:length(synonyms)){
  if (str_detect(toString(synonyms[i]), "GBP6")){
    print(paste(synonyms[i], ", ", gene_labels[i]))
  }
}

rownames(inp.sig_genes) <- combo_labs

inp.sig_genes <- cbind(inp.sig_genes, synonyms)

write.csv(inp.sig_genes, "../../data/ex_7/gene_tb_od_sig_factors_alpha0.5.csv")

# If just list of labels
#for (i in 1:nrow(inp.gene_labels)){
#  ind <- match(as.numeric(substring(inp.gene_labels[i,],2)), ref_id)
#  ind.gene_labels <- c(ind.gene_labels, ind)
#  gene_labels <- c(gene_labels, df.probe_labels.probes$ILMN_Gene[ind])
#}

#write.csv(gene_labels, "../../data/ex_7/gene_labels.csv")

# If entire csv
#for (i in 1:length(labels.gene_data)){
#  ind <- match(as.numeric(substring(labels.gene_data[i],2)), ref_id)
#  ind.gene_labels <- c(ind.gene_labels, ind)
  #gene_labels <- c(gene_labels, df.probe_labels.probes$ILMN_Gene[ind])
#}

#gene_labels <- df.probe_labels.probes$ILMN_Gene[ind.gene_labels]

#colnames(inp.gene_data)[8:ncol(inp.gene_data)] <- gene_labels

#write.csv(inp.gene_data, "../data/gene_ex2/labd_gene_data_ex.csv", row.names=TRUE)

#inp.gene_data.body.scaled <- scale(inp.gene_data[,-c(1:7)])

#inp.gene_data.scaled <- cbind(inp.gene_data[,c(1:7)],inp.gene_data.body.scaled)

#write.csv(inp.gene_data.scaled, "../data/gene_ex2/labd_gene_data_ex_scaled.csv", row.names=TRUE)
