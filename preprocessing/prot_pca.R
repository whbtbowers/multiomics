library(stats)
library(ggfortify)
stats::prcomp

# Load CSVs
df.prot_data.all <- read.csv("../data/protein_extraction_retry/protein_data_hiv_all_filtimp_k10.csv", header=TRUE, row.names = 1)
df.prot_data.hiv_pos <- read.csv("../data/protein_extraction_retry/protein_data_hiv_pos_filtimp_k10.csv", header=TRUE, row.names = 1)
df.prot_data.hiv_neg <- read.csv("../data/protein_extraction_retry/protein_data_hiv_neg_filtimp_k10.csv", header=TRUE, row.names = 1)

df.prot_data.all$inf.status <- as.character(df.prot_data.all$inf.status)
df.prot_data.hiv_pos$inf.status <- as.character(df.prot_data.hiv_pos$inf.status)
df.prot_data.hiv_neg$inf.status <- as.character(df.prot_data.hiv_neg$inf.status)

# Extract data
df.prot_data.all.body <- df.prot_data.all[,-c(1,2)]
df.prot_data.hiv_pos.body <- df.prot_data.hiv_pos[,-c(1,2)]
df.prot_data.hiv_neg.body <- df.prot_data.hiv_neg[,-c(1,2)]


#All variables, first 2 principal components

png("../data/extraction_retry/img/prot_pca_pc1_pc2_hiv_all.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.all.body[,-8]), data = df.prot_data.all, colour = "inf.status")
dev.off()

#All variables, 3rd and 4th principal components

png("../data/extraction_retry/img/prot_pca_pc3_pc4_hiv_all.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.all.body), x=3, y=4, data = df.prot_data.all, colour = "inf.status")
dev.off()

# HIV+ variables, first 2 principal components

png("../data/extraction_retry/img/prot_pca_pc1_pc2_hiv_pos.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.hiv_pos.body), data = df.prot_data.hiv_pos, colour = "inf.status")
dev.off()

# HIV+ variables, 3rd and 4th principal components

png("../data/extraction_retry/img/prot_pca_pc3_pc4_hiv_pos.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.hiv_pos.body), x=3, y=4, data = df.prot_data.hiv_pos, colour = "inf.status")
dev.off()

# HIV- variables, first 2 principal components

png("../data/extraction_retry/img/prot_pca_pc1_pc2_hiv_neg.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.all.body), data = df.prot_data.all, colour = "inf.status")
dev.off()

# HIV- variables, 3rd and 4th principal components

png("../data/extraction_retry/img/prot_pca_pc3_pc4_hiv_neg.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.hiv_neg.body), x=3, y=4, data = df.prot_data.hiv_neg, colour = "inf.status")
dev.off()

