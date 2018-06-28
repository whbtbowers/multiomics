library(stats)
library(ggfortify)

# Load CSVs
df.prot_data.all <- read.csv("../data/prot_ex2/protein_data_hiv_all_filtimp_k10.csv", header=TRUE, row.names = 1)

# Extract data
df.prot_data.all.body <- df.prot_data.all[,-c(1,2)]

#All variables, first 2 principal components

png("img/prot_pca_pc1_pc2_hiv_all.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.all.body), data = df.prot_data.all, variance_percentage = TRUE, colour = "inf.status")
dev.off()

#All variables, 3rd and 4th principal components

png("img/prot_pca_pc3_pc4_hiv_all.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.prot_data.all.body), x=3, y=4, variance_percentage = TRUE, data = df.prot_data.all, colour = "inf.status")
dev.off()

df.body.scaled <- scale(df.prot_data.all.body)

#All variables, scaled, first 2 principal components

png("img/prot_pca_pc1_pc2_hiv_all_scaled.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.body.scaled), data = df.prot_data.all, variance_percentage = TRUE, colour = "inf.status")
dev.off()

#All variables, scaled, 3rd and 4th principal components

png("img/prot_pca_pc3_pc4_hiv_all_scaled.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

autoplot(prcomp(df.body.scaled), x=3, y=4, variance_percentage = TRUE, data = df.prot_data.all, colour = "inf.status")
dev.off()
