library(heatmap3)

# Load CSVs
df.prot_data.all <- read.csv("../data/protein_extraction_retry/protein_data_hiv_all_filtimp_k10.csv", header=TRUE, row.names = 1)
df.prot_data.hiv_pos <- read.csv("../data/protein_extraction_retry/protein_data_hiv_pos_filtimp_k10.csv", header=TRUE, row.names = 1)
df.prot_data.hiv_neg <- read.csv("../data/protein_extraction_retry/protein_data_hiv_neg_filtimp_k10.csv", header=TRUE, row.names = 1)

# Extract data
df.prot_data.all.body <- df.prot_data.all[,-c(1,2)]
df.prot_data.hiv_pos.body <- df.prot_data.hiv_pos[,-c(1,2)]
df.prot_data.hiv_neg.body <- df.prot_data.hiv_neg[,-c(1,2)]

#Plot all
png("../data/extraction_retry/img/prot_hmap_all.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap3(df.prot_data.all.body[,-8],
         Rowv=NA,
         Colv=NA,
         margins=c(10,5),
         xlab="Proteins",
         ylab="Patients",
         )
dev.off()

#Plot positive
png("../data/extraction_retry/img/prot_hmap_hiv_pos.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap3(df.prot_data.hiv_pos.body,
         Rowv=NA,
         Colv=NA,
         margins=c(10,5),
         xlab="Proteins",
         ylab="Patients",
)
dev.off()

#Plot negative
png("../data/extraction_retry/img/prot_hmap_hiv_neg.png",
    width = 5*300,        # 5 x 300 pixels
    height = 5*300,
    res = 300,            # 300 pixels per inch
    pointsize = 8)        # smaller font size

heatmap3(df.prot_data.hiv_neg.body,
         Rowv=NA,
         Colv=NA,
         margins=c(10,5),
         xlab="Proteins",
         ylab="Patients",
)
dev.off()