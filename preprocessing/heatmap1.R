library("ggplot2")

df.prot_data <- read.csv("../data/protein_data_filtimp_k10.csv")

#Remove patient numbers
df.noX_prot_data <- df.prot_data[,-1]

# Patient numbers and protien names
pat_num <- paste(rep("PATIENT", 20), df.prot_data$X, sep="_")
prot_name <- colnames(df.noX_prot_data)

df.prot_data$X <- pat_num

# Melt dataframe
prot_heatmap <- melt(df.prot_data, id.vars = "X")
names(prot_heatmap)[2:3] <- c("proteins", "expression_level")

head(prot_heatmap)

ggplot(prot_heatmap, aes(proteins, X)) +
  geom_tile(aes(fill = expression_level), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ylab("Patients") +
  xlab("Proteins") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Expression level")