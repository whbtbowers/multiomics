library(reshape2)
library(ggplot2)

# Row and column labels
name_genes <- paste(rep("GEN", 20), LETTERS[1:20], sep="_") # rows
name_patients <- paste(rep("PATIENT", 20), seq(1,20,1), sep="_") # columns

# Generation of dataframe
value_expression <- data.frame(genes = name_genes, 
                               matrix(rnorm(400, 2, 1.8),nrow = 20, ncol = 20))
names(value_expression)[2:21] <- name_patients

# Melt dataframe
df_heatmap <- melt(value_expression, id.vars = "genes")
names(df_heatmap)[2:3] <- c("patient", "expression_level")

head(df_heatmap)

# Elaboration of heatmap (white - steelblue)

ggplot(df_heatmap, aes(patient, genes )) +
  geom_tile(aes(fill = expression_level), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  ylab("List of genes ") +
  xlab("List of patients") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(fill = "Expression level")
