setwd("/home/whb17/Documents/project3/project_files/preprocessing/ex_10/")

library(ggplot2)
library(stats)
library(ggfortify)
library(heatmap3)


# Load in meta and protein files
df.meta <- read.csv("../../data/ex_9/gp_data_meta.csv", header=TRUE, row.names = 1)
df.meta$group <- as.character(df.meta$group)

df.prot <- read.csv("../../data/ex_8/prot_data_body.csv", header=TRUE, row.names = 1)
df.prot.l2t.qn <- read.csv("../../data/ex_9/prot_data_body_l2t_qn.csv", header=TRUE, row.names = 1)

# To direct to the correct folder
date <- "2018-08-16/"
ex_dir <- "ex_10/"

# Sets to run through automated graphing

#Select HIV- patients

ind.hiv_neg <- c()

for (i in 1:nrow(df.prot)){
  if((df.meta$group[i] == 1) || (df.meta$group[i] == 3) || (df.meta$group[i] == 6)){
    ind.hiv_neg <- c(ind.hiv_neg, i)
  }
}

df.meta.hiv_neg <- df.meta[ind.hiv_neg,]
df.prot.hiv_neg <- df.prot[ind.hiv_neg,]
df.prot.l2t.qn.hiv_neg <- df.prot.l2t.qn[ind.hiv_neg,]

# Create melt table across seperate data and metadata for initial data

Label <- c()
variable <- c()
value <- c()

for (i in 1:nrow(df.prot.hiv_neg)){
  for (j in 1:ncol(df.prot.hiv_neg)){
    Label <- c(Label, as.character(df.meta.hiv_neg$tb.status[i]))
    variable <- c(variable, colnames(df.prot.hiv_neg)[j])
    value <- c(value, df.prot.hiv_neg[i,j])
  }
}

prot.hiv_neg.melt <- data.frame(Label, variable, value)

prot.hiv_neg.melt$value <- prot.hiv_neg.melt$value + 5e3

prot.hiv_neg.melt.nohap <- prot.hiv_neg.melt[-which(prot.hiv_neg.melt$variable %in% "Haptoglobin..ng.ml."),]

prot.hiv_neg.melt.nohap.noa2m <- prot.hiv_neg.melt.nohap[-which(prot.hiv_neg.melt.nohap$variable %in% "A2M..ng.ml."),]

prot.hiv_neg.melt.lowex <- prot.hiv_neg.melt[which(prot.hiv_neg.melt$variable %in% c("IFNa2..pg.ml.", "IFNg..pg.ml.","IL.1RA..pg.ml.","IP.10..pg.ml.","Procalcitonin..pg.ml.","TGF.a..pg.ml.","TNFa..pg.ml.","VEGF..pg.ml.")),]

prot.hiv_neg.melt.lowerex <- prot.hiv_neg.melt[which(prot.hiv_neg.melt$variable %in% c("IFNa2..pg.ml.", "IFNg..pg.ml.", "IL.1RA..pg.ml.", "TGF.a..pg.ml.", "TNFa..pg.ml.")),]

# Create melt table across seperate data and metadata for normalised data

Label.l2t.qn <- c()
variable.l2t.qn <- c()
value.l2t.qn <- c()

for (i in 1:nrow(df.prot.l2t.qn.hiv_neg)){
  for (j in 1:ncol(df.prot.l2t.qn.hiv_neg)){
    Label.l2t.qn <- c(Label.l2t.qn, as.character(df.meta.hiv_neg$tb.status[i]))
    variable.l2t.qn <- c(variable.l2t.qn, colnames(df.prot.l2t.qn.hiv_neg)[j])
    value.l2t.qn <- c(value.l2t.qn, df.prot.l2t.qn.hiv_neg[i,j])
  }
}

prot.l2t.qn.hiv_neg.melt <- data.frame(Label.l2t.qn, variable.l2t.qn, value.l2t.qn)

#Boxplots

#Initial
plot.melt <- ggplot(data = prot.hiv_neg.melt.lowerex, 
       aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=Label)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="Distribution of protein expression in initial data",
        x ="Proteins", 
       y = "Expression (log2 scaled)"
       ) +
  coord_trans(y = "log2")

ggsave(paste("../../img/", ex_dir, date, "/", "prot_bplot_l2ax_lowerex", sep=""),
       plot = plot.melt,
       device = "png",
       width = 10,        
       height = 5,
       dpi = 300,
       limitsize = FALSE
       )

# Transformed and Normed

plot.melt.l2t.qn <- ggplot(data = prot.l2t.qn.hiv_neg.melt, 
                    aes(x=variable.l2t.qn, y=value.l2t.qn)) + 
  geom_boxplot(aes(fill=Label.l2t.qn)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title="Distribution of protein expression in Log2 transformed and quantile normalised data",
       x ="Proteins", 
       y = "Expression"
  )

ggsave(paste("../../img/", ex_dir, date, "/", "prot_l2t_qn_bplot.png", sep=""),
       plot = plot.melt.l2t.qn,
       device = "png",
       width = 10,        
       height = 5,
       dpi = 300,
       limitsize = FALSE
)