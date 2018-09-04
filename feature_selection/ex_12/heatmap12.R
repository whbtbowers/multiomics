library(heatmap3)
library(RColorBrewer)
library(gplots)
library(SNFtool)
library(genefilter)

ind.tb_ltbi <- c()

for (j in 1:nrow(df.meta.val)){
  if ((df.meta.val$group[j] == 1) || (df.meta.val$group[j] == 3)){
    ind.tb_ltbi <- c(ind.tb_ltbi, j)
  }
}

df.meta.val.tb_ltbi <- df.meta.val[ind.tb_ltbi,]
df.gene.val.tb_ltbi <- as.data.frame(t(df.gene.val[ind.tb_ltbi, match(sel.gene.tb_ltbi$features, colnames(df.gene.val))]))



dist.tb_ltbi <- dist(as.matrix(df.gene.val.tb_ltbi))
clust.tb_ltbi <- hclust(dist.tb_ltbi)

plot(clust.tb_ltbi)

groups <- df.meta.val.tb_ltbi$group

groupSideCols <- c()
for (i in 1:length(groups)){
  if (groups[i]==1){
    groupSideCols <- c(groupSideCols, "steelblue")
  } else {
    groupSideCols <- c(groupSideCols, "red")
  }
}



heatmap.2(as.matrix(df.gene.val.tb_ltbi),
          labCol = df.gene.val.tb_ltbi$group,
          dendrogram="row",
          ColSideColors = groupSideCols,
          trace="none",
          col=colorRampPalette(c("green", "black","red"))
          )

