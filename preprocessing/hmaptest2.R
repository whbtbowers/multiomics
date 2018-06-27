nba <- read.csv("http://datasets.flowingdata.com/ppg2008.csv")

nba$Name <- with(nba, reorder(Name, PTS))

library("ggplot2")
library("reshape2")

nba.m <- melt(nba)
nba.m <- ddply(nba.m, .(variable), transform, rescale = rescale(value))

p <- ggplot(nba.m, aes(variable, Name)) + geom_tile(aes(fill = rescale),
  colour = "white") + scale_fill_gradient(low = "white",
  high = "steelblue")

base_size <- 9
p + theme_grey(base_size = base_size) + labs(x = "",
  y = "") + scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) + opts(legend.position = "none",
  axis.ticks = theme_blank(), axis.text.x = theme_text(size = base_size *
  0.8, angle = 330, hjust = 0, colour = "grey50"))
