library(pROC)
library(ggplot2)
library(tiger)

data("aSAH")

# Create smooth ROC curve
rocobj1 <- roc(aSAH$outcome, aSAH$s100b)#, smooth=TRUE)
rocobj2 <- roc(aSAH$outcome, aSAH$wfns)#, smooth=TRUE)
plot.roc(rocobj1)
plot.roc(rocobj2)

rand <-roc(c(0,1), c(0,1))

g <- ggroc(rocobj1)

g2 <- ggroc(list(s100b=rocobj1, wfns=rocobj2, random=rand))

multiroc <- multiclass.roc(aSAH$gos6, aSAH$s100b)
plot(multiroc)
