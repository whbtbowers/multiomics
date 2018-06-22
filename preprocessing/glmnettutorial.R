library(glmnet)

load("QuickStartExample.RData")

fit = glmnet(x, y)
