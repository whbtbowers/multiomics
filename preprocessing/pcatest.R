
library(stats)
library(ggfortify)
stats::prcomp

df <- iris[c(1, 2, 3, 4)]
#autoplot(prcomp(df))

autoplot(prcomp(df), data = iris, colour = 'Species')

#autoplot(prcomp(df), data = iris, colour = 'Species', label = TRUE, label.size = 3)

#autoplot(prcomp(df), data = iris, colour = 'Species', shape = FALSE, label.size = 3)
