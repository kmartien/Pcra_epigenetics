library(dplyr)
library(ggplot2)
library(gridExtra)
library(viridis)
source("R/functions/plotting.funcs.R")

load("results/glmnet.sex.cluster.rda")
load("data/age.data.rda")
load("data/color.palettes.rda")

min.cov <- 100
meth.type <- "logit" # choose "pct", "pct.no.zero", or "logit"
weighted <- TRUE
description <- paste0(meth.type, ".min.cov.",min.cov)
training.min.CR <- 3
if (weighted) description <- paste0("weighted.", description)

# Age distribution of training data set
jpeg(file = paste0("results/glmnet.age.distribution.histogram.minCR",training.min.CR,".jpg"))
p.age.distribution(age.data, min.CR = training.min.CR)
dev.off()

# Bar chart of social cluster membership by sex
clust.dat <- data.complete

jpeg(file = paste0("results/glmnet.cluster.and.sex.distribution.histogram.jpg"))
p.cluster.distribution(clust.dat)
dev.off()

# Regression plots of glmnet cv.fit results and LOOCV results
#cv.training.regression <- lm(predicted.age~age.point, data = filter(age.sum, confidence >= training.min.CR))
plots <- lapply(2:4, function(x){
  description <- paste0(meth.type, ".min.cov.",min.cov)
  if (weighted) description <- paste0(".weighted.", description)
  load(paste0("results/glmnet.absolute.age.minCR", x, description, ".results.rda"))
  load(paste0("results/glmnet.LOOV.minCR", x, description, ".results.rda"))
  names(loov.res)[which(names(loov.res) %in% c("age.point","confidence"))] <- c("age.best","age.confidence")
  plot.loov.res(loov.res, min.CR = x)
})

jpeg(file = "results-raw/glmnet.regression.plots.compare.minCR.jpg", width = 500, height = 1200)
plots$ncol <- 1
do.call(grid.arrange,plots)
dev.off()

# Distribution of deviations for glmnet
load(paste0("results/glmnet.absolute.age.minCR", training.min.CR, ".results.rda"))
dat <- age.sum
names(dat)[which(names(dat) %in% c("age.point","confidence"))] <- c("age.best","age.confidence")

jpeg(file = "results/deviation.plot.glmnet.jpg", width = 800, height = 600)
plot.deviation(dat, training.min.CR)
dev.off()

# Histogram of age errors from LOOCV results
load(paste0("results/glmnet.LOOV.minCR", training.min.CR, ".results.rda"))
names(loov.res)[which(names(loov.res) %in% c("age.point","confidence"))] <- c("age.best","age.confidence")
jpeg(file = paste0("results/loov.error.histogram.minCR", training.min.CR, ".jpg"))
loov.hist(loov.res)
dev.off()


