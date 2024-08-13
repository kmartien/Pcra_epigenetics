library(dplyr)
library(ggplot2)
library(gridExtra)
library(viridis)
source("R/functions/plotting.funcs.R")
source("R/data.prep/combine.age.and.meth.data.R")

#load("results/glmnet.sex.cluster.rda")
load("data/color.palettes.rda")

min.cov <- 100
meth.type <- "logit" # choose "pct", "pct.no.zero", or "logit"
weighted <- TRUE
description <- paste0(meth.type, ".min.cov.",min.cov)
training.min.CR <- 1
if (weighted) description <- paste0("weighted.", description)

dat <- combine.age.and.meth.data(description)
data.complete <- select(dat$all.samples, c(id, age.best, age.confidence, social.cluster, sex, all_of(dat$site.names))) %>%
  rename(cluster.louvain = social.cluster)
NAs.by.sample <- sapply(1:nrow(data.complete), function(s){
  length(which(is.na(data.complete[s,4:ncol(data.complete)])))
})
data.complete <- data.complete[-which(NAs.by.sample>0),]

# Age distribution of training data set
jpeg(file = paste0("results-raw/glmnet.age.distribution.histogram.minCR",training.min.CR,".jpg"))
p.age.distribution(data.complete, min.CR = training.min.CR)
dev.off()

age.dist.table <- data.frame(do.call(bind_rows, lapply(2:5, function(cr){
  x <- filter(data.complete, age.confidence == cr)
  return(c(age.confidence = cr, n = nrow(x), median = median(x$age.best), min = min(x$age.best), max = max(x$age.best)))
})) %>% add_row(age.confidence = 6, n = nrow(data.complete), 
                            median = median(data.complete$age.best), min = min(data.complete$age.best), max = max(data.complete$age.best)))
age.dist.table$age.confidence <- as.character(age.dist.table$age.confidence)
age.dist.table$age.confidence[5] <- "All"
write.csv(age.dist.table, file = "results-raw/Table.age.distribution.csv", row.names = FALSE)

# Bar chart of social cluster membership by sex
jpeg(file = paste0("results-raw/cluster.and.sex.distribution.histogram.jpg"))
p.cluster.distribution(data.complete)
dev.off()

# Regression plots of glmnet cv.fit results and LOOCV results
#cv.training.regression <- lm(predicted.age~age.point, data = filter(age.sum, confidence >= training.min.CR))
loov.sum <- lapply(2:4, function(x){
  description <- paste0(meth.type, ".min.cov.",min.cov)
  if (weighted) description <- paste0(".weighted.", description)
  load(paste0("results/glmnet.absolute.age.minCR", x, description, ".results.rda"))
  load(paste0("results/glmnet.LOOV.minCR", x, description, ".results.rda"))
  names(loov.res)[which(names(loov.res) %in% c("age.point","confidence"))] <- c("age.best","age.confidence")
  plot.loov.res(loov.res, min.CR = x)
})

plots <- lapply(loov.sum, function(x){x[[1]]})
jpeg(file = paste0("results-raw/glmnet.regression.plots.", description, ".jpg"), width = 500, height = 1200)
plots$ncol <- 1
do.call(grid.arrange,plots)
dev.off()

fit.sum <- do.call(rbind, lapply(2:4, function(min.cr){
  cbind(min.cr, loov.sum[[min.cr-1]]$fit.sum)
}))
write.csv(fit.sum, file =paste0("results-raw/glmnet.fit.summary.", description, ".csv"))

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


