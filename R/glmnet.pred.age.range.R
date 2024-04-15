library(glmnet)
library(dplyr)
library(tidyverse)

min.cov <- 100
meth.type <- "logit" # "logit", "pct", or "pct.no.zero"
description <- paste0(meth.type, ".min.cov.",min.cov)

load(paste0("results/glmnet.stability.results.", description, ".rda"))
site.names <- names(calibration.set.complete)[15:ncol(calibration.set.complete)] 

nreps <- length(cvfit.list)
training.min.CR <- min(calibration.set.complete$confidence)

x.meth <- as.matrix(select(calibration.set.complete, matches(site.names)))
y.age <- as.matrix(calibration.set.complete$age.point)

best.alpha <- filter(chosen.model, lambda.type == 1) %>% 
  filter(error.measure == "median.error") 

range.pred.age <- data.frame(do.call('rbind', lapply(1:nreps, function(i){
  a <- filter(best.alpha, iter == i) %>% select(alpha)
  cvfit <- cvfit.list[[i]][[a$alpha[1] * 10]]
  pred.age <- predict(cvfit, x.meth, s = "lambda.min")
  range(pred.age)
})))
names(range.pred.age) <- c("min", "max")

write.csv(range.pred.age, file = paste0("results/range.predicted.ages.", description, ".rda"))
