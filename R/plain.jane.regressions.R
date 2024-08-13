library(glmnet)
library(dplyr)
library(tidyverse)
library(randomForest)
library(rfPermute)
library(e1071)
source("R/functions/LOOV.funcs.R")
source("R/functions/plotting.funcs.R")
source("R/data.prep/combine.age.and.meth.data.R")
load("data/color.palettes.rda")

description <- "logit_mincov100" 
age.transform <- "ln" # "none" "ln" "sqrt"
dat <- combine.age.and.meth.data(description)

weight.type <- "none"
training.min.CR <- 4

site.names <- dat$site.names
all.samples <- dat$all.samples
if (weight.type == "linear") all.samples$wt <- all.samples$age.confidence/5
if (weight.type == "none") all.samples$wt <- 1

calibration.set.complete <- subset(all.samples, subset=(all.samples$age.confidence >= training.min.CR))

#alternately remove sites and samples to eliminate all NAs
NAs.by.sample <- sapply(1:nrow(calibration.set.complete), function(s){
  length(which(is.na(calibration.set.complete[s,dat$first.meth.col:ncol(calibration.set.complete)])))
})

calibration.set.complete <- calibration.set.complete[which(NAs.by.sample == 0),]

################################################################################
# Elastic Net Regression

x.meth <- as.matrix(select(calibration.set.complete, matches(site.names)))
y.age <- as.matrix(calibration.set.complete$age.best)
if(age.transform == "sqrt") y.age <- sqrt(y.age + 1)
if(age.transform == "ln") y.age <- log(y.age)
alpha.vals <- seq(from = 0.1, to = 0.9, by = 0.1)

test.alpha <- lapply(alpha.vals, function(a){
  cvfit <- cv.glmnet(x.meth, y.age, alpha = a)
  corr.coef <- coef(cvfit, s = "lambda.min")
  predicted.age <- predict(cvfit, x.meth, s = "lambda.min")
  age.error <- y.age - predicted.age
  r2 <- summary(lm(predicted.age ~ y.age))$adj.r.squared
  if(age.transform == "sqrt") {
    predicted.age <- predicted.age^2 - 1
    age.error <- (y.age^2 - 1) - predicted.age
  }
  if(age.transform == "ln") {
    predicted.age <- exp(predicted.age)
    age.error <- exp(y.age) - predicted.age
  }
  num.predictors <- length(which(corr.coef != 0)) - 1 #don't count the y-intercept
  return(list(cvfit=cvfit, corr.coef=corr.coef,age.error=age.error, predicted.age=predicted.age,
              mean.age.error=mean(abs(age.error)), median.age.error = 
                median(abs(age.error)), r2 = r2, num.predictors = num.predictors))
})

age.error.sum <- data.frame(t(sapply(test.alpha, function(a){
  c(Mean = a$mean.age.error, Median = a$median.age.error, Rsquared = a$r2, num.predictors = a$num.predictors)
})))
age.error.sum <- cbind("alpha" = seq(.1, .9, by = .1), age.error.sum) %>% arrange(Median)

best.alpha <- age.error.sum$alpha[1]

loov.res <- glmnet.loov(calibration.set.complete, site.names = site.names, alpha = best.alpha, age.transform = age.transform)
loov.res <- left_join(loov.res, select(calibration.set.complete, c(id,sex,age.confidence, wt)), by = "id")
loov.error.sum <- c(mean = mean(abs(loov.res$error)), median = median(abs(loov.res$error)))

save(test.alpha, calibration.set.complete, best.alpha, site.names,loov.res, loov.error.sum, 
     file = paste0("results/glmnet.LOOV.plainjane-", age.transform, ".results.rda"))

# repeat loov with alpha = 0.5
loov.res <- glmnet.loov(calibration.set.complete, site.names = site.names, alpha = 0.5)
loov.res <- left_join(loov.res, select(calibration.set.complete, c(id,sex,age.confidence, wt)), by = "id")
loov.error.sum <- c(mean = mean(abs(loov.res$error)), median = median(abs(loov.res$error)))

save(test.alpha, calibration.set.complete, best.alpha, site.names,loov.res, loov.error.sum, 
     file = paste0("results/glmnet.LOOV.plainjane.alphaHalf-", age.transform, ".results.rda"))

################################################################################
# Random Forest Regression

rf.data.age <- select(calibration.set.complete, c(age.best, all_of(dat$site.names)))
wt <- calibration.set.complete$wt
if(age.transform == "sqrt") rf.data.age$age.best <- sqrt(rf.data.age$age.best + 1)
if(age.transform == "ln") rf.data.age$age.best <- log(rf.data.age$age.best)

date()
rf.age <- rfPermute(
  age.best ~ .,
  rf.data.age,
  ntree = 2000,
  num.rep = 2,
  weights = wt
)
date()
if(age.transform == "sqrt") rf.age$rf$predicted <- rf.age$rf$predicted^2 - 1
if(age.transform == "ln") rf.age$rf$predicted <- exp(rf.age$rf$predicted)
rf.res <- bind_cols(calibration.set.complete, predicted.age = rf.age$rf$predicted) %>% 
  select(c(id, age.best, wt, predicted.age)) %>% 
  left_join(select(all.samples, c(id, sex, age.confidence))) %>% 
  mutate(error = predicted.age - age.best)

save(rf.age, rf.res, file = paste0("results/rf.regression-plainjane-", age.transform, ".rda"))

################################################################################
# SVM

data.age <- select(calibration.set.complete, c(id, age.best, all_of(dat$site.names)))
#data.age$untransformed.age.best <- data.age$age.best
if(age.transform == "sqrt") data.age$age.best <- sqrt(data.age$age.best + 1)
if(age.transform == "ln") data.age$age.best <- log(data.age$age.best)

date()
tune.obj <- tune(svm, 
                 age.best ~ .,
                 data = select(data.age, -id),
                 ranges = list(
                   cost = 10^(seq(-4, 5, 0.1)),
                   gamma = 10^(seq(-5, 4, 0.1))),
                 tunecontrol = tune.control(sampling = "cross"),
                 cross = 10)
date()

loov.age <- svm.loov(
  dat = data.age,
  cost = tune.obj$best.parameters$cost,
  gamma = tune.obj$best.parameters$gamma,
  age.transform = age.transform
) 
loov.age <- left_join(loov.age, select(all.samples, c(id,sex,age.confidence, wt)), by = "id")

loov.error.sum <- c(mean = mean(abs(loov.age$error)), median = median(abs(loov.age$error)))

save(tune.obj, loov.age, file = paste0("results/svm.regression-plainjane-", age.transform, ".rda"))

