library(glmnet)
library(dplyr)
library(tidyverse)
library(randomForest)
source("R/functions/glmnet.LOOV.R")
source("R/functions/plotting.funcs.R")
source("R/data.prep/combine.age.and.meth.data.R")

min.cov <- 100
meth.type <- "pct" # choose "pct", "pct.no.zero", or "logit"
weighted <- TRUE
description <- paste0(meth.type, ".min.cov.",min.cov)
training.min.CR <- 2

dat <- combine.age.and.meth.data(description)
all.samples <- dat$all.samples
if (weighted) description <- paste0("weighted.", description)
if(!weighted) all.samples$wt <- 1

calibration.set <- subset(all.samples, subset=(all.samples$age.confidence >= training.min.CR))

alpha.vals <- seq(from = 0.1, to = 0.9, by = 0.1)

#alternately remove sites and samples to eliminate all NAs

calibration.set.complete <- calibration.set

NAs.by.sample <- sapply(1:nrow(calibration.set.complete), function(s){
  length(which(is.na(calibration.set.complete[s,dat$first.meth.col:ncol(calibration.set.complete)])))
})

calibration.set.complete <- calibration.set.complete[which(NAs.by.sample == 0),]

#################################################
# Cross-validated linear age fit of age over a range of alphas

x.meth <- as.matrix(select(calibration.set.complete, matches(dat$site.names)))
y.age <- as.matrix(calibration.set.complete$age.best)
wt <- calibration.set.complete$wt

test.alpha <- lapply(alpha.vals, function(a){
  cvfit <- cv.glmnet(x.meth,y.age, weights = wt, alpha = a)
  corr.coef <- coef(cvfit, s = "lambda.min")
  predicted.age <- predict(cvfit, x.meth, s = "lambda.min")
  age.error <- y.age - predicted.age
  r2 <- summary(lm(predicted.age ~ y.age))$adj.r.squared
  num.predictors <- length(which(corr.coef != 0)) - 1 #don't count the y-intercept
  return(list(cvfit=cvfit, corr.coef=corr.coef,age.error=age.error, predicted.age=predicted.age,
              mean.age.error=mean(abs(age.error)), median.age.error = 
                median(abs(age.error)), r2 = r2, num.predictors = num.predictors))
})

age.error.sum <- data.frame(t(sapply(test.alpha, function(a){
  c(a$mean.age.error, a$median.age.error, a$r2, a$num.predictors)
})))
names(age.error.sum) <- c("Mean","Median","Rsquared","num.predictors")
age.error.sum <- cbind("alpha" = seq(.1, .9, by = .1), age.error.sum)
#write.csv(age.error.sum, file = paste0("results-raw/age.error.by.alpha.", description, ".csv"))

best.alpha.index <- which(age.error.sum$Median == min(age.error.sum$Median))
best.alpha <- alpha.vals[best.alpha.index]
pred.age <- predict(test.alpha[[best.alpha.index]]$cvfit, as.matrix(select(all.samples, matches(site.names))), s = "lambda.min")
age.sum <- cbind(all.samples[,1:13], pred.age)
names(age.sum)[14] <- "predicted.age"
#plot(age.sum$age.best, age.sum$predicted.age, pch=(age.sum$numeric.sex+24), 
#     col = c("red","green","blue", "black"), bg= c("red","green","blue", "black"), 
#     xlim = c(0,45), ylim = c(0,45))
#lines(c(0,50),c(0,50))

model.coeff <- test.alpha[[best.alpha.index]]$corr.coef
model.coeff <- data.frame(model.coeff[which(model.coeff != 0),])
#write.csv(model.coeff, file = paste0("results/model.coefficients.bestalpha.", description, ".csv"))

save(test.alpha, calibration.set.complete,age.error.sum,best.alpha,age.sum, site.names,
     file = paste0("results/glmnet.absolute.age.", description, ".results.rda"))

dat <- select(calibration.set.complete, c(id, age.best, matches(site.names)))
loov.res <- glmnet.loov(dat, alpha = best.alpha)
loov.res <- left_join(loov.res, select(age.sum, c(id,sex,age.confidence)), by = "id")

loov.error.sum <- c(mean = mean(abs(loov.res$error)), median = median(abs(loov.res$error)))
#loov.regression <- lm(predicted.age~age.best, data = loov.res)

save(calibration.set.complete, best.alpha, site.names,loov.res, loov.error.sum, 
     file = paste0("results/glmnet.LOOV.", description, ".results.rda"))

load("data/color.palettes.rda")
jpeg(file = paste0("results-raw/glmnet.regression.plots.minCR.",training.min.CR, ".", description,".jpg"), width = 500, height = 400)
plot.loov.res(loov.res, min.CR = training.min.CR)
dev.off()
