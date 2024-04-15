library(glmnet)
library(dplyr)
library(randomForest)

min.cov <- 100
meth.type <- "pct.no.zero" # "logit" or "pct"
weighted <- TRUE
description <- paste0(meth.type, ".min.cov.",min.cov)
label <- paste0("model.6.", meth.type) #"CR4.5

dat <- combine.age.and.meth.data(description)
all.samples <- dat.list$all.samples
if (weighted) description <- paste0("weighted.", description)
if(!weighted) all.samples$wt <- 1

calibration.set <- subset(all.samples, subset=(all.samples$type == "Training"))

################################################################################
# glmnet analyses
################################################################################

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
y.age <- as.matrix(calibration.set.complete$age.point)

test.alpha <- lapply(alpha.vals, function(a){
  cvfit <- cv.glmnet(x.meth,y.age, alpha = a)
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
write.csv(age.error.sum, file = paste0("results/age.error.by.alpha.half.", label, ".csv"))

best.alpha.index <- which(age.error.sum$Median == min(age.error.sum$Median))
best.alpha <- alpha.vals[best.alpha.index]
pred.age <- predict(test.alpha[[best.alpha.index]]$cvfit, as.matrix(select(all.samples, matches(dat$site.names))), s = "lambda.min")
age.sum <- cbind(all.samples[,1:13], pred.age)
names(age.sum)[14] <- "predicted.age"

save(test.alpha, calibration.set.complete,age.error.sum,best.alpha,age.sum, dat$site.names,
     file = paste0("results/glmnet.absolute.age.", label, ".results.rda"))

#################################################
# Cross-validated multinomial fit of decade over a range of alphas

y.decade <- as.matrix(calibration.set.complete$decade)

test.alpha.decade <- lapply(alpha.vals, function(a){
  cvfit <- cv.glmnet(x.meth,y.decade, alpha = a, family="multinomial")
  corr.coef <- coef(cvfit, s = "lambda.min")
  predicted.decade <- predict(cvfit, x.meth, s = "lambda.min", type="class")
  decade.error <- y.decade - as.numeric(predicted.decade)
  return(list(corr.coef=corr.coef,decade.error=decade.error, 
              mean.decade.error=mean(abs(decade.error)), median.decade.error = median(abs(decade.error))))
})

decade.error.sum <- t(sapply(test.alpha.decade, function(a){
  c(a$mean.decade.error, a$median.decade.error)
}))

