library(glmnet)
library(dplyr)
library(randomForest)

#samps.2.exclude <- c("z0132662","z0190864")
min.cov <- 100
description <- paste0("min.cov.",min.cov)
training.min.CR <- 3

load(paste0("data/corrected.pct.meth.",description,".Rdata"))
load("data/age.data.rda")
load("data/sites.and.inds.from.Eric.rda")
age.data <- filter(age.data, swfsc.labid %in% ids.to.keep)
sites.2.keep <- gsub("_00", replacement = ".", sites.to.keep)
sites.2.keep <- gsub("_0", replacement = ".", sites.2.keep)
sites.2.keep <- gsub("_", replacement = ".", sites.2.keep)
corrected.pct.meth <- select(corrected.pct.meth, sites.2.keep)

site.names <- names(corrected.pct.meth)
corrected.pct.meth <- cbind(id=rownames(corrected.pct.meth),corrected.pct.meth)
names(age.data)[which(names(age.data) %in% c("swfsc.labid","age.best","age.confidence"))] <- c("id","age.point","confidence")
age.data$numeric.sex <- 0
age.data$numeric.sex[which(age.data$sex=="Female")] <- 1
age.data$numeric.sex[which(age.data$sex=="Male")] <- 0
age.data$decade <- floor(age.data$age.point/10)
age.data$decade[which(age.data$decade>3)] <- 3
first.meth.col <- dim(age.data)[2] + 1

all.samples <- left_join(age.data,corrected.pct.meth)
calibration.set <- subset(all.samples, subset=(all.samples$confidence >= training.min.CR))

################################################################################
# Pearson rank correlations

sex.corr.coeff <- do.call('rbind',lapply(first.meth.col:length(all.samples), function(site){
  temp <- cor.test(all.samples[[site]],all.samples$numeric.sex,method="pearson", use="na.or.complete")
  return(c(corr.coeff=temp$estimate, p.val = temp$p.value))
}))
rownames(sex.corr.coeff) <- site.names

age.corr.coeff <- do.call('rbind',lapply(first.meth.col:length(calibration.set), function(site){
  temp <- cor.test(calibration.set[[site]],calibration.set$age.point,method="pearson", use="na.or.complete")
  return(c(corr.coeff=temp$estimate, p.val = temp$p.value))
}))
rownames(age.corr.coeff) <- site.names

################################################################################
# glmnet analyses
################################################################################

alpha.vals <- seq(from = 0.1, to = 0.9, by = 0.1)

#alternately remove sites and samples to eliminate all NAs

calibration.set.complete <- calibration.set

NAs.by.sample <- sapply(1:nrow(calibration.set.complete), function(s){
  length(which(is.na(calibration.set.complete[s,first.meth.col:ncol(calibration.set.complete)])))
})

calibration.set.complete <- calibration.set.complete[-which(NAs.by.sample>0),]

#################################################
# Cross-validated linear age fit of age over a range of alphas

x.meth <- as.matrix(select(calibration.set.complete, matches(site.names)))
y.age <- as.matrix(calibration.set.complete$age.point)

test.alpha <- lapply(alpha.vals, function(a){
  cvfit <- cv.glmnet(x.meth,y.age, alpha = a)
  corr.coef <- coef(cvfit, s = "lambda.min")
  predicted.age <- predict(cvfit, x.meth, s = "lambda.min")
  age.error <- y.age - predicted.age
  r2 <- summary(lm(predicted.age ~ y.age))$adj.r.squared
  num.predictors <- length(which(corr.coef != 0)) - 1 #don't count the y-intercept
  return(list(corr.coef=corr.coef,age.error=age.error, predicted.age=predicted.age,
              mean.age.error=mean(abs(age.error)), median.age.error = 
                median(abs(age.error)), r2 = r2, num.predictors = num.predictors))
})

age.error.sum <- data.frame(t(sapply(test.alpha, function(a){
  c(a$mean.age.error, a$median.age.error, a$r2, a$num.predictors)
})))
names(age.error.sum) <- c("Mean","Median","Rsquared","num.predictors")
age.error.sum <- cbind("alpha" = seq(.1, .9, by = .1), age.error.sum)
write.csv(age.error.sum, file = paste0("results/age.error.by.alpha.minCR", training.min.CR, ".csv"))

best.alpha.index <- which(age.error.sum$Median == min(age.error.sum$Median))
#best.alpha.index <- 3
best.alpha <- alpha.vals[best.alpha.index]
age.sum <- cbind(calibration.set.complete[,1:13], test.alpha[[best.alpha.index]]$predicted.age)
names(age.sum)[14] <- "predicted.age"
plot(age.sum$age.point, age.sum$predicted.age, pch=(age.sum$numeric.sex+24), 
     col = c("red","green","blue"), bg= c("red","green","blue"), 
     xlim = c(0,45), ylim = c(0,45))
lines(c(0,50),c(0,50))

model.coeff <- test.alpha[[best.alpha.index]]$corr.coef
model.coeff <- data.frame(model.coeff[which(model.coeff != 0),])
write.csv(model.coeff, file = paste0("results/model.coefficients.bestalpha.minCR", training.min.CR, ".csv"))

save(test.alpha, calibration.set.complete,age.error.sum,best.alpha,age.sum, site.names,
     file = paste0("results/glmnet.absolute.age.minCR", training.min.CR, ".results.rda"))

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

