library(glmnet)
library(dplyr)
library(randomForest)
library(ggplot2)

training.min.CR <- 3
load(paste0("results/glmnet.absolute.age.minCR", training.min.CR, ".results.rda"))

loov.res <- do.call(rbind,lapply(calibration.set.complete$id, function(i){
  
  dat.loov <- calibration.set.complete[-which(calibration.set.complete$id == i),]
  x.meth <- as.matrix(select(dat.loov, matches(site.names)))
  y.age <- as.matrix(dat.loov$age.point)
  
  cvfit <- cv.glmnet(x.meth,y.age, alpha = best.alpha)
#  corr.coef <- coef(cvfit, s = "lambda.min")
  oob.meth <- as.matrix(filter(calibration.set.complete, id == i) %>%
    select(matches(site.names)))
  oob.age <- calibration.set.complete$age.point[which(calibration.set.complete$id == i)]
  predicted.age <- predict(cvfit, oob.meth, s = "lambda.min")
  return(c("age.point" = oob.age, "predicted.age" = predicted.age, "error" = (predicted.age - oob.age)))
})) %>% data.frame()

loov.res <- cbind(labid = calibration.set.complete$id, loov.res)
loov.res$confidence <- age.sum$confidence
loov.res$sex <- age.sum$sex
loov.error.sum <- c(mean = mean(abs(loov.res$error)), median = median(abs(loov.res$error)))
loov.regression <- lm(predicted.age~age.point, data = loov.res)

save(calibration.set.complete, best.alpha, site.names,loov.res, loov.error.sum, 
     file = paste0("results/glmnet.LOOV.minCR", training.min.CR, ".results.rda"))
