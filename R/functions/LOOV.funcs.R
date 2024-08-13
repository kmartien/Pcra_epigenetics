glmnet.loov <- function(dat, site.names, alpha, age.transform = "none"){
  loov.res <- do.call(rbind,lapply(dat$id, function(i){
    
    dat.loov <- filter(dat, id != i)
    x.meth <- as.matrix(select(dat.loov, matches(site.names)))
    y.age <- as.matrix(dat.loov$age.best)
    if(age.transform == "sqrt") y.age <- sqrt(y.age + 1)
    if(age.transform == "ln") y.age <- log(y.age)
    wt <- dat.loov$wt
    
    cvfit <- cv.glmnet(x.meth,y.age, weights = wt, alpha = alpha)
    oob.meth <- as.matrix(filter(dat, id == i) %>%
                            select(matches(site.names)))
    oob.age <- dat$age.best[which(dat$id == i)]
    predicted.age <- predict(cvfit, oob.meth, s = "lambda.min")
    if(age.transform == "sqrt") predicted.age <- predicted.age ^2 - 1
    if(age.transform == "ln") predicted.age <- exp(predicted.age)
    return(c("age.best" = oob.age, "predicted.age" = predicted.age, "error" = (predicted.age - oob.age)))
  })) %>% data.frame()
  
  loov.res <- cbind(id = dat$id, loov.res)
}

svm.loov <- function(dat, cost, gamma, age.transform = "none"){
  loov.res <- do.call(rbind,lapply(dat$id, function(i){
    
    dat.loov <- filter(dat, id != i) %>% select(-id)
    svm.age <- svm(
      age.best ~ ., data = dat.loov, cost = cost, gamma = gamma)
    oob.meth <- as.matrix(filter(dat, id == i) %>%
                            select(-c(id, age.best)))
    oob.age <- dat$age.best[which(dat$id == i)]
    predicted.age <- unname(predict(svm.age, oob.meth))
    if(age.transform == "sqrt") {
#      oob.age <- dat$untransformed.age.best[which(dat$id == i)]
      oob.age <- as.integer(round(oob.age^2 - 1))
      predicted.age <- predicted.age^2 - 1
    }
    if(age.transform == "ln") {
      predicted.age <- exp(predicted.age)
#      oob.age <- dat$untransformed.age.best[which(dat$id == i)]
      oob.age <- as.integer(round(exp(oob.age)))
    }
    return(c("age.best" = oob.age, "predicted.age" = predicted.age, "error" = (predicted.age - oob.age)))
  })) %>% data.frame()
  
  loov.res <- cbind(id = dat$id, loov.res)
}
