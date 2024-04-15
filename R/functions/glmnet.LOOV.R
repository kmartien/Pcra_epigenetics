glmnet.loov <- function(dat, alpha){
  loov.res <- do.call(rbind,lapply(dat$id, function(i){
    
    dat.loov <- filter(dat, id != i)
    x.meth <- as.matrix(select(dat.loov, -c(id, age.best, wt)))
    y.age <- as.matrix(dat.loov$age.best)
    wt <- dat.loov$wt
    
    cvfit <- cv.glmnet(x.meth,y.age, weights = wt, alpha = alpha)
    oob.meth <- as.matrix(filter(dat, id == i) %>%
                            select(-c(id, age.best, wt)))
    oob.age <- dat$age.best[which(dat$id == i)]
    predicted.age <- predict(cvfit, oob.meth, s = "lambda.min")
    return(c("age.best" = oob.age, "predicted.age" = predicted.age, "error" = (predicted.age - oob.age)))
  })) %>% data.frame()
  
  loov.res <- cbind(id = dat$id, loov.res)
}
