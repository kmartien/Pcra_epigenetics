library(dplyr)
load("data/model.params.rda")

loov.age.res.sum <- lapply(1:nrow(model.params), function(i){
  params <- model.params[i,]
  description <- paste0("model", params$model, "-logit_mincov100")
  load(file = paste0("results/glmnet.LOOV.model", params$model, "-logit_mincov100.results.rda"))
  min.cr <- params$training.cr
  mae <- do.call(rbind, lapply(min.cr:5, function(cr){
    x <- filter(loov.res, age.confidence == cr)
    age.error <- mutate(x, age.error = abs(error)) %>%
      select(age.error) 
    cor.coeff <- round(cor.test(x$age.best, x$predicted.age, method = "pearson")$estimate,2)
    return(c(CR = cr, MAE = median(age.error$age.error), corr = cor.coeff))
  }))
  x <- filter(loov.res, age.confidence >= 4)
  mae.4.5 <- summarise(x, median(abs(error)))
  cor.coeff.4.5 <- round(cor.test(x$age.best, x$predicted.age, method = "pearson")$estimate,2)
  mae <- rbind(mae, c(4.5, as.numeric(mae.4.5), cor.coeff.4.5))
  return(list(mae = mae, loov.error = loov.res))
})

mae.sum <- do.call(rbind, lapply(1:nrow(model.params), function(i){
  cbind(i, loov.age.res.sum[[i]]$mae)
})) %>% data.frame()
mae.sum$CR <- as.character(mae.sum$CR)
mae.sum$CR[which(mae.sum$CR == "4.5")] <- "4&5"
write.csv(mae.sum, file = "results-raw/glmnet.age.mae.summary.csv", row.names = FALSE)

loov.age.error.sum <- do.call(rbind, lapply(1:nrow(model.params), function(i){
  cbind(model = model.params$model[i], loov.age.res.sum[[i]]$loov.error)
})) %>% data.frame()
write.csv(loov.age.error.sum, file = "results-raw/glmnet.age.loov.error.summary.csv", row.names = FALSE)

