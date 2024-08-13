library(dplyr)
load("data/model.params.rda")

svm.age.res.sum <- lapply(1:nrow(model.params), function(i){
  params <- model.params[i,]
  description <- paste0("model", params$model, "-logit_mincov100")
  load(file = paste0("results/svm.regression-", description, ".rda"))
  #lapply(2:4, function(min.cr){
  #load(file = paste0("results/svm.results.minCR", min.cr, ".", description, ".rda"))
  x <- filter(loov.age, age.confidence >= 4)
  mae.4.5 <- summarise(x, median(abs(error)))
  cor.coeff.4.5 <- round(cor.test(x$age.best, x$predicted.age, method = "pearson")$estimate,2)
  return(list(best.params = tune.obj$best.parameters, mae = rbind(mae, c(4.5, as.numeric(mae.4.5), cor.coeff.4.5)), loov.age = loov.age))
})

mae.sum <- do.call(rbind, lapply(1:nrow(model.params), function(i){
  cbind(i, svm.age.res.sum[[i]]$mae)
})) %>% data.frame()
mae.sum$CR <- as.character(mae.sum$CR)
mae.sum$CR[which(mae.sum$CR == "4.5")] <- "4&5"
write.csv(mae.sum, file = paste0("results-raw/svm.age.mae.summary.csv"), row.names = FALSE)

loov.age.error.sum <- do.call(rbind, lapply(1:nrow(model.params), function(i){
  cbind(model = model.params$model[i], svm.age.res.sum[[i]]$loov.age)
  })) %>% data.frame()
write.csv(loov.age.error.sum, file = "results-raw/svm.age.loov.error.summary.csv", row.names = FALSE)

param.sum <- do.call(rbind, lapply(1:nrow(model.params), function(i){
  cbind(i, svm.age.res.sum[[i]]$best.params)
}))
write.csv(param.sum, file = paste0("results-raw/svm.age.best.param.summary.csv"), row.names = FALSE)
