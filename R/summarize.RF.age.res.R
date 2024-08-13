library(dplyr)
load("data/model.params.rda")

select.sites <-  "glmnet" #"glmnet" or "rf"

rf.age.res.sum <- lapply(1:nrow(model.params), function(i){
  params <- model.params[i,]
  description <- paste0("model", params$model, "-logit_mincov100")
  load(file = paste0("results/rf.regression-model", params$model, "-logit_mincov100.rda"))
  mae <- read.csv(file = paste0("results-raw/rf.mae-model", params$model, "-logit_mincov100.csv"))[-1]
  x <- filter(rf.res, age.confidence >= 4)
  mae.4.5 <- summarise(x, median(abs(error)))
  cor.coeff.4.5 <- round(cor.test(x$age.best, x$predicted.age, method = "pearson")$estimate,2)
  mae <- rbind(mae, c(4.5, as.numeric(mae.4.5), cor.coeff.4.5))
  return(list(mae = mae, loov.error = rf.res))
})

mae.sum <- do.call(rbind, lapply(1:nrow(model.params), function(i){
  cbind(model.params$model[i], rf.age.res.sum[[i]]$mae)
})) %>% data.frame()
mae.sum$CR <- as.character(mae.sum$CR)
mae.sum$CR[which(mae.sum$CR == "4.5")] <- "4&5"
write.csv(mae.sum, file = paste0("results-raw/rf.age.mae.summary.csv"), row.names = FALSE)

rf.age.error.sum <- do.call(rbind, lapply(1:nrow(model.params), function(i){
  cbind(model = model.params$model[i], rf.age.res.sum[[i]]$loov.error)
})) %>% data.frame()
write.csv(rf.age.error.sum, file = "results-raw/rf.age.loov.error.summary.csv", row.names = FALSE)

