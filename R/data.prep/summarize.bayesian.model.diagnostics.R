library(tidyverse)
library(dplyr)
library(modeest)

bayesian.diag.files <- as.list(list.files(path = "data/bayesian.model.diagnostics"))
names(bayesian.diag.files) <- paste("model", c(1:4,6), sep=".")

model.smry <- lapply(1:5, function(i){
  model <- names(bayesian.diag.files)[i]
  load(paste0("data/bayesian.model.diagnostics/",bayesian.diag.files[[i]]))
  pred.mode <- smry$age.pred |>
    group_by(swfsc.id) |>
    summarize(pred.mode = modeest::venter(pred.age), .groups = 'drop')
  return(right_join(distinct(select(smry$age.pred, c("swfsc.id", "type"))), pred.mode))
})
names(model.smry) <- names(bayesian.diag.files)

save(model.smry, file = "data/bayesian.models.smry.rdata")
