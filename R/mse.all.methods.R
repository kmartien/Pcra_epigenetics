library(dplyr)
library(tidyverse)
library(ggplot2)
library(gridExtra)
source("R/functions/plotting.funcs.R")
#load("data/model.params.rda")
load("data/color.palettes.rda")

age.transform <- "ln"
training.min.CR <- 2
#model <- 30
#params <- model.params[model,]
#description <- paste0("model", params$model, "-logit_mincov100")

dat <- list()
load(paste0("results/glmnet.LOOV.plainjane-", age.transform, "-minCR", training.min.CR, ".results.rda"))
dat$ENR <- loov.res
load(paste0("results/glmnet.LOOV.plainjane.alphaHalf-", age.transform, "-minCR", training.min.CR, ".results.rda"))
dat$ENR0.5 <- loov.res
load(paste0("results/rf.regression-plainjane-", age.transform, "-minCR", training.min.CR, ".rda"))
dat$RF <- rf.res
load(paste0("results/svm.regression-plainjane-", age.transform, "-minCR", training.min.CR, ".rda"))
dat$SVM <- loov.age

mse <- do.call(rbind, lapply(1:length(dat), function(i){
  bind_cols(method = names(dat)[i],
            bind_rows(
              (dat[[i]] %>% 
                 group_by(age.confidence) %>% 
                 mutate(sq.error = error^2) %>% 
                 summarise(mse = sum(sq.error) / n())),
              (dat[[i]] %>% 
                 mutate(sq.error = error^2) %>% 
                 summarise(mse = sum(sq.error) / n()))
            )
  )
}))
save(mse, file = paste0("results/mse_by_confidence-all_methods-minCR", training.min.CR, ".rda"))
write.csv(mse, file = paste0("results-raw/mse_by_confidence-all_methods-minCR", training.min.CR, ".csv"))
