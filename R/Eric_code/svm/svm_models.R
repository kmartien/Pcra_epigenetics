library(tidyverse)
library(mgcv)
library(e1071)
source('R/Eric_code/misc_funcs.R')
load("R/Eric_code/age_and_methylation_data.rdata")

sites.2.use <- "RFsites" #"All" or "RFsites"

sites <- sites.to.keep
if(sites.2.use == "RFsites"){
  # select important sites from Random Forest
  sites <- readRDS('R/Eric_code/random forest/rf_site_importance.rds') |> 
    filter(pval <= 0.1) |> 
    pull('loc.site')
}

age.df <- age.df |>  
  filter(swfsc.id %in% ids.to.keep)

model.df <- age.df |> 
  left_join(
    logit.meth.normal.params |> 
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  ) 

nrep <- 1000
ncores <- 10

# Best age and methylation estimates --------------------------------------

train.df <- filter(model.df, age.confidence %in% 4:5)
# tune.obj <- tune(svm, 
#                  age.best ~ .,
#                  data = select(train.df, c(age.best, all_of(sites))),
#                  ranges = list(
#                    cost = 10^(seq(-4, 5, 0.1)),
#                    gamma = 10^(seq(-5, 4, 0.1))),
#                  tunecontrol = tune.control(sampling = "cross"),
#                  cross = 10)
#save(tune.obj, file = paste0("R/Eric_code/svm/tune.obj.", sites.2.use, ".rda"))
load(paste0("R/Eric_code/svm/tune.obj.", sites.2.use, ".rda"))

predictAllIDsSVM(train.df, model.df, sites, 'age.best', tune.obj) |> 
  saveRDS(paste0('R/Eric_code/svm/svm_best_', sites.2.use, '.rds'))


# Random age and best methylation estimates -------------------------------

parallel::mclapply(1:nrep, function(j) {
  # random sample of ages and methylation - only use random age
  ran.df <- model.df |> 
    left_join(
      sampleAgeMeth(age.df, logit.meth.normal.params) |> 
        select(swfsc.id, age.ran),
      by = 'swfsc.id'
    )
  
  train.df <- filter(ran.df, age.confidence %in% 4:5)
  predictAllIDsSVM(train.df, ran.df, sites, 'age.ran', tune.obj)
}, mc.cores = ncores) |> 
  bind_rows() |> 
  saveRDS(paste0('R/Eric_code/svm/svm_ran_age_', sites.2.use, '.rds'))


# Random age and random methylation estimates -----------------------------

parallel::mclapply(1:nrep, function(j) {
  # random sample of ages and methylation
  ran.df <- sampleAgeMeth(age.df, logit.meth.normal.params) 
  
  train.df <- filter(ran.df, age.confidence %in% 4:5)
  predictAllIDsSVM(train.df, ran.df, sites, 'age.ran', tune.obj)
}, mc.cores = ncores) |> 
  bind_rows() |> 
  saveRDS(paste0('R/Eric_code/svm/svm_ran_age_meth_', sites.2.use, '.rds'))
