rm(list = ls())
library(tidyverse)
library(mgcv)
library(randomForest)
source('R/Eric_code/misc_funcs.R')
load("R/Eric_code/age_and_methylation_data.rdata")
load("R/Eric_code/random forest/rf_optim_params.rda")

sites.2.use <- "RFsites" #"All" or "RFsites"

sites <- sites.to.keep
if(sites.2.use == "RFsites"){
  # select important sites from Random Forest
  sites <- readRDS('R/Eric_code/random forest/rf_site_importance.rds') |> 
    filter(pval <= 0.1) |> 
    pull('loc.site')
}

#rf.params <- param.df |>
#  filter(mse == min(param.df$mse)) |>
#  select(c(mtry, sampsize))
rf.params <- tibble(mtry = floor(length(sites)/3), sampsize = floor(50 * 0.632)) 
 
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
predictAllIDsRF(train.df, model.df, sites, 'age.best', rf.params) |> 
  saveRDS(paste0('rf_best_', sites.2.use, '.rds'))


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
  predictAllIDsRF(train.df, ran.df, sites, 'age.ran', rf.params)
}, mc.cores = ncores) |> 
  bind_rows() |> 
  saveRDS(paste0('rf_ran_age_', sites.2.use, '.rds'))


# Random age and random methylation estimates -----------------------------

parallel::mclapply(1:nrep, function(j) {
  # random sample of ages and methylation
  ran.df <- sampleAgeMeth(age.df, logit.meth.normal.params) 
  
  train.df <- filter(ran.df, age.confidence %in% 4:5)
  predictAllIDsRF(train.df, ran.df, sites, 'age.ran', rf.params)
}, mc.cores = ncores) |> 
  bind_rows() |> 
  saveRDS(paste0('rf_ran_age_meth_', sites.2.use, '.rds'))
