rm(list = ls())
library(tidyverse)
library(randomForest)
library(rfPermute)

load("R/Eric_code/age_and_methylation_data.rdata")

sites.2.use <- "All" # All or RFsites

model.df <- age.df |>
  filter(swfsc.id %in% ids.to.keep & age.confidence >= 4) |>
  left_join(
    logit.meth.normal.params |>
      filter(loc.site %in% sites.to.keep) |>
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  ) |> 
  column_to_rownames('swfsc.id') |> 
  mutate(var.wt = 1 / age.var)

params <- expand.grid(
  sampsize = seq(24, 50, by = 2), 
  mtry = seq(30, 70, by = 2), 
  KEEP.OUT.ATTRS = FALSE
)

#' run randomForest over sampsize and mtry grid 
#' and report deviance stats
param.df <- parallel::mclapply(1:nrow(params), function(i){
  rf <- randomForest(
    y = model.df$age.best,
    x = model.df[, sites.to.keep],
    mtry = params$mtry[i],
    ntree = 10000,
    sampsize = params$sampsize[i],
    replace = FALSE
  )
  data.frame( 
    sampsize = params$sampsize[i],
    mtry = params$mtry[i],
    mse = rf$mse[length(rf$mse)],
    rsq = rf$rsq[length(rf$rsq)],
    pct.var = 100 * rf$rsq[length(rf$rsq)]
  )
}, mc.cores = 6) |> 
  bind_rows() |> 
  as.data.frame()

param.df |> 
  ggplot() +
  geom_tile(aes(sampsize, mtry, fill = mse)) +
  scale_fill_distiller(palette = 'YlOrRd')

rf.params <- filter(param.df, mse = min(param.df$mse)) |>
  select(c(mtry, sapsize))
save(rf.params, file = paste0('R/Eric_code/random forest/rf_optim_params_', sites.2.use, '_CR4&5.rda'))

rp <- rfPermute(
  y = model.df$age.best,
  x = model.df[, sites.to.keep],
  mtry = rf.params$mtry,
  ntree = 10000,
  weights = model.df$var.wt,
  sampsize = rf.params$sampsize,
  importance = TRUE,
  replace = FALSE,
  num.rep = 1000,
  num.cores = 10
)

# save site importance scores and p-values
rp |> 
  importance() |> 
  as.data.frame() |> 
  rownames_to_column('loc.site') |> 
  rename(
    incMSE = '%IncMSE',
    pval = '%IncMSE.pval'
  ) |> 
  select(loc.site, incMSE, pval) |> 
  saveRDS('R/Eric_code/random forest/rf_site_importance_CR4&5.rds')

print(rp)
plotTrace(rp)
plotInbag(rp)
plotImportance(rp)
plotImportance(rp, sig = TRUE)
