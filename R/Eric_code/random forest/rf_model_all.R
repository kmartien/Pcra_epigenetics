rm(list = ls())
library(tidyverse)
library(randomForest)
library(rfPermute)

load("R/Eric_code/age_and_methylation_data.rdata")

model.df <- age.df |>
  filter(swfsc.id %in% ids.to.keep) |>
  left_join(
    logit.meth.normal.params |>
      filter(loc.site %in% sites.to.keep) |>
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  ) |> 
  column_to_rownames('swfsc.id') |> 
  mutate(inv.var = 1 / age.var)

# grid of sampsize and mtry to optimize MSE over
params <- expand.grid(
  sampsize = seq(30, 65, by = 1), 
  mtry = seq(40, 70, by = 1), 
  KEEP.OUT.ATTRS = FALSE
)

#' run inverse variance weighted randomForest over sampsize and mtry grid 
#' and report deviance stats
param.df <- parallel::mclapply(1:nrow(params), function(i){
  rf <- randomForest(
    y = model.df$age.best,
    x = model.df[, sites.to.keep],
    mtry = params$mtry[i],
    ntree = 10000,
    weights = model.df$inv.var,
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
save(param.df, file = "R/Eric_code/random forest/RF.grid.search.rda")

# plot heatmap of MSE across grid
param.df |> 
  ggplot() +
  geom_tile(aes(sampsize, mtry, fill = mse)) +
  scale_fill_distiller(palette = 'RdYlBu', direction = 1)


# run inverse variance weighted rfPermute at sampsize and mtry with minimum MSE
rp <- rfPermute(
  y = model.df$age.best,
  x = model.df[, sites.to.keep],
  mtry = 55,
  ntree = 20000,
  weights = model.df$inv.var,
  sampsize = 46,
  importance = TRUE,
  replace = FALSE,
  num.rep = 1000,
  num.cores = 10
)

save.image('rf_model_all.rdata')


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
  saveRDS('rf_site_importance.rds')



# model diagnostics
print(rp)
plotTrace(rp)
plotImportance(rp, imp = '%IncMSE', alpha = 0.2)
plotImportance(rp, sig = TRUE)

model.df |> 
  mutate(
    age.est = rp$rf$predicted,
    age.confidence = factor(age.confidence)
  ) |> 
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(age.best, age.est, color = age.confidence), size = 3) +
  geom_segment(aes(x = age.min, xend = age.max, y = age.est, yend = age.est, color = age.confidence), alpha = 0.5) +
  scale_color_manual(values = conf.colors) +
  theme(legend.position = 'top')

age.pred <- rp$rf$predicted |> 
  enframe(name = 'swfsc.id', value = 'age.pred') |> 
  left_join(age.df, by = 'swfsc.id') |> 
  mutate(mse = (age.best - age.pred) ^ 2)

age.pred |> 
  mutate(age.confidence = as.character(age.confidence)) |> 
  group_by(age.confidence) |> 
  summarize(mean.mse = mean(mse), .groups = 'drop') |> 
  bind_rows(
    data.frame(age.confidence = 'All', mean.mse = mean(age.pred$mse))
  )