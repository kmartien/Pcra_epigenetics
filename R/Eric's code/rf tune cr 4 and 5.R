rm(list = ls())
library(tidyverse)
library(randomForest)
library(rfPermute)

load("../age_and_methylation_data.rdata")

model.ages <- age.df |>
  filter(swfsc.id %in% ids.to.keep) |>
  arrange(age.confidence, age.best, age.min, age.max) |> 
  left_join(
    logit.meth.normal.params |>
      filter(loc.site %in% sites.to.keep) |>
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  ) |> 
  filter(age.confidence >= 4) |> 
  column_to_rownames('swfsc.id')

params <- expand.grid(sampsize = 35:49, mtry = 30:150, KEEP = FALSE)
param.df <- do.call(rbind, parallel::mclapply(1:nrow(params), function(i){
  x <- randomForest(
    y = model.ages$age.best,
    x = model.ages[, sites.to.keep],
    mtry = params$mtry[i],
    ntree = 10000,
    sampsize = params$sampsize[i],
    replace = FALSE
  )
  c( 
    sampsize = params$sampsize[i],
    mtry = params$mtry[i],
    mean.mse = x$mse[length(x$mse)],
    pct.var = 100 * x$rsq[length(x$rsq)]
  )
}, mc.cores = 5)) |> 
  as.data.frame()

param.df |> 
  group_by(sampsize, mtry) |> 
  summarize(pct.var = mean(pct.var), .groups = 'drop') |> 
  ggplot() +
  geom_tile(aes(sampsize, mtry, fill = pct.var)) +
  scale_fill_distiller(palette = 'YlOrRd', direction = 1)




rp <- rfPermute(
  y = model.ages$age.best,
  x = model.ages[, sites.to.keep],
  mtry = 85,
  ntree = 5000,
  importance = TRUE,
  sampsize = 48,
  replace = FALSE,
  num.rep = 1000,
  num.cores = 10
)
print(rp)
plotTrace(rp)

plotInbag(rp)

plotImportance(rp)
plotImportance(rp, sig = TRUE)

model.ages |> 
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


test.ages <- age.df |>
  filter(swfsc.id %in% ids.to.keep) |>
  arrange(age.confidence, age.best, age.min, age.max) |> 
  left_join(
    logit.meth.normal.params |>
      filter(loc.site %in% sites.to.keep) |>
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  ) |> 
  filter(age.confidence < 4) |> 
  column_to_rownames('swfsc.id')


test.ages$age.est <- predict(rp$rf, test.ages[, sites.to.keep])

model.ages |> 
  mutate(age.est = rp$rf$predicted) |> 
  bind_rows(test.ages) |> 
  mutate(age.confidence = factor(age.confidence)) |> 
  ggplot() +
  geom_abline(intercept = 0, slope = 1) +
  geom_point(aes(age.best, age.est, color = age.confidence), size = 3) +
  geom_segment(aes(x = age.min, xend = age.max, y = age.est, yend = age.est, color = age.confidence), alpha = 0.5) +
  scale_color_manual(values = conf.colors) +
  theme(legend.position = 'top')






rp |> 
  importance() |> 
  as.data.frame() |> 
  rownames_to_column('loc.site') |> 
  rename(pval = '%IncMSE.pval') |> 
  filter(pval <= 0.05) |> 
  pull('loc.site') |> 
  saveRDS('rf_imp_sites_cr_4_and_5.rds')

save.image('random_forest_cr_4_and_5.rdata')

