rm(list = ls())
library(tidyverse)
library(randomForest)

load("../age_and_methylation_data.rdata")
ran.df <- readRDS('../age_meth_sample.rds')
sites <- readRDS('rf_imp_sites.rds')

age.est.rf <- parallel::mclapply(1:length(ran.df), function(i) {
  model.df <- column_to_rownames(ran.df[[i]], 'swfsc.id')
  randomForest(
    x = model.df[, sites],
    y = model.df$ran.age,
    mtry = 3,
    ntree = 20000,
    sampsize = 82,
    replace = FALSE,
    weights = model.df$confidence.wt
  )
}, mc.cores = 10) 

save(ran.df, sites, age.est.rf, file = 'rf_replicates.rdata')


pct.var <- sapply(age.est.rf, function(rf) 100 * rf$rsq[length(rf$rsq)])
hist(pct.var)
swfscMisc::distSmry(pct.var, method = 'venter')

age.est.rf |> 
  lapply(function(rf) enframe(rf$predicted, name = 'swfsc.id', value = 'age.est')) |> 
  bind_rows() |> 
  filter(age.est <= 80 & age.est >= 0) |> 
  group_by(swfsc.id) |> 
  summarize(
    mode.age.est = modeest::venter(age.est),
    lower = unname(HDInterval::hdi(age.est)['lower']),
    upper = unname(HDInterval::hdi(age.est)['upper']),
    .groups = 'drop'
  ) |> 
  left_join(age.df, by = 'swfsc.id') |> 
  mutate(age.confidence = factor(age.confidence)) |> 
  ggplot() + 
  geom_abline(intercept = 0, slope = 1) +
  geom_segment(aes(x = age.min, xend = age.max, y = mode.age.est, yend = mode.age.est, color = age.confidence), alpha = 0.4) +
  geom_segment(aes(x = age.best, xend = age.best, y = lower, yend = upper, color = age.confidence), alpha = 0.4) +
  geom_point(aes(age.best, mode.age.est, fill = age.confidence), color = 'white', size = 4, shape = 21) +
  scale_color_manual(values = conf.colors) +
  scale_fill_manual(values = conf.colors) +
  labs(x = 'CRC Age', y = 'Estimated Age') +
  theme(legend.position = 'top')
