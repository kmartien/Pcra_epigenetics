rm(list = ls())
library(tidyverse)
load("R/Eric_code/age_and_methylation_data.rdata")

sites.2.use <- "RFsites" #"All" or "RFsites"

pred <- bind_rows(
  readRDS(paste0('R/Eric_code/random forest/rf_best_', sites.2.use, '.rds')) |> 
    mutate(type = 'best'),
  readRDS(paste0('R/Eric_code/random forest/rf_ran_age_', sites.2.use, '.rds')) |> 
    mutate(type = 'ran.age'),
  readRDS(paste0('R/Eric_code/random forest/rf_ran_age_meth_', sites.2.use, '.rds')) |> 
    mutate(type = 'ran.age.meth')
)
  
errSmry <- function(df) {
  df |> 
    summarize(
      mean.dev = mean(dev),
      median.dev = median(dev),
      mode.dev = modeest::venter(dev),
      mean.se = mean(sq.err),
      median.se = median(sq.err),
      mode.se = modeest::venter(sq.err), 
      cor.age = cor(age.pred, age.resp),
      .groups = 'drop'
    )
}

smry <- pred |> 
  left_join(
    select(age.df, c(swfsc.id, age.confidence)),
    by = 'swfsc.id'
  ) |> 
  mutate(age.confidence = as.character(age.confidence)) |> 
  group_by(type, age.confidence) |> 
  errSmry() |> 
  bind_rows(
    pred |> 
      group_by(type) |> 
      errSmry() |> 
      mutate(age.confidence = 'All')
  ) |> 
  arrange(type, age.confidence)
print(smry)


# Squared error distribution ----------------------------------------------

pred |> 
  left_join(age.df, by = 'swfsc.id') |> 
  mutate(age.confidence = factor(age.confidence)) |> 
  ggplot() +
  geom_histogram(aes(sq.err, fill = age.confidence)) +
  scale_fill_manual(values = conf.colors) +
  scale_x_log10() +
  facet_grid(type ~ age.confidence, scales = 'free') + 
  theme(legend.position = 'top')


# Predicted vs best age ---------------------------------------------------

pred |> 
  group_by(type, swfsc.id) |> 
  summarize(
    age.pred = modeest::venter(age.pred)#,
#    lci = modeest::venter(lci),
#    uci = modeest::venter(uci)
  ) |> 
  left_join(age.df, by = 'swfsc.id') |> 
  mutate(age.confidence = factor(age.confidence)) |> 
  ggplot() + 
  geom_abline(intercept = 0, slope = 1) +  
#  geom_segment(aes(x = age.best, xend = age.best, y = lci, yend = uci, color = age.confidence), alpha = 0.5) +
  geom_segment(aes(x = age.min, xend = age.max, y = age.pred, yend = age.pred, color = age.confidence), alpha = 0.5) +
  geom_point(aes(x = age.best, y = age.pred, color = age.confidence), size = 3) +
  scale_color_manual(values = conf.colors) +
  labs(x = 'CRC age', y = 'GAM predicted age') +
  facet_grid(age.confidence ~ type) +
  theme(legend.position = 'none')


# Residual vs best age ----------------------------------------------------

pred |> 
  group_by(type, swfsc.id) |> 
  summarize(
    age.pred = modeest::venter(age.pred),
#    lci = modeest::venter(lci),
#    uci = modeest::venter(uci),
#    resid.lci = unname(quantile(resid, 0.025)),
#    resid.uci = unname(quantile(resid, 0.975)),
    resid = modeest::venter(resid)
  ) |> 
  left_join(age.df, by = 'swfsc.id') |>
  mutate(age.confidence = factor(age.confidence)) |> 
  ggplot() + 
  geom_hline(yintercept = 0) + 
#  geom_segment(aes(x = age.best, xend = age.best, y = resid.lci, yend = resid.uci, color = age.confidence), alpha = 0.5) +
  geom_segment(aes(x = age.min, xend = age.max, y = resid, yend = resid, color = age.confidence), alpha = 0.5) +
  geom_point(aes(x = age.best, y = resid, color = age.confidence), size = 3) +
  scale_color_manual(values = conf.colors) +
  labs(x = 'CRC age', y = 'Residual') +
  facet_grid(age.confidence ~ type, scales = 'free_y') +
  theme(legend.position = 'none')
