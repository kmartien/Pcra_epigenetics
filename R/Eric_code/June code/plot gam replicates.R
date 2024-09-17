rm(list = ls())
library(tidyverse)
library(mgcv)
source('../misc_funcs.R')
load("../age_and_methylation_data.rdata")
load('gam_replicates_20240715_2136.rdata')

ids.to.use <- age.df |> 
  filter(swfsc.id %in% ids.to.keep & age.confidence %in% 4:5) |> 
  pull('swfsc.id')

test.ids <- setdiff(age.df$swfsc.id, ids.to.use)

pred.age.df <- lapply(gam.rep, function(x) {
  parallel::mclapply(x$gam.fit, function(fit) {
    if(is.null(fit)) return(NULL)
    df <- sampleAgeMeth(age.df, logit.meth.normal.params, 1, 1)[[1]] |> 
      filter(swfsc.id %in% c(x$cv.id, test.ids)) |> 
      column_to_rownames('swfsc.id')
    predict(fit, df, type = 'response') |> 
      enframe(name = 'swfsc.id', value = 'pred.age')
  }, mc.cores = 10) |> 
    bind_rows()
}) |>
  bind_rows()
save(pred.age.df, file = 'gam_pred_age.rdata')


pred.age.df |> 
  filter(pred.age <= 80 & pred.age >= 0) |> 
  group_by(swfsc.id) |> 
  summarize(
    median.pred.age = median(pred.age),
    mode.pred.age = modeest::venter(pred.age),
    lower = unname(HDInterval::hdi(pred.age)['lower']),
    upper = unname(HDInterval::hdi(pred.age)['upper']),
    .groups = 'drop'
  ) |> 
  left_join(age.df, by = 'swfsc.id') |> 
  mutate(age.confidence = factor(age.confidence)) |> 
  ggplot(aes(x = age.best)) + 
  geom_abline(intercept = 0, slope = 1) +
  geom_segment(aes(xend = age.best, y = lower, yend = upper, color = age.confidence), alpha = 0.5) +
  geom_segment(aes(x = age.min, xend = age.max, y = mode.pred.age, yend = mode.pred.age, color = age.confidence), alpha = 0.5) +
  geom_point(aes(y = mode.pred.age, color = age.confidence), size = 3) +
  scale_color_manual(values = conf.colors) +
  labs(x = 'CRC best age', y = 'GAM predicted age') +
  theme(
    legend.position = 'inside',
    legend.position.inside = c(1, 1),
    legend.justification = c(1, 1)
  )

pred.age.df |> 
  filter(pred.age <= 80 & pred.age >= 0) |> 
  group_by(swfsc.id) |> 
  summarize(
    mean.pred.age = mean(pred.age),
    median.pred.age = median(pred.age),
    mode.pred.age = modeest::venter(pred.age),
    .groups = 'drop'
  ) |> 
  left_join(age.df, by = 'swfsc.id') |> 
  mutate(
    mean.sq.err = (mean.pred.age - age.best) ^ 2,
    median.sq.err = (median.pred.age - age.best) ^ 2,
    mode.sq.err = (mode.pred.age - age.best) ^ 2,
    age.confidence = as.character(age.confidence)
  ) |> 
  ggplot() + 
  geom_histogram(aes(mean.sq.err, fill = age.confidence)) +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 1000), labels = c('0.01', '0.1', '1', '10', '1000')) +
  scale_fill_manual(values = conf.colors) +
  labs(x = "MSE") +
  facet_wrap(~ age.confidence, ncol = 1) +
  theme(legend.position = 'top')

pred.age.df |> 
  filter(pred.age <= 80 & pred.age >= 0) |> 
  group_by(swfsc.id) |> 
  summarize(
    mean.pred.age = mean(pred.age),
    median.pred.age = median(pred.age),
    mode.pred.age = modeest::venter(pred.age),
    .groups = 'drop'
  ) |> 
  left_join(age.df, by = 'swfsc.id') |> 
  mutate(
    mean.sq.err = (mean.pred.age - age.best) ^ 2,
    median.sq.err = (median.pred.age - age.best) ^ 2,
    mode.sq.err = (mode.pred.age - age.best) ^ 2,
    age.confidence = as.character(age.confidence)
  ) |> 
  ggplot() +
  geom_point(aes(x = age.best, y = mean.sq.err, color = age.confidence), size = 3) +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10, 1000), labels = c('0.01', '0.1', '1', '10', '1000')) +
  scale_color_manual(values = conf.colors) +
  labs(x = 'CRC best age', y = 'MSE') +
  facet_wrap(~ age.confidence) +
  theme(legend.position = 'top')


pred.age.df |> 
  filter(pred.age <= 80 & pred.age >= 0) |> 
  group_by(swfsc.id) |> 
  summarize(
    mean.pred.age = mean(pred.age),
    median.pred.age = median(pred.age),
    mode.pred.age = modeest::venter(pred.age),
    .groups = 'drop'
  ) |> 
  left_join(age.df, by = 'swfsc.id') |> 
  mutate(
    mean.sq.err = (mean.pred.age - age.best) ^ 2,
    median.sq.err = (median.pred.age - age.best) ^ 2,
    mode.sq.err = (mode.pred.age - age.best) ^ 2,
  ) |> 
  summarize(
    mean.sq.err = mean(mean.sq.err),
    median.sq.err = mean(median.sq.err),
    mode.sq.err = mean(mode.sq.err)
  ) |> 
  mutate(age.confidence = 'All') |> 
  bind_rows(
    pred.age.df |> 
      filter(pred.age <= 80 & pred.age >= 0) |> 
      group_by(swfsc.id) |> 
      summarize(
        mean.pred.age = mean(pred.age),
        median.pred.age = median(pred.age),
        mode.pred.age = modeest::venter(pred.age),
        .groups = 'drop'
      ) |> 
      left_join(age.df, by = 'swfsc.id') |> 
      mutate(
        mean.sq.err = (mean.pred.age - age.best) ^ 2,
        median.sq.err = (median.pred.age - age.best) ^ 2,
        mode.sq.err = (mode.pred.age - age.best) ^ 2,
        age.confidence = as.character(age.confidence)
      ) |> 
      group_by(age.confidence) |> 
      summarize(
        mean.sq.err = mean(mean.sq.err),
        median.sq.err = mean(median.sq.err),
        mode.sq.err = mean(mode.sq.err),
        .groups = 'drop'
      )
  ) |> 
  select(age.confidence, everything())
