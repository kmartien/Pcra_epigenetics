rm(list = ls())
library(tidyverse)

load("glmnet.absolute.age.minCR4.results.rda")
load("../results/model.1b.conf.4-5.no.sex.diagnostics.231206_0530.rdata")

df <- smry$age.pred |> 
  group_by(swfsc.id, type) |> 
  summarize(posterior.mode = modeest::venter(pred.age), .groups = "drop") |> 
  left_join(
    age.sum |> 
      rename(glmnet.age = predicted.age),
    by = c(swfsc.id = 'id')
  ) |> 
  mutate(
    bayesian.diff = posterior.mode - age.point,
    glmnet.diff = glmnet.age - age.point,
    bayesian.glmnet.diff = bayesian.diff - glmnet.diff,
    age.range = age.max - age.min,
    confidence = factor(confidence)
  ) 

df |> 
  select(swfsc.id, type, age.range, confidence, bayesian.diff, glmnet.diff) |> 
  pivot_longer(cols = c(bayesian.diff, glmnet.diff), names_to = 'model', values_to = 'diff') |> 
  ggplot() +
  geom_histogram(aes(diff, fill = confidence)) +
  geom_vline(xintercept = 0) +
  facet_grid(model ~ type) 

df |> 
  select(swfsc.id, type, age.range, confidence, bayesian.diff, glmnet.diff) |> 
  pivot_longer(cols = c(bayesian.diff, glmnet.diff), names_to = 'model', values_to = 'diff') |> 
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_point(aes(age.range, diff, color = confidence), size = 3) +
  facet_grid(model ~ type) 

df |> 
  select(swfsc.id, type, age.point, confidence, bayesian.diff, glmnet.diff) |> 
  pivot_longer(cols = c(bayesian.diff, glmnet.diff), names_to = 'model', values_to = 'diff') |> 
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_point(aes(age.point, diff, color = confidence), size = 3) +
  facet_grid(model ~ type) 


df |> 
  ggplot() +
  geom_vline(xintercept = 0) +
  geom_histogram(aes(bayesian.glmnet.diff, fill = confidence)) +
  facet_wrap(~ type)

df |> 
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_point(aes(age.range, bayesian.glmnet.diff, color = confidence), size = 3) +
  facet_wrap(~ type)

df |> 
  ggplot() +
  geom_hline(yintercept = 0) +
  geom_point(aes(age.point, bayesian.glmnet.diff, color = confidence), size = 3) +
  facet_wrap(~ type)

