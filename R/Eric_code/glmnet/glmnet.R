rm(list = ls())
library(tidyverse)
library(glmnet)
load("R/Eric_code/age_and_methylation_data.rdata")

model.df <- age.df |>
  filter(swfsc.id %in% ids.to.keep) |>
  arrange(age.confidence, age.best, age.min, age.max) |> 
  left_join(
    logit.meth.normal.params |>
      filter(loc.site %in% sites.to.keep) |>
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  ) |> 
  column_to_rownames('swfsc.id') |> 
  mutate(inv.var = 1 / age.var)


# standard inverse variance weighted LOO cv.glmnet run at given alpha
ageCVglmnetLOO <- function(df, sites, alpha) {
  cv.glmnet(
    x = as.matrix(df[, sites]),
    y = df$age.best,
    weights = df$inv.var,
    alpha = alpha,
    foldid = 1:nrow(df)
  ) 
}

# return median cvm at 1se for alpha and nrep replicates
median.cvm.1se <- function(alpha, df, nrep) {
  parallel::mclapply(1:nrep, function(i) {
    cv.fit <- tryCatch({
      ageCVglmnetLOO(df, sites.to.keep, alpha)
    }, error = function(e) NULL)
    if(is.null(cv.fit)) NA else {
      cv.fit$cvm[cv.fit$lambda == cv.fit$lambda.1se]
    } 
  }, mc.cores = 6) |> 
    unlist() |> 
    median(na.rm = TRUE)
}

# find optimal alpha (lowest median cvm at 1se)
optimum.alpha <- optim(
  par = c(alpha = 0.3),
  fn = median.cvm.1se,
  method = "Brent",
  lower = 0.0001,
  upper = 0.5,
  df = model.df,
  nrep = 1000
)
optimum.alpha

# multiple LOO runs of the glmnet model with alpha at lowest cvm
glmnet.fit <- parallel::mclapply(1:1000, function(i) {
  tryCatch({
    ageCVglmnetLOO(model.df, sites.to.keep, optimum.alpha$par)
  }, error = function(e) NULL)
}, mc.cores = 10)

#save.image('glmnet.rdata')



# distribution of lambda at 1se across model replicates
lambda.1se <- sapply(glmnet.fit, function(x) {
  if(is.list(x)) x$lambda.1se else NA
})
hist(lambda.1se)
swfscMisc::distSmry(lambda.1se, method = 'venter')


# distribution of predicted age across model replicates
pred.age <- lapply(glmnet.fit, function(x) {
  predict(
    x, 
    as.matrix(model.df[, sites.to.keep]),
    s = 'lambda.1se'
  ) |> 
    as.data.frame() |> 
    rownames_to_column('swfsc.id') |> 
    rename(pred.age = 'lambda.1se')
}) |> 
  bind_rows() 

pred.age |> 
  left_join(
    age.df |> 
      select(swfsc.id, age.best, age.confidence),
    by = 'swfsc.id'
  ) |> 
  mutate(age.confidence = as.character(age.confidence)) |> 
  group_by(age.confidence) |> 
  summarize(mse = mean((age.best - pred.age) ^ 2), .groups = 'drop') |> 
  bind_rows(
    pred.age |> 
      left_join(age.df, by = 'swfsc.id') |> 
      summarize(mse = mean((age.best - pred.age) ^ 2)) |> 
      mutate(age.confidence = 'All')
  )


pred.age |>   
  group_by(swfsc.id) |> 
  summarize(
    median.pred.age = median(pred.age),
    lci = quantile(pred.age, 0.025),
    uci = quantile(pred.age, 0.975),
    .groups = 'drop'
  ) |> 
  left_join(
    age.df,
    by = 'swfsc.id'
  ) |>
  mutate(age.confidence = factor(age.confidence)) |> 
  ggplot(aes(x = age.best)) + 
  geom_abline(intercept = 0, slope = 1) +  
  geom_segment(aes(x = age.best, xend = age.best, y = lci, yend = uci, color = age.confidence), alpha = 0.5) +
  geom_segment(aes(x = age.min, xend = age.max, y = median.pred.age, yend = median.pred.age, color = age.confidence), alpha = 0.5) +
  geom_point(aes(y = median.pred.age, color = age.confidence), size = 3) +
  scale_color_manual(values = conf.colors) +
  labs(x = 'CRC best age', y = 'GLMNET predicted age') +
  theme(
    legend.position = 'inside',
    legend.position.inside = c(1, 1),
    legend.justification = c(1, 1)
  )
