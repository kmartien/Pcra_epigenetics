library(tidyverse)
library(mgcv)
library(e1071)
#load("../age_and_methylation_data.rdata")
#source('../misc_funcs.R')
load("R/Eric_code/age_and_methylation_data.rdata")
source('R/Eric_code/misc_funcs.R')

sites.2.use <- "All" #"All" or "RFsites"

# predict testing data from fitted model
predictTest <- function(fit, cv.df, sites, resp){
  if(is.null(fit)) return(NULL)
  
  pred <- predict(
    fit, 
    select(cv.df, all_of(sites))
  )

  tibble(
    swfsc.id = cv.df$swfsc.id,
    age.resp = cv.df[[resp]],
    age.pred = unname(ifelse(pred < 0, 0, pred)),
  ) |> 
    mutate(
      sq.err.resp = (age.resp - age.pred) ^ 2,
      dev.resp = abs(age.resp - age.pred),
      sq.err.best = (cv.df$age.best - age.pred) ^ 2,
      dev.best = abs(cv.df$age.best - age.pred)
    )
}

sites <- sites.to.keep
if(sites.2.use == "RFsites"){
  # select important sites from Random Forest
  sites <- readRDS('R/Eric_code/random forest/rf_site_importance.rds') |> 
    filter(pval <= 0.10) |> 
    pull('loc.site')
}

# combine age data with mean logit(Pr.methylation)
model.df <- age.df |>
  filter(swfsc.id %in% ids.to.keep) |>
  left_join(
    logit.meth.normal.params |> 
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  ) 

# CR 4 - 5 samples for model fitting
cr45.ids <- model.df |> 
  filter(age.confidence %in% 4:5) |> 
  pull('swfsc.id')

# number of random age sample replicates
nrep <- 1000

# Cross-validation model for CR 4 & 5 -------------------------------------

tune.obj <- tune(svm, 
                  age.best ~ .,
                  data = select(filter(model.df, swfsc.id %in% cr45.ids), c(age.best, all_of(sites))),
                  ranges = list(
                    cost = 10^(seq(-4, 5, 0.1)),
                    gamma = 10^(seq(-5, 4, 0.1))),
                  tunecontrol = tune.control(sampling = "cross"),
                  cross = 2)
save(tune.obj, file = paste0("R/Eric_code/svm/tune.obj.", sites.2.use, ".rda"))
load(paste0("R/Eric_code/svm/tune.obj.", sites.2.use, ".rda"))

cv.pred.4_5 <- lapply(cr45.ids, function(cv.id) {
  cat(cv.id, '\n', sep = '')
    
  # fit to best age and predict test sample age
  fit <- svm(
    age.best ~ ., 
    data = select(filter(model.df, swfsc.id %in% setdiff(cr45.ids, cv.id)), c(age.best, all_of(sites))),
    cost = tune.obj$best.parameters$cost,
    gamma = tune.obj$best.parameters$gamma)
  
  pred.best <- predictTest(fit, filter(model.df, swfsc.id == cv.id), sites, 'age.best')

  # nrep draws of age from CWSN
   #pred.ran <- parallel::mclapply(1:nrep, function(j) {
   pred.ran <- lapply(1:nrep, function(j) {
      #   # randomly sample ages, fit model and predict test sample age
     model.df$age.ran <- sapply(1:nrow(model.df), ranAge, df = model.df)
     fit <- svm(
       age.ran ~ ., 
       data = select(filter(model.df, swfsc.id %in% setdiff(cr45.ids, cv.id)), c(age.ran, all_of(sites))),
       cost = tune.obj$best.parameters$cost,
       gamma = tune.obj$best.parameters$gamma)
     predictTest(fit, filter(model.df, swfsc.id == cv.id), sites, 'age.ran')
   #}, mc.cores = 8) |> 
   }) |> 
  bind_rows()
   
   list(pred.best = pred.best, pred.ran = pred.ran)
}) |> 
  list_transpose()
cv.pred.4_5$pred.ran <- bind_rows(cv.pred.4_5$pred.ran)


# Full CR 4 & 5 model to predict CR 2 & 3 ---------------------------------

# fit to best age and predict CR 2 and 3
full.fit <- fitTrain(filter(model.df, age.confidence %in% 4:5), sites, 'age.best')
full.pred.2_3 <- list(
  pred.best = predictTest(full.fit, filter(model.df, age.confidence %in% 2:3), 'age.best')
)

# nrep draws of age from CWSN
full.pred.2_3$pred.ran <- parallel::mclapply(1:nrep, function(j) {
  # randomly sample ages, fit model and predict CR 2 and 3 ages
  model.df$age.ran <- sapply(1:nrow(model.df), ranAge, df = model.df)
  fit <- fitTrain(filter(model.df, age.confidence %in% 4:5), sites, 'age.ran')
  predictTest(fit, filter(model.df, age.confidence %in% 2:3), 'age.ran')
}, mc.cores = 8) |> 
  bind_rows()


# Combine and save --------------------------------------------------------

#pred.best <- bind_rows(cv.pred.4_5$pred.best, full.pred.2_3$pred.best)
#pred.ran <- bind_rows(cv.pred.4_5$pred.ran, full.pred.2_3$pred.ran)
pred.best <- cv.pred.4_5$pred.best
pred.ran <- cv.pred.4_5$pred.ran
save(sites, nrep, pred.best, pred.ran, tune.obj, file = paste0("R/Eric_code/svm/svm_", sites.2.use, "_replicate_model_results.rdata"))
