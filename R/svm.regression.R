library(dplyr)
library(tidyverse)
library(e1071)
source("R/data.prep/combine.age.and.meth.data.R")
source("R/functions/LOOV.funcs.R")
source("R/select.RF.important.sites.R")
source("R/select.glmnet.chosen.sites.R")
load("data/model.params.rda")

min.cov <- 100
meth.type <- "logit" # choose "pct", "pct.no.zero", or "logit"
select.sites <- "rf" # or "rf"
description <- paste0(meth.type, "_mincov",min.cov)
dat <- combine.age.and.meth.data(description)

lapply(1:nrow(model.params), function(i){
  print(i)
  params <- model.params[i,]
  description <- paste0("model", params$model, "-", meth.type, "_mincov",min.cov)
  weight.type <- params$weight
  training.min.CR <- params$training.cr
  
  if(params$site.select.regr.meth != "none"){
    if (params$site.select.regr.meth == "glmnet") {
      chosen.sites <- select.glmnet.chosen.sites(params = params, incl.prob.threshold = 50)
    } else chosen.sites <- select.RF.important.sites(params = params, incl.prob.threshold = 50)
    sites2remove <- dat$site.names[-which(dat$site.names %in% chosen.sites$site)]
    dat$all.samples <- select(dat$all.samples, -all_of(sites2remove))
    dat$site.names <- as.character(chosen.sites$site)
  }
  site.names <- dat$site.names
  all.samples <- dat$all.samples
  if (weight.type == "linear") all.samples$wt <- all.samples$age.confidence/5
  if (weight.type == "none") all.samples$wt <- 1
  
  all.samples <- dat$all.samples %>% filter(age.confidence >= training.min.CR)
  all.samples <- mutate(all.samples, age.range = age.max - age.min, .after = age.max) %>% 
    mutate(age.midpoint = age.range/2 + age.min, .after = age.range) %>%
    mutate(wt.range = (1 - .2*floor(age.range/10)), .after = wt) %>%
    rename(cluster.louvain = social.cluster)
  
  data.complete <- select(all.samples, c(id, age.best, numeric.sex, cluster.louvain, wt, all_of(dat$site.names)))
  NAs.by.sample <- sapply(1:nrow(data.complete), function(s){
    length(which(is.na(data.complete[s,6:ncol(data.complete)])))
  })
  if(sum(NAs.by.sample) > 0) data.complete <- data.complete[-which(NAs.by.sample>0),]
  
  # age prediction
  data.age <- select(data.complete, -c(numeric.sex, cluster.louvain, wt))
  
  date()
  tune.obj <- tune(svm, 
                   age.best ~ .,
                   data = select(data.age, -id),
                   ranges = list(
                     cost = 10^(seq(-4, 5, 0.1)),
                     gamma = 10^(seq(-5, 4, 0.1))),
                   tunecontrol = tune.control(sampling = "cross"),
                   cross = 10)
  date()
  
  loov.age <- svm.loov(
    dat = data.age,
    cost = tune.obj$best.parameters$cost,
    gamma = tune.obj$best.parameters$gamma
  ) 
  loov.age <- left_join(loov.age, select(all.samples, c(id,sex,age.confidence, wt)), by = "id")
  
  loov.error.sum <- c(mean = mean(abs(loov.age$error)), median = median(abs(loov.age$error)))
  
  mae <- do.call(rbind, lapply(training.min.CR:5, function(cr){
    x <- filter(loov.age, age.confidence == cr)
    age.error <- mutate(x, age.error = abs(error)) %>%
      select(age.error) 
    cor.coeff <- round(cor.test(x$age.best, x$predicted.age, method = "pearson")$estimate,2)
    return(c(CR = cr, MAE = median(age.error$age.error), corr = cor.coeff))
  }))
  write.csv(mae, file =paste0("results-raw/svm.mae-", description, ".csv"))
  
  save(tune.obj, loov.age, mae, file = paste0("results/svm.regression-", description, ".rda"))

})
