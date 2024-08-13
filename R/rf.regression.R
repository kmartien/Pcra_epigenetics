library(dplyr)
library(tidyverse)
library(randomForest)
library(rfPermute)
library(ggplot2)
source("R/functions/plotting.funcs.R")
source("R/data.prep/combine.age.and.meth.data.R")
source("R/select.RF.important.sites.R")
source("R/select.glmnet.chosen.sites.R")
#load("data/training_set_params.rda")
load("data/model.params.rda")
load("data/color.palettes.rda")

min.cov <- 100
meth.type <- "logit" # choose "pct", "pct.no.zero", or "logit"
select.sites <- "glmnet" # or "rf"
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
  
  all.samples <- mutate(all.samples, age.range = age.max - age.min, .after = age.max) %>% 
    mutate(age.midpoint = age.range/2 + age.min, .after = age.range) %>%
    mutate(wt.range = (1 - .2*floor(age.range/10)), .after = wt)
  
  data.complete <- filter(all.samples, age.confidence >= training.min.CR) %>%
    select(c(id, age.best, wt, all_of(dat$site.names)))
  NAs.by.sample <- sapply(1:nrow(data.complete), function(s){
    length(which(is.na(data.complete[s,4:ncol(data.complete)])))
  })
  if(sum(NAs.by.sample) > 0) data.complete <- data.complete[-which(NAs.by.sample>0),]
  
  rf.data.age <- select(data.complete, -c(id, wt))
  wt <- data.complete$wt
  
  date()
  rf.age <- rfPermute(
    age.best ~ .,
    rf.data.age,
    ntree = 2000,
    num.rep = 2,
    weights = wt
  )
  date()
  data.complete$predicted.age <- rf.age$rf$predicted
  rf.res <- select(data.complete, c(id, age.best, wt, predicted.age)) %>% 
    left_join(select(all.samples, c(id, sex, age.confidence))) %>% 
    mutate(error = predicted.age - age.best)
  
  mae <- do.call(rbind, lapply(training.min.CR:5, function(cr){
    x <- filter(rf.res, age.confidence == cr)
    age.error <- mutate(x, age.error = abs(error)) %>%
      select(age.error) 
    cor.coeff <- round(cor.test(x$age.best, x$predicted.age, method = "pearson")$estimate,2)
    return(c(CR = cr, MAE = median(age.error$age.error), corr = cor.coeff))
  }))
  
  write.csv(mae, file =paste0("results-raw/rf.mae-", description, ".csv"))
  
  save(rf.age, rf.res, mae, file = paste0("results/rf.regression-", description, ".rda"))
  
  # ImpData <- as.data.frame(importance(rf.age$rf)) %>% filter(`%IncMSE` > 0)
  # ImpData$Var.Names <- row.names(ImpData)
  # 
  # ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  #   geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  #   geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  #   theme_light() +
  #   coord_flip() +
  #   theme(
  #     legend.position="bottom",
  #     panel.grid.major.y = element_blank(),
  #     panel.border = element_blank(),
  #     axis.ticks.y = element_blank()
  #   )
  # 
  # error.hist <- loov.hist(rf.res, min.cr = 2)
  # jpeg(filename = paste0("results-raw/rf.regression.error.histogram", description, ".jpg"))
  # error.hist
  # dev.off()
  # 
  # jpeg(filename = paste0("results-raw/rf.regression.plot.", description, ".jpg"))
  # plot.loov.res(rf.res, min.CR = 2)
  # dev.off()
  # 
  #imp <- rownames_to_column(data.frame(rf.age$rf$importance), var = "site") %>% select(-IncNodePurity)
  #jpeg(filename = paste0("results-raw/rf.regression.site.importance.", description, ".jpg"))
  #plot.site.importance(imp)
  #dev.off()
})

