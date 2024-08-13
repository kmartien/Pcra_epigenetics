library(dplyr)
library(tidyverse)
library(randomForest)
library(rfPermute)
library(ggplot2)
source("R/functions/plotting.funcs.R")
source("R/data.prep/combine.age.and.meth.data.R")
load("data/color.palettes.rda")
load("data/training_set_params.rda")

nreps <- 100 # reps in stability test
num.reps <- 100 # reps used by rfPermute to calculate p-values
min.cov <- 100
meth.type <- "logit" # "logit" or "pct"

training.sets <- filter(training_set_params, sites2use == "all")
lapply(1:nrow(training.sets), function(t){
  description <- paste0("train", training.sets$id[t], "_", meth.type, "_mincov",min.cov)
  weight.type <- training.sets$weight[t]
  training.min.CR <- training.sets$minCR[t]

  dat <- combine.age.and.meth.data(paste0(meth.type, "_mincov",min.cov))
  all.samples <- dat$all.samples
  if (weight.type == "linear") all.samples$wt <- all.samples$age.confidence/5
  if (weight.type == "none") all.samples$wt <- 1
  

  data.complete <- select(all.samples, c(id, age.best, wt, all_of(dat$site.names)))
  NAs.by.sample <- sapply(1:nrow(data.complete), function(s){
    length(which(is.na(data.complete[s,5:ncol(data.complete)])))
  })
  data.complete <- data.complete[-which(NAs.by.sample>0),]
  
  rf.data.age <- select(data.complete, -c(id,wt))
  wt <- data.complete$wt
  
  # Run RandomForest nrep times and save raw results
  rf.list <- lapply(1:nreps, function(i){
    print(i)
    print(date())
    rf.age <- rfPermute(
      age.best ~ .,
      rf.data.age,
      weights = wt,
      ntree = 2500,
      num.rep = num.reps,
      num.cores = 4
    )
  })
  
  # Summarize site importance
  importance.sum <- do.call(rbind, lapply(1:nreps, function(i){
    pctIncMSE <- rf.list[[i]]$rf$importance[,1] %>% 
      data.frame() %>% rownames_to_column(var = "site")
    p.val <- importance(rf.list[[i]]) %>% data.frame() %>% 
      select(X.IncMSE.pval) %>% rownames_to_column(var = "site")
    cbind(full_join(pctIncMSE, p.val), i)
  }))
  names(importance.sum) <- c("site", "pctIncMSE", "pval", "rep")
  
  importance.sum <- bind_cols(do.call(rbind, strsplit(importance.sum$site, split = "[.]")), importance.sum)
  names(importance.sum) <- c("locus", "pos", "site", "pctIncMSE", "pval", "rep")

  description <- paste0("train", training.sets$training.set[t], "_", meth.type, "_mincov",min.cov)
  load(file = paste0("results/rf_stability_results-", description, ".rda"))  
  imp.pvals <- do.call(cbind, lapply(rf.list, function(i){
    i$pval[,1,2]
  })) %>% data.frame()
  imp.pvals$prob <- do.call(rbind, lapply(1:nrow(imp.pvals), function(i){
    length(which(imp.pvals[i,] <= 0.05))
  }))
  
  save(rf.list, importance.sum, imp.pvals, file = paste0("results/rf_stability_results-", description, ".rda"))
  
  # bar plot of frequency site had significant p-value
  site.incl.sum <- select(importance.sum, c(locus, pos, site, pval)) %>%
    filter(pval <= 0.05) %>% group_by(site) %>% summarise(prob = length(pval))
  g.site.incl.bar <- 
    ggplot(site.incl.sum, aes(x = site, y = (prob/nreps))) +
    geom_bar(stat = "identity") +
    scale_x_discrete(guide = guide_axis(angle = 90))
  jpeg(file = paste0("results-raw/rf.site.inclusion.", description, ".jpg"), width = 900, height = 500)
  g.site.incl.bar
  dev.off()
  
})

# # boxplot of importance for each site across runs
# importance.boxplot <- ggplot(importance.sum) +
#   geom_boxplot(aes(x = pos, y = pctIncMSE)) +
#   facet_wrap(~ locus, scales = "free_x")
# jpeg(file = paste0("results-raw/rf.importance.boxplot.fixed.scale.",description, ".jpg"), width = 900, height = 1200)
# importance.boxplot
# dev.off()

###################################################################
# # Summarize age prediction errors
# age.error.sum <- do.call(rbind, lapply(1:nreps, function(i){
#   age.error <- data.frame(data.complete$id, rf.list[[i]]$rf$predicted - rf.data.age$age.best) 
#   cbind(age.error, i)
# }))
# names(age.error.sum) <- c("id", "age.error", "rep")
# age.error.sum <- right_join(select(dat$all.samples, c(id,sex,age.best,age.confidence)), age.error.sum) %>%
#   arrange(age.best)
# 
# # boxplot of age errors for each CR = 5 individual across runs; individuls
# # ordered along x-axis from youngest to oldest
# cr5.inds <- filter(age.error.sum, age.confidence == 5) %>% select(c(id, age.best)) %>% unique()
# age.error.boxplot <- ggplot(age.error.sum) +
#   geom_boxplot(aes(x = id, y = age.error)) +
#   facet_wrap(~age.confidence, scales = "free_x", ncol = 2)
#   #  scale_x_discrete(limits = cr5.inds$id)
# jpeg(file = paste0("results-raw/rf.age.error.boxplot.", description, ".jpg"), width = 1500, height = 1000)
# age.error.boxplot
# dev.off()
# 
# # calculate median age error for CR = 5 animals from each run
# mae <- filter(age.error.sum, age.confidence == 5) %>% 
#   group_by(rep) %>% summarise(mae = median(abs(age.error)))
# 
# # histogram of mae across replicate runs
# mae.hist <- ggplot(mae) + 
#   geom_histogram(aes(x = mae), col = "grey25") +
# #  scale_x_continuous(breaks = seq(0,18,3)) +
#   ggtitle("MAE across replicate rf regressions")
# 
# 
##################################################################