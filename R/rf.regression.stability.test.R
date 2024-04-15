library(dplyr)
library(tidyverse)
library(randomForest)
library(rfPermute)
library(ggplot2)
source("R/functions/plotting.funcs.R")
source("R/data.prep/combine.age.and.meth.data.R")
load("data/color.palettes.rda")

nreps <- 10
min.cov <- 100
meth.type <- "logit" # choose "pct", "pct.no.zero", or "logit"
weighted <- FALSE
description <- paste0(meth.type, ".min.cov.",min.cov)
training.min.CR <- 2

dat <- combine.age.and.meth.data(description)
all.samples <- dat$all.samples
if (weighted) description <- paste0("weighted.", description)
if(!weighted) all.samples$wt <- 1

data.complete <- select(all.samples, c(id, age.best, wt, all_of(dat$site.names)))
NAs.by.sample <- sapply(1:nrow(data.complete), function(s){
  length(which(is.na(data.complete[s,5:ncol(data.complete)])))
})
data.complete <- data.complete[-which(NAs.by.sample>0),]

rf.data.age <- select(data.complete, -c(id, wt))

# Run RandomForest nrep times and save raw results
rf.list <- lapply(1:nreps, function(i){
  print(i)
  print(date())
  rf.age <- rfPermute(
    age.best ~ .,
    rf.data.age,
    ntree = 2500,
    num.rep = 10,
    num.cores = 4
  )
})
save(rf.list, file = "results/rf.regression.stability.results.rda")

###################################################################
# Summarize age prediction errors
age.error.sum <- do.call(rbind, lapply(1:nreps, function(i){
  age.error <- data.frame(data.complete$id, rf.list[[i]]$predicted - rf.data.age$age.best) 
  cbind(age.error, i)
}))
names(age.error.sum) <- c("id", "age.error", "rep")
age.error.sum <- right_join(select(dat$all.samples, c(id,sex,age.best,age.confidence)), age.error.sum) %>%
  arrange(age.best)

# boxplot of age errors for each CR = 5 individual across runs; individuls
# ordered along x-axis from youngest to oldest
cr5.inds <- filter(age.error.sum, age.confidence == 5) %>% select(c(id, age.best)) %>% unique()
age.error.boxplot <- ggplot(age.error.sum) +
  geom_boxplot(aes(x = id, y = age.error)) +
  facet_wrap(~age.confidence, scales = "free_x", ncol = 2)
  #  scale_x_discrete(limits = cr5.inds$id)
jpeg(file = "results-raw/rf.age.error.boxplot.jpg", width = 1500, height = 1000)
age.error.boxplot
dev.off()

# calculate median age error for CR = 5 animals from each run
mae <- filter(age.error.sum, age.confidence == 5) %>% 
  group_by(rep) %>% summarise(mae = median(abs(age.error)))

# histogram of mae across replicate runs
mae.hist <- ggplot(mae) + 
  geom_histogram(aes(x = mae), col = "grey25") +
#  scale_x_continuous(breaks = seq(0,18,3)) +
  ggtitle("MAE across replicate rf regressions")


##################################################################
# Summarize site importance
importance.sum <- do.call(rbind, lapply(1:nreps, function(i){
  pctIncMSE <- rf.list[[i]]$importance[,1] %>% 
    data.frame() %>% rownames_to_column(var = "site")
  cbind(pctIncMSE, i)
}))
names(importance.sum) <- c("site", "pctIncMSE", "rep")

importance.sum <- bind_cols(do.call(rbind, strsplit(importance.sum$site, split = "[.]")), importance.sum)
names(importance.sum) <- c("locus", "pos", "site", "pctIncMSE", "rep")

# boxplot of importance for each site across runs
importance.boxplot <- ggplot(importance.sum) +
  geom_boxplot(aes(x = pos, y = pctIncMSE)) +
  facet_wrap(~ locus, scales = "free")
jpeg(file = "results-raw/rf.importance.boxplot.jpg", width = 900, height = 1200)
importance.boxplot
dev.off()

