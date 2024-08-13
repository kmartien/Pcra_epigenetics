library(dplyr)
library(tidyverse)
library(ggplot2)
library(gridExtra)
source("R/functions/plotting.funcs.R")
load("data/model.params.rda")
load("data/color.palettes.rda")

age.transform <- "ln"
model <- 30
params <- model.params[model,]
description <- paste0("model", params$model, "-logit_mincov100")

dat <- list()
load(paste0("results/glmnet.LOOV.plainjane-", age.transform, ".results.rda"))
dat$ENR <- loov.res
load(paste0("results/glmnet.LOOV.plainjane.alphaHalf-", age.transform, ".results.rda"))
dat$ENR0.5 <- loov.res
load(paste0("results/rf.regression-plainjane-", age.transform, ".rda"))
dat$RF <- rf.res
load(paste0("results/svm.regression-plainjane-", age.transform, ".rda"))
dat$SVM <- loov.age
dat$SVM$age.best <- as.integer(dat$SVM$age.best)

plots <- lapply(1:length(dat), function(i){
  p <- plot.loov.res(dat[[i]], min.CR = 4)
  p$p.loov$labels$title <- paste0(names(dat)[i])
  return(p$p.loov)
})

plots[[1]]$labels$title <- "ENR optimized alpha = 0.1"
plots[[2]]$labels$title <- "ENR alpha = 0.5"
plots$nrow <- 2
jpeg(file = paste0("results-raw/plainjane.regression.plots-",age.transform,".jpg"), width = 1000, height = 800)
do.call(grid.arrange, plots[-2])
dev.off()

# box plots
age.errors.long <- do.call(bind_rows, lapply(1:length(dat), function(i){
  bind_cols(method = names(dat)[i], dat[[i]])
}))
age.errors.long$abs.error <- abs(age.errors.long$error)
age.errors.long$CR <- as.factor(age.errors.long$age.confidence)
box.plot <- 
  ggplot(age.errors.long) +
  geom_boxplot(aes(x = CR, y = abs.error)) +
  labs(x = "Confidence Rating",
       y = "Absolute age error (yrs)") +  
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 30, hjust=1)) +
  #  facet_wrap(~training.cr, nrow = 3, labeller = labeller(training.cr = training.cr.labs))  
  facet_wrap(~method, nrow = 2)  

# MAE for CR=5 males vs. females
MAE.by.sex <- do.call(rbind, lapply(1:length(dat), function(i){
  filter(dat[[i]], age.confidence == 5) %>% group_by(sex) %>% 
    summarise(MAE = median(abs(error))) %>% bind_cols(method = names(dat)[i])
}))

# MAE by age and CR
breaks <- c(0,10,25,40)
MAE.by.age <- do.call(cbind, lapply(1:length(dat), function(i){
  do.call(rbind, lapply(1:(length(breaks)-1), function(a){
    filter(dat[[i]], age.best >= breaks[a]) %>% filter(age.best < breaks[a+1]) %>% 
      group_by(age.confidence) %>% summarise(MAE = median(abs(error))) %>% bind_cols(age_bins = breaks[a])
  }))
}))
write.csv(MAE.by.age, file = paste0("results-raw/MAE.by.age.across.methods-", age.transform, ".csv"))

# compare duplicates
dupe.ids <- calibration.set.complete$crc.id[which(duplicated(calibration.set.complete$crc.id))]
dupe.sum <- do.call(bind_rows, lapply(c(1,3,4), function(m){
  dupes <- select(calibration.set.complete, c(crc.id, id, date.biopsy, age.confidence)) %>%
    left_join(dat[[m]]) %>% filter(crc.id %in% dupe.ids)
  do.call(rbind, lapply(unique(dupes$crc.id), function(i){
    inds <- filter(dupes, crc.id == i) %>% arrange(age.best)
    actual.diff <- difftime(as.Date(inds$date.biopsy[2], "%d-%b-%y"), as.Date(inds$date.biopsy[1], "%d-%b-%y"), units = "days")/365
#    actual.diff <- inds$age.best[2] - inds$age.best[1]
    predicted.diff <- inds$predicted.age[2] - inds$predicted.age[1]
    return(data.frame(method = names(dat)[m], crc.id = i, actual.diff = actual.diff, predicted.diff = predicted.diff))
  }))
}))
pair.plot <- ggplot(dupe.sum) +
  geom_point(aes(x = actual.diff, y = predicted.diff), size = 2) +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Actual age difference", y = "Predicted age difference") +
  theme(text = element_text(size = 24)) +
  facet_wrap(~method, nrow = 2)
jpeg(filename = paste0("results-raw/pair.plot-", age.transform, ".jpg"), width = 960, height = 960)
pair.plot
dev.off()

