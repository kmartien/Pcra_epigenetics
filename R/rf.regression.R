library(dplyr)
library(tidyverse)
library(randomForest)
library(rfPermute)
library(ggplot2)
source("R/functions/plotting.funcs.R")
source("R/data.prep/combine.age.and.meth.data.R")
load("data/color.palettes.rda")

min.cov <- 100
meth.type <- "logit" # choose "pct", "pct.no.zero", or "logit"
weighted <- TRUE
description <- paste0(meth.type, ".min.cov.",min.cov)
training.min.CR <- 2

dat <- combine.age.and.meth.data(description)
all.samples <- dat$all.samples
all.samples <- mutate(all.samples, age.range = age.max - age.min, .after = age.max) %>% 
  mutate(age.midpoint = age.range/2 + age.min, .after = age.range) %>%
  mutate(wt.range = (1 - .2*floor(age.range/10)), .after = wt)
if (weighted) description <- paste0("weighted.", description)
if(!weighted) all.samples$wt <- 1
description <- paste0("age.midpoint.", description)

data.complete <- select(all.samples, c(id, age.best, age.midpoint, wt.range, all_of(dat$site.names)))
NAs.by.sample <- sapply(1:nrow(data.complete), function(s){
  length(which(is.na(data.complete[s,5:ncol(data.complete)])))
})
data.complete <- data.complete[-which(NAs.by.sample>0),]

rf.data.age <- select(data.complete, -c(id, age.midpoint, wt.range))

rf.age <- rfPermute(
  age.best ~ .,
  rf.data.age,
  ntree = 2500,
  num.rep = 10
)

data.complete$predicted.age <- rf.age$rf$predicted
rf.res <- select(data.complete, c(id, age.best, wt, predicted.age)) %>% 
  left_join(select(all.samples, c(id, age.confidence))) %>% 
  mutate(error = predicted.age - age.best)

mae <- do.call(rbind, lapply(2:5, function(cr){
  age.error <- filter(rf.res, age.confidence == cr) %>% mutate(age.error = abs(error)) %>%
    select(age.error) 
  return(c(CR = cr, MAE = median(age.error$age.error)))
}))

save(rf.age, rf.res, mae, file = paste0("results/rf.regression.", description, ".rda"))

ImpData <- as.data.frame(importance(rf.age$rf)) %>% filter(`%IncMSE` > 0)
ImpData$Var.Names <- row.names(ImpData)

ggplot(ImpData, aes(x=Var.Names, y=`%IncMSE`)) +
  geom_segment( aes(x=Var.Names, xend=Var.Names, y=0, yend=`%IncMSE`), color="skyblue") +
  geom_point(aes(size = IncNodePurity), color="blue", alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme(
    legend.position="bottom",
    panel.grid.major.y = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank()
  )

error.hist <- loov.hist(rf.res, min.cr = 2)
jpeg(filename = paste0("results-raw/rf.regression.error.histogram", description, ".jpg"))
error.hist
dev.off()

jpeg(filename = paste0("results-raw/rf.regression.plot.", description, ".jpg"))
plot.loov.res(rf.res, min.CR = 2)
dev.off()

imp <- rownames_to_column(data.frame(rf.age$rf$importance), var = "site") %>% select(-IncNodePurity)
jpeg(filename = paste0("results-raw/rf.regression.site.importance.", description, ".jpg"))
plot.site.importance(imp)
dev.off()


