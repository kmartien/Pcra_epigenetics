library(dplyr)
library(tidyverse)
library(e1071)
source("R/data.prep/combine.age.and.meth.data.R")

min.cov <- 100
meth.type <- "logit" # choose "pct", "pct.no.zero", or "logit"
description <- paste0(meth.type, ".min.cov.",min.cov)

dat <- combine.age.and.meth.data(description)
all.samples <- dat$all.samples
all.samples <- rename(all.samples, cluster.louvain = social.cluster)
all.samples$wt <- 1

data.complete <- select(all.samples, c(id, numeric.sex, cluster.louvain, wt, all_of(dat$site.names)))
meth.dat <- select(data.complete, all_of(dat$site.names))
NAs.by.sample <- sapply(1:nrow(meth.dat), function(s){
  length(which(is.na(meth.dat[s,])))
})
if(sum(NAs.by.sample) > 0) data.complete <- data.complete[-which(NAs.by.sample>0),]

# sex classification
data.sex <- select(data.complete, -c(cluster.louvain, wt))
data.sex$numeric.sex <- as.factor(data.sex$numeric.sex)

#Nakamura use cost = 10^(seq(-4,5,0.1)) and gamma = 10^(seq(-5,4,0.1)) for tuning
date()
sex.tune.obj <- tune(svm, 
                 numeric.sex ~ .,
                 data = select(data.sex, -id),
                 ranges = list(
                   cost = 10^(seq(-4, 5, 0.1)),
                   gamma = 10^(seq(-5, 4, 0.1))),
                 tunecontrol = tune.control(sampling = "cross"),
                 cross = 10)
date()

loov.sex <- do.call(rbind,lapply(data.sex$id, function(i){
  
  dat.loov <- filter(data.sex, id != i) %>% select(-id)
  svm.sex <- svm(numeric.sex ~ ., data = dat.loov, 
    cost = sex.tune.obj$best.parameters$cost,
    gamma = sex.tune.obj$best.parameters$gamma)
  oob.meth <- as.matrix(filter(data.sex, id == i) %>%
                          select(-c(id, numeric.sex)))
  oob.sex <- data.sex$numeric.sex[which(data.sex$id == i)]
  predicted.sex <- unname(predict(svm.sex, oob.meth))
  return(c("sex" = oob.sex, "predicted.sex" = predicted.sex))
})) %>% data.frame()

sex.confusion.mat <- table(pred = loov.sex$predicted.sex, true = loov.sex$sex)
write.csv(sex.confusion.mat, file = paste0("results-raw/svm.sex.confusion.mat.", description, ".csv"))

save(data.complete, sex.tune.obj, loov.sex, file = paste0("results/svm.sex.clust.",description, ".rda"))

# social cluster classification
data.clust <- select(data.complete, -c(numeric.sex, wt))
data.clust$cluster.louvain <- as.factor(data.clust$cluster.louvain)

date()
clust.tune.obj <- tune(svm, 
                 cluster.louvain ~ .,
                 data = select(data.clust, -id),
                 ranges = list(
                   cost = 10^(seq(-4, 5, 0.1)),
                   gamma = 10^(seq(-5, 4, 0.1))),
                 tunecontrol = tune.control(sampling = "cross"),
                 cross = 10)
date()

loov.clust <- do.call(rbind,lapply(data.clust$id, function(i){
  
  dat.loov <- filter(data.clust, id != i) %>% select(-id)
  svm.clust <- svm(cluster.louvain ~ ., data = dat.loov, 
                 cost = clust.tune.obj$best.parameters$cost,
                 gamma = clust.tune.obj$best.parameters$gamma)
  oob.meth <- as.matrix(filter(data.clust, id == i) %>%
                          select(-c(id, cluster.louvain)))
  oob.clust <- data.clust$cluster.louvain[which(data.clust$id == i)]
  predicted.clust <- unname(predict(svm.clust, oob.meth))
  return(c("clust" = oob.clust, "predicted.clust" = predicted.clust))
})) %>% data.frame()

clust.confusion.mat <- table(pred = loov.clust$predicted.clust, true = loov.clust$clust)
write.csv(clust.confusion.mat, file = paste0("results-raw/svm.clust.confusion.mat.", description, ".csv"))

save(data.complete, sex.tune.obj, loov.sex, clust.tune.obj, loov.clust, file = paste0("results/svm.sex.clust.",description, ".rda"))
