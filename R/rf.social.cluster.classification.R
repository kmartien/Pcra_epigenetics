library(tidyverse)
library(dplyr)
library(randomForest)
library(rfPermute)
library(glmnet)

samps.2.exclude <- c("z0132662","z0190864")
min.cov <- 100
description <- paste0("min.cov.",min.cov)

load(paste0("data/corrected.pct.meth.",description,".Rdata"))
load("data/age.data.rda")
load("data/sites.and.inds.from.Eric.rda")
age.data <- filter(age.data, swfsc.labid %in% ids.to.keep)
sites.2.keep <- gsub("_00", replacement = ".", sites.to.keep)
sites.2.keep <- gsub("_0", replacement = ".", sites.2.keep)
sites.2.keep <- gsub("_", replacement = ".", sites.2.keep)
corrected.pct.meth <- select(corrected.pct.meth, sites.2.keep)
soc.clust.assignments <- read.csv("/Users/Shared/KKMDocuments/Documents/Github.Repos/Pcra/Pcra.database.data/data-raw/social.cluster.assignments.csv")
names(soc.clust.assignments) <- c("crc.id","cluster.louvain")

site.names <- names(corrected.pct.meth)
corrected.pct.meth <- corrected.pct.meth[-which(rownames(corrected.pct.meth) %in% samps.2.exclude),]
corrected.pct.meth <- cbind(id=rownames(corrected.pct.meth),corrected.pct.meth)
names(age.data)[which(names(age.data) %in% c("swfsc.labid","age.best","age.confidence"))] <- c("id","age.point","confidence")
age.data$numeric.sex <- 0
age.data$numeric.sex[which(age.data$sex=="Female")] <- 1
age.data$numeric.sex[which(age.data$sex=="Male")] <- 0
age.data$decade <- floor(age.data$age.point/10)
age.data$decade[which(age.data$decade>3)] <- 3
first.meth.col <- dim(age.data)[2] + 1
age.data <- left_join(age.data, soc.clust.assignments)

all.samples <- left_join(select(age.data, c(id, sex, cluster.louvain, decade)),corrected.pct.meth)
data.complete <- all.samples#[-which(is.na(all.samples$cluster.louvain)),]
data.complete$cluster.louvain <- factor(data.complete$cluster.louvain)
data.complete$sex <- factor(data.complete$sex)
data.complete$decade <- factor(data.complete$decade)

#alternately remove sites and samples to eliminate all NAs

NAs.by.sample <- sapply(1:nrow(data.complete), function(s){
  length(which(is.na(data.complete[s,5:ncol(data.complete)])))
})

data.complete <- data.complete[-which(NAs.by.sample>0),]

# glmnet 10-fold cross volidation sex prediction
alpha.vals <- seq(from = 0.1, to = 0.9, by = 0.1)

data.sex <- select(data.complete, -c(cluster.louvain, decade))
data.sex$numeric.sex <- 0
data.sex$numeric.sex[which(data.sex$sex=="Female")] <- 1

x.meth <- as.matrix(select(data.sex, matches(site.names)))
y.sex <- as.matrix(data.sex$numeric.sex)

sex.alpha <- sapply(alpha.vals, function(a){
  cvfit <- cv.glmnet(x.meth, y.sex, alpha = a, family = "binomial", type.measure = "class")
  error.rate.at.lambdamin <- cvfit$cvm[cvfit$index[1]]
  return(c(a, error.rate.at.lambdamin))
})
sex.best.alpha <- sex.alpha[1,which(sex.alpha[2,] == min(sex.alpha[2,]))]


# glmnet LOOV Sex
loov.sex.res <- do.call(rbind,lapply(data.sex$id, function(i){
  
  dat.loov <- data.sex[-which(data.sex$id == i),]
  x.meth <- as.matrix(select(dat.loov, matches(site.names)))
  y.sex <- as.matrix(dat.loov$numeric.sex)
  
  cvfit <- cv.glmnet(x.meth,y.sex, alpha = sex.best.alpha, family = "binomial", type.measure = "class")
  #  corr.coef <- coef(cvfit, s = "lambda.min")
  oob.meth <- as.matrix(filter(data.sex, id == i) %>%
                          select(matches(site.names)))
  oob.sex <- data.sex$numeric.sex[which(data.sex$id == i)]
  predicted.sex <- predict(cvfit, oob.meth, type = "class", s = "lambda.min")
  return(c("sex" = oob.sex, "predicted.sex" = predicted.sex, "error" = (as.numeric(predicted.sex) - as.numeric(oob.sex))))
})) %>% data.frame()
loov.sex.correct.assignment <- nrow(y.sex) - sum(abs(as.numeric(loov.sex.res$error)))

# glmnet 10-fold cross volidation cluster prediction
data.cluster <- select(data.complete, -c(sex, decade)) %>% filter(!is.na(cluster.louvain))

x.meth <- as.matrix(select(data.cluster, matches(site.names)))
y.cluster <- as.matrix(data.cluster$cluster.louvain)

cluster.alpha <- sapply(alpha.vals, function(a){
  cvfit <- cv.glmnet(x.meth, y.cluster, alpha = a, family = "multinomial", type.measure = "class")
  error.rate.at.lambdamin <- cvfit$cvm[cvfit$index[1]]
  return(c(a, error.rate.at.lambdamin))
})
cluster.best.alpha <- cluster.alpha[1,which(cluster.alpha[2,] == min(cluster.alpha[2,]))]


# glmnet LOOV Cluster
loov.clust.res <- do.call(rbind,lapply(data.cluster$id, function(i){
  
  dat.loov <- data.cluster[-which(data.cluster$id == i),]
  x.meth <- as.matrix(select(dat.loov, matches(site.names)))
  y.clust <- as.matrix(dat.loov$cluster.louvain)
  
  cvfit <- cv.glmnet(x.meth,y.clust, alpha = cluster.best.alpha, family = "multinomial", type.measure = "class")
  #  corr.coef <- coef(cvfit, s = "lambda.min")
  oob.meth <- as.matrix(filter(data.cluster, id == i) %>%
                          select(matches(site.names)))
  oob.cluster <- data.cluster$cluster.louvain[which(data.cluster$id == i)]
  predicted.cluster <- predict(cvfit, oob.meth, type = "class", s = "lambda.min")
  return(c("cluster" = oob.cluster, "predicted.cluster" = predicted.cluster))
})) %>% data.frame()
loov.clust.errors <- length(which(loov.clust.res$cluster != loov.clust.res$predicted.cluster))
loov.clust.confusion <- table(loov.clust.res)
write.csv(loov.clust.confusion, file = "results/loov.cluster.confusion.matrix.csv")

save(data.complete, sex.alpha, loov.sex.res, cluster.alpha, loov.clust.res, file = "results/glmnet.sex.cluster.rda")

# Random Forest sex classification
data.sex <- select(data.complete, -c(id, cluster.louvain, decade))
freq <- table(data.sex$sex)
sampsize <- rep(ceiling(min(freq / 2)), length(freq))

rf.sex <- rfPermute(
  sex ~ .,
  data.sex,
  sampsize = sampsize,
  replace = FALSE,
  importance = TRUE,
  ntree = 100000,
  keep.forest = FALSE,
  nrep = 1000
)
confusionMatrix(rf.sex)

# Random Forest with predictors identified as significant (p < 0.05) by rfPermute
imp <- importance(rf.sex)
vars2keep <- c("sex",rownames(imp)[which(imp$MeanDecreaseGini.pval < 0.05)])

data.important.predictors <- select(data.sex, all_of(vars2keep))
rf.sex.imp <- rfPermute(
  sex ~ .,
  data.important.predictors,
  sampsize = sampsize,
  replace = FALSE,
  importance = TRUE,
  ntree = 100000,
  keep.forest = FALSE,
  nrep = 1000
)

# Random Forest with predictors chosen by glmnet
load("results/glmnet.sex.model.rda")
glmnet.age.vars[1] <- "sex"
data.glmnet.predictors <- select(data.sex, all_of(glmnet.age.vars))
rf.sex.glmnet.vars <- rfPermute(
  sex ~ .,
  data.glmnet.predictors,
  sampsize = sampsize,
  replace = FALSE,
  importance = TRUE,
  ntree = 100000,
  keep.forest = FALSE,
  nrep = 1000
)

# Random Forest social cluster classification
data.cluster <- select(data.complete, -c(id, sex, decade))
freq <- table(data.cluster$cluster.louvain)
sampsize <- rep(ceiling(min(freq / 2)), length(freq))

rf.cluster <- randomForest(
  cluster.louvain ~ .,
  data.cluster,
  sampsize = sampsize,
  replace = FALSE,
  importance = TRUE,
  ntree = 100000,
  keep.forest = TRUE
)

# Random Forest age decade classification
data.decade <- select(data.complete, -c(id, sex, cluster.louvain))
freq <- table(data.decade$decade)
sampsize <- rep(ceiling(min(freq / 2)), length(freq))

rf.decade <- randomForest(
  decade ~ .,
  data.decade,
  sampsize = sampsize,
  replace = FALSE,
  importance = TRUE,
  ntree = 100000,
  keep.forest = TRUE
)

save(rf.cluster, rf.decade, rf.sex, file = "results/RF.sex.cluster.decade.classification.rda")
