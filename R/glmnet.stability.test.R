library(glmnet)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(viridis)
source("R/data.prep/combine.age.and.meth.data.R")
load("data/training.set.params.rda")

nreps <- 100
min.cov <- 100
meth.type <- "logit" # "logit" or "pct"

training.sets <- filter(training.set.params, sites2use == "all")
lapply(1:nrow(training.sets), function(t){
  description <- paste0("train", training.sets$id[t], "_", meth.type, "_mincov",min.cov)
  weight.type <- training.sets$weight[t]
  training.min.CR <- training.sets$minCR[t]
  #weighted <- TRUE
  #weight.type <- "CRC" # "CRC" or "linear"
  #training.min.CR <- 4
  
  dat <- combine.age.and.meth.data(paste0(meth.type, "_mincov",min.cov))
  site.names <- dat$site.names
  all.samples <- dat$all.samples
  if (weight.type == "linear") all.samples$wt <- all.samples$age.confidence/5
  if (weight.type == "none") all.samples$wt <- 1
  
  calibration.set <- subset(all.samples, subset=(all.samples$age.confidence >= training.min.CR))
  
  alpha.vals <- seq(from = 0.1, to = 0.9, by = 0.1)
  
  #alternately remove sites and samples to eliminate all NAs
  
  calibration.set.complete <- calibration.set
  
  NAs.by.sample <- sapply(1:nrow(calibration.set.complete), function(s){
    length(which(is.na(calibration.set.complete[s,dat$first.meth.col:ncol(calibration.set.complete)])))
  })
  
  calibration.set.complete <- calibration.set.complete[which(NAs.by.sample == 0),]
  
  #################################################
  # Cross-validated linear age fit of age over a range of alphas
  
  x.meth <- as.matrix(select(calibration.set.complete, matches(site.names)))
  y.age <- as.matrix(calibration.set.complete$age.best)
  wt <- calibration.set.complete$wt

  print(paste0("Training set", training.sets$id[t]))  
  cvfit.list <- lapply(1:nreps, function(i){
    print(i)
    lapply(alpha.vals, function(a){
      cv.glmnet(x.meth, y.age, weights = wt, alpha = a)
    })
  })
  lambda.name <- c("lambda.min", "lambda.1se")
  date()
  alpha.test <- do.call('rbind', lapply(1:nreps, function(i){
    test.alpha <- data.frame(do.call('rbind', lapply(1:length(alpha.vals), function(j){
      #    cvfit <- cv.glmnet(x.meth,y.age, alpha = a)
      cvfit <- cvfit.list[[i]][[j]]
      a <- alpha.vals[j]
      res <- do.call('rbind', lapply(1:2, function(l){
        lambda <- lambda.name[l]
        lambda.val <- cvfit$lambda[cvfit$index[l]]
        corr.coef <- coef(cvfit, s = lambda)
        included.sites <- as.array(corr.coef)
        names(included.sites) <- corr.coef@Dimnames[[1]]
        glmnet.predicted.age <- predict(cvfit, x.meth, s = lambda)
        regr.model <- lm(glmnet.predicted.age ~ y.age)
        regr.predicted.age <- predict.lm(regr.model)
        mean.error <- mean(abs(y.age - glmnet.predicted.age))
        median.error <- median(abs(y.age - glmnet.predicted.age))
        mse <- cvfit$cvm[cvfit$index[l]]
        num.predictors <- length(which(corr.coef != 0)) - 1 #don't count the y-intercept
        return(c(iter = i, lambda.type = l, lambda.val = lambda.val, mean.error = mean.error, median.error = median.error, mse = mse, num.predictors = num.predictors, included.sites))
      }))
      res <- cbind(alpha = a, res)
      return(res = res)
    })))
  }))
  
  chosen.model <- do.call('rbind', lapply(1:nreps, function(i){
    test.alpha <- filter(alpha.test, iter == i)
    res <- data.frame(do.call('rbind', lapply(1:2, function(l){
      data.frame(do.call('rbind',lapply(c("mean.error", "median.error", "mse"), function(e){
        temp <- filter(test.alpha, lambda.type == l) %>% select(c(lambda.type, all_of(e), alpha, num.predictors))
        temp <- temp[which(temp[,2] == min(temp[,2])),]
        temp$error.measure <- e
        names(temp)[2] <- "error"
        return(temp)
      })))
    })))
    res$iter <- i
    return(res)
  }))
  date()
  
  # site inclusion
  chosen.alphas <- filter(chosen.model, lambda.type == 1) %>% 
    filter(error.measure == "median.error") %>%
    select(c(iter, alpha))
  site.incl.counts <- 
    filter(alpha.test, lambda.type == 1) %>%
    right_join(chosen.alphas) %>%
    select(all_of(site.names)) %>%
    summarise_all(~sum(. != 0)) %>%
    t() %>% data.frame() %>% 
    rownames_to_column(var = "site")
  site.incl.counts <- bind_cols(site.incl.counts, 
                                do.call(rbind, strsplit(site.incl.counts$site, split = "[.]")))
  names(site.incl.counts) <- c("site", "prob", "locus", "pos")
  site.incl.counts <- filter(site.incl.counts, prob > 0)
  x.labs <- do.call(rbind, lapply(unique(site.incl.counts$locus), function(x){
    indices <- which(site.incl.counts$locus == x)
    breaks <- site.incl.counts$site[min(indices)]
    mid.points <- site.incl.counts$site[indices[(floor(length(indices)/2) + 1)]]
    return(c(breaks = breaks, mid.points = mid.points))
  })) %>% data.frame()
  site.incl.counts$site <- fct_inorder(site.incl.counts$site)
  
  save(cvfit.list, alpha.test, chosen.model, calibration.set.complete, site.incl.counts, file = paste0("results/glmnet_stability_results-", description,".rda"))
  
  g.site.incl.bar <- 
    ggplot(site.incl.counts, aes(x = site, y = (prob/nreps))) +
    geom_bar(stat = "identity") +
    scale_x_discrete(breaks = x.labs$breaks, 
                     labels = c("DIRAS3", "FLJ0945", "GRIA2", "HMG20B", "KCNC4", 
                                "PRDM12","TET2", "VGF"), 
                     guide = guide_axis(angle = 90)) +
    labs(x = "Site", y = "Inclusion Probability")
  jpeg(file = paste0("results-raw/glmnet_site_inclusion-", description, ".jpg"), width = 900, height = 500)
  g.site.incl.bar
  dev.off()
})
# 
# best.alpha.hist <- lapply(1:2, function(l){
#   ggplot(filter(chosen.model, lambda.type == l)) + 
#     facet_wrap(~error.measure) +
#     geom_histogram(aes(x = alpha, fill = error.measure), binwidth = 0.1, col = "grey25") +
#     ggtitle(lambda.name[l])
# })
# best.alpha.hist$nrow = 2
# 
# num.predictors.hist <- lapply(1:2, function(l){
#   ggplot(filter(chosen.model, lambda.type == l)) + 
#     facet_wrap(~error.measure) +
#     geom_histogram(aes(x = num.predictors, fill = error.measure), col = "grey25") +
#     ggtitle(lambda.name[l])
# })
# 
# num.predictors.hist$nrow = 2
# 
# dat <- filter(chosen.model, error.measure %in% c("mean.error", "median.error"))
# mae.hist <- lapply(1:2, function(l){
#   #  ggplot(filter(chosen.model, lambda.type == l)) + 
#   ggplot(filter(dat, lambda.type == l)) + 
#     facet_wrap(~error.measure) +
#     geom_histogram(aes(x = error, fill = error.measure), col = "grey25") +
#     ggtitle(lambda.name[l])
# })
# mae.hist$nrow = 2
# 
# jpeg(paste0("results-raw/Best.alpha.histograms.minCR.", training.min.CR, ".", description, ".jpg"), height = 1000, width = 700)
# do.call(grid.arrange, best.alpha.hist)
# dev.off()
# 
# jpeg(paste0("results-raw/num.predictors.histograms.minCR.", training.min.CR, ".", description, ".jpg"), height = 1000, width = 700)
# do.call(grid.arrange, num.predictors.hist)
# dev.off()
# 
# jpeg(paste0("results-raw/mae.histograms.zoomed.minCR.", training.min.CR, ".", description, ".jpg"), height = 1000, width = 700)
# do.call(grid.arrange, mae.hist)
# dev.off()
# 
# # for a given alpha, look at how lambda-min, mse, and 
# # num.predictors varies between runs
# 
# dat <- filter(alpha.test, lambda.type == 1) 
# pdf(paste0("results-raw/alpha_variability.minCR.", training.min.CR, ".", description, ".pdf"), width = 6)
# for (a in alpha.vals){
#   lambda.hist <- ggplot(filter(dat, alpha == a)) + 
#     geom_histogram(aes(x = lambda.val), fill = "#F8766D", col = "grey25") +
#     ggtitle(paste0("Alpha = ", a))
#   mse.hist <- ggplot(filter(dat, alpha == a)) + 
#     geom_histogram(aes(x = mse), fill = "#00BA38", col = "grey25") 
#   predictors.hist <- ggplot(filter(dat, alpha == a)) + 
#     geom_histogram(aes(x = num.predictors),fill = "#619CFF", col = "grey25")
#   #  jpeg(paste0("results-raw/alpha_", a, "_variability.jpg"), height = 1000, width = 700)
#   plots <- list(lambda.hist, mse.hist, predictors.hist)
#   plots$nrow <- 3
#   do.call(grid.arrange, plots)
#   #  dev.off()
# }
# dev.off()
# 
# x.labs <- c("DIRAS3.102", "FLJ0945.102", "GRIA2.102", "HMG20B.104", "KCNC4.100", 
#             "PRDM12.105","TET2.107", "VGF.107")
# site.incl.counts <- do.call(rbind, lapply(alpha.vals, function(a){
#   filter(alpha.test, lambda.type == 1) %>%
#     filter(alpha == a) %>%
#     select(all_of(site.names)) %>%
#     summarise_all(~sum(. != 0)) %>%
#     add_column(alpha = a) %>% relocate(alpha)
# }))
# site.incl.counts <- site.incl.counts[,which(colSums(site.incl.counts)>0)]
# save(cvfit.list, alpha.test, chosen.model, site.incl.counts, calibration.set.complete, file = paste0("results/glmnet.stability.results.minCR",training.min.CR, ".", description,".rda"))
# 
# 
# pdf(paste0("results-raw/site.inclusion.by.alpha.minCR.", training.min.CR, ".", description, ".pdf"), width = 10, height = 7)
# g.site.incl.bar <- pivot_longer(site.incl.counts, cols = -alpha) %>%
#   ggplot(., aes(x = name, y = value, fill = as.character(alpha))) +
#   geom_bar(stat = "identity", position = position_dodge()) +
#   scale_x_discrete(guide = guide_axis(angle = 90)) +
#   labs(fill = "Alpha", x = "Site", y = "Count")
# grid.arrange(g.site.incl.bar, ncol = 1)
# dev.off()
# 
# dat <- pivot_longer(select(calibration.set.complete, c("id","age.point", all_of(site.names))),
#                     all_of(site.names), names_to = "site", values_to = "meth")
# temp <- data.frame(t(site.incl.counts[5,-1])) 
# temp <- rownames_to_column(temp)
# names(temp) <- c("site","prob.chosen")
# dat <- left_join(dat, temp) %>% filter(prob.chosen > 0)
# conf.palette <- colorRampPalette(viridis::magma(30, direction = -1))
# g.for.eric <- ggplot(dat, aes(x = age.point, y = meth, col = prob.chosen)) +
#   geom_point(position = "jitter") +
#   ylab(paste0(meth.type, ".meth")) +
#   scale_color_viridis()
# jpeg(paste0("results-raw/plot.eric.asked.for.minCR.",training.min.CR, ".", description,".jpg"), width = 2000, height = 1000)
# g.for.eric 
# dev.off()
# 
# g.site.incl.bar <- filter(pivot_longer(corrected.meth, cols = -id), name %in% names(site.incl.counts)) %>%
#   ggplot(., aes(x = name, y = value)) +
#   geom_bar(stat = "identity")
