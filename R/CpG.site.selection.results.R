library(tidyverse)
library(ggplot2)
library(gridExtra)
library(grid)
source("R/select.RF.important.sites.R")
source("R/select.glmnet.chosen.sites.R")
#load("data/training_set_params.rda")
load("data/color.palettes.rda")


site.selection.params <- read.csv("data-raw/site.selection.params.csv")

site.selection.results <- lapply(1:nrow(site.selection.params), function(i){
  params <- site.selection.params[i,]
  if (params$regression.method == "glmnet") {
    chosen.sites <- select.glmnet.chosen.sites(params = params, incl.prob.threshold = 50)
  } else chosen.sites <- select.RF.important.sites(params = params, incl.prob.threshold = 50)
  return(chosen.sites)
})

site.selection.params$num.CpGs <- do.call(rbind, lapply(site.selection.results, nrow))
write.csv(site.selection.params, file = "results-raw/sites.selection.results.csv")

save(site.seleciton.results, file = "results/site.selection.results.rda")

params2plot <- site.selection.params
params2plot$site.selection.param.set <- as.character(params2plot$site.selection.param.set)
params2plot$min.CR <- as.character(params2plot$min.CR) %>% as.factor()
site.incl.sum <- do.call(rbind, lapply(1:nrow(params2plot), function(i){
  if(params2plot$regression.method[i] == "glmnet") {
    load(paste0("results/glmnet_stability_results-sites", params2plot$site.selection.param.set[i], "_logit_mincov100.rda"))
    site.incl.counts <- select(site.incl.counts, c(site, prob)) %>% bind_cols(params2plot$site.selection.param.set[i])
    site.incl.counts$site <- as.character(site.incl.counts$site)
    return(site.incl.counts)
  } else {
    load(paste0("results/rf_stability_results-sites", params2plot$site.selection.param.set[i], "_logit_mincov100.rda"))
    temp <- filter(importance.sum, pval <= 0.05) %>% group_by(site) %>% 
      summarise(prob = n()) %>% bind_cols(params2plot$site.selection.param.set[i])
  }
}))
site.incl.sum <- filter(site.incl.sum, prob > 10)
names(site.incl.sum)[3] <- "site.selection.param.set"
temp <- pivot_wider(site.incl.sum,id_cols = site, names_from = site.selection.param.set, values_from = prob, values_fill = 0) 
incl.sum <- pivot_longer(temp, cols = !site, names_to = "site.selection.param.set", values_to = "prob") %>%
  left_join(params2plot) 

# bar plot of frequency site had significant p-value
nreps <- 100
g.site.incl.bar.minCR <- lapply(unique(incl.sum$regression.method), function(m){
  dat <- filter(incl.sum, weighting.scheme == "linear") %>% filter(regression.method == m)
  sites2keep <- dat %>% 
    group_by(site) %>% summarise(tot = sum(prob)) %>% filter(tot > 0)
  ggplot(filter(dat, site %in% sites2keep$site), aes(x = site, y = (prob/nreps), fill = min.CR)) +
  geom_bar(stat = "identity", position = "dodge", colour = "black") +
  labs(title = m) +
  scale_fill_manual(values = training.set.palette, name = "Training sample set", 
                    labels = c("ALL", "CR3+", "CR4+")) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
    theme(text = element_text(size = 20))
})

g.site.incl.bar.weight <- lapply(unique(incl.sum$regression.method), function(m){
  dat <- filter(incl.sum, min.CR == 4) %>% filter(regression.method == m)
  sites2keep <- dat %>% 
    group_by(site) %>% summarise(tot = sum(prob)) %>% filter(tot > 0)
  ggplot(filter(dat, site %in% sites2keep$site), aes(x = site, y = (prob/nreps), fill = weighting.scheme)) +
  geom_bar(stat = "identity", position = "dodge", colour = "white") +
  labs(title = m) +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
    theme(text = element_text(size = 20))
})

plots <- c(g.site.incl.bar.minCR, g.site.incl.bar.weight) %>% map(~.x + labs(x=NULL, y=NULL))
bottom <- textGrob("CpG site", gp = gpar(fontsize = 26))
left <- textGrob("Probability included in site selection model", rot = 90, gp = gpar(fontsize = 26))
jpeg(file = paste0("results-raw/site.selection.bar.chart.jpg"), width = 1900, height = 1500)
grid.arrange(grobs = plots, nrow = 4, left = left, bottom = bottom)
dev.off()

g.site.incl.hist <- 
  ggplot(dat) + geom_histogram(aes(x = prob)) +
  facet_wrap(vars(regression.method, min.CR))

# This plots the probability of inclusion for each site, ordered from lowest to highest, from a single run
g.one.run.example <-
  ggplot(temp, aes(x = reorder(site, prob), y = (prob/nreps))) + geom_bar(stat = "identity")
