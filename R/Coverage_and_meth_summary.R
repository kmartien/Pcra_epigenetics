library(tidyverse)
library(dplyr)
library(ggplot2)
library(gridExtra)

load("results/From.Eric/age_and_methylation_data.rdata")

meth.cov.reduced <- meth.cov[which(rownames(meth.cov) %in% ids.to.keep),]
meth.cov.reduced <- meth.cov.reduced[,which(colnames(meth.cov.reduced) %in% sites.to.keep)]

#median coverage per individual
med.cov.ind <- data.frame(sapply(1:nrow(meth.cov.reduced), function(i){
  median(meth.cov.reduced[i,], na.rm = TRUE)
}))
range(med.cov.ind)
names(med.cov.ind) <- "coverage"

#median coverage per locus
med.cov.site <- data.frame(sapply(1:ncol(meth.cov.reduced), function(i){
  median(meth.cov.reduced[,i], na.rm = TRUE)
}))
names(med.cov.site) <- "coverage"
range(med.cov.site)

# Plot distribution of cover over individuals
p.cov.ind <- ggplot(med.cov.ind, aes(x = coverage)) +
  geom_histogram(color = "black", fill = "grey") +
  labs(x = "Median Coverage", y = "Count", title = "90 Samples", tag = "A") +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5)
  )
p.cov.site <- ggplot(med.cov.site, aes(x = coverage)) +
  geom_histogram(color = "black", fill = "grey") +
  labs(x = "Median Coverage", y = "Count", title = "184 CpG Sites", tag = "B") +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5)
  )
jpeg(file = "results/coverage.hist.90samps.184sites.jpg", width = 900, height = 450)
plots <- list(p.cov.ind, p.cov.site)
plots$nrow <- 1
do.call(grid.arrange,plots)
dev.off()

# Conversion efficiency
range(non.cpg$pct.conv)

# Error rate
epi.df.reduced <- filter(epi.df, id %in% ids.to.keep) %>% filter(loc.site %in% sites.to.keep)
tot.errors <- sum(epi.df.reduced$errors)
tot.reads <- sum(epi.df.reduced$coverage)
error.rate <- tot.errors/tot.reads

# write summaries
cov.mat <- select(epi.df.reduced, c(id, loc.site, coverage)) %>% pivot_wider(names_from = loc.site, values_from = coverage)
write.csv(cov.mat, file = "results/coverage.90samps.184sites.csv")

meth.mat <- select(epi.df.reduced, c(id, loc.site, freq.meth)) %>% pivot_wider(names_from = loc.site, values_from = freq.meth)
write.csv(meth.mat, file = "results/freq.meth.90samps.184sites.csv")

error.mat <- select(epi.df.reduced, c(id, loc.site, errors)) %>% pivot_wider(names_from = loc.site, values_from = errors)
write.csv(error.mat, file = "results/errors.90samps.184sites.csv")

select(non.cpg, c(id, pct.conv)) %>% filter(id %in% ids.to.keep) %>% write.csv(file = "results/pct.conversion.csv")
