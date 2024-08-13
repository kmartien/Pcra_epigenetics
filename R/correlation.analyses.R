library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)
source("R/data.prep/combine.age.and.meth.data.R")

min.cov <- 100
meth.type <- "logit" # choose "pct", "pct.no.zero", or "logit"
description <- paste0(meth.type, ".min.cov.",min.cov)

dat <- combine.age.and.meth.data(description)
all.samples <- dat$all.samples
NAs.by.sample <- sapply(1:nrow(all.samples), function(s){
  length(which(is.na(all.samples[s,dat$first.meth.col:ncol(all.samples)])))
})
all.samples <- all.samples[which(NAs.by.sample == 0),]

meth.dat <- select(all.samples, c(id, all_of(dat$site.names)))

loci <- unique(sapply(strsplit(dat$site.names, split = "[.]"), function(i){i[1]}))

# correlation with age for each site
age.corr.coeff <- do.call('rbind',lapply(dat$first.meth.col:length(all.samples), function(site){
  temp <- cor.test(all.samples[[site]],all.samples$age.best,method="pearson", use="na.or.complete")
  return(c(corr.coeff=temp$estimate, p.val = temp$p.value))
})) %>% data.frame()
rownames(age.corr.coeff) <- dat$site.names
age.corr.hist <- hist(abs(age.corr.coeff$corr.coeff.cor), breaks = seq(0, 1, by = 0.1), plot = FALSE)
corr.sum.table <- data.frame(cbind(abs.val.r = age.corr.hist$breaks[1:10], age = age.corr.hist$counts))

# correlation with sex for each site
sex.corr.coeff <- do.call('rbind',lapply(dat$first.meth.col:length(all.samples), function(site){
  temp <- cor.test(all.samples[[site]],all.samples$numeric.sex,method="pearson", use="na.or.complete")
  return(c(corr.coeff=temp$estimate, p.val = temp$p.value))
})) %>% data.frame()
rownames(sex.corr.coeff) <- dat$site.names
sex.corr.hist <- hist(abs(sex.corr.coeff$corr.coeff.cor), breaks = seq(0, 1, by = 0.1), plot = FALSE)
corr.sum.table <- cbind(corr.sum.table, sex = sex.corr.hist$counts)

# correlation with social cluster for each site
all.samples$social.cluster <- as.integer(all.samples$social.cluster)
clust.corr.coeff <- do.call('rbind',lapply(dat$first.meth.col:length(all.samples), function(site){
  temp <- cor.test(all.samples[[site]],all.samples$social.cluster,method="pearson", use="na.or.complete")
  return(c(corr.coeff=temp$estimate, p.val = temp$p.value))
})) %>% data.frame()
rownames(clust.corr.coeff) <- dat$site.names
clust.corr.hist <- hist(abs(clust.corr.coeff$corr.coeff.cor), breaks = seq(0, 1, by = 0.1), plot = FALSE)
corr.sum.table <- cbind(corr.sum.table, social.cluster = clust.corr.hist$counts)

# correlations between sites within a locus
cor.by.loc <- lapply(loci, function(l){
  temp <- cor(meth.dat[-1, grep(l, names(meth.dat))]) %>% data.frame() %>%
    rownames_to_column(var = "loc1")
  temp <- pivot_longer(temp, cols = 2:ncol(temp), names_to = "loc2")
  temp$site.1 <- fct_inorder(do.call(rbind, strsplit(temp$loc1, split = "[.]"))[,2])
  temp$site.2 <- fct_inorder(do.call(rbind, strsplit(temp$loc2, split = "[.]"))[,2])
  return(temp)
})
names(cor.by.loc) <- loci

cor.by.loc.table <- do.call(bind_cols, lapply(cor.by.loc, function(loc){
  corr.hist <- hist(abs(loc$value), breaks = seq(0, 1, by = 0.1), plot = FALSE)
  corr.hist$counts[10] <- corr.hist$counts[10] - sqrt(nrow(loc))
  corr.hist$counts <- corr.hist$counts / (nrow(loc) - sqrt(nrow(loc)))
  return(data.frame(corr.hist$counts))
}))
names(cor.by.loc.table) <- names(cor.by.loc)
corr.sum.table <- cbind(corr.sum.table, cor.by.loc.table)
write.csv(corr.sum.table, file = paste0("results-raw/correlation.summary.", description, ".csv"))

save(age.corr.coeff, sex.corr.coeff, clust.corr.coeff, cor.by.loc, file = paste0("results/correlation.results.", description, ".rda"))

plots <- lapply(1:length(cor.by.loc), function(l){
  ggplot(cor.by.loc[[l]]) + 
  geom_tile(aes(x = site.1, y = site.2, fill = abs(value))) +
  scale_fill_gradientn(colours=brewer.pal(9,"Greys"), breaks = c(0.01,1),
                       labels = c("0", "1"),
                       name = "Abs(correlation coefficient)") +
#  scale_fill_viridis(option = "magma") +
  theme(text = element_text(size = 30),
        legend.key.size = unit(2, 'cm'),
        legend.title = element_text(size = 30),
        axis.text.x=element_blank(), axis.text.y = element_blank())+
  labs(title=loci[l], x="",y="")+
  coord_fixed()
})

plots$nrow <- 2
plots$ncol <- 4
plots$common.legend <- TRUE
plots$legend <- "bottom"
jpeg(file = paste0("results-raw/correlations.by.locus.", description, ".jpg"), width = 1800, height = 1100)
do.call(ggarrange, plots)
dev.off()

