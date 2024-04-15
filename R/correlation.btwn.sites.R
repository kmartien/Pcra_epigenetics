library(dplyr)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)
source("R/data.prep/combine.age.and.meth.data.R")

min.cov <- 100
meth.type <- "pct" # choose "pct", "pct.no.zero", or "logit"
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
}))
rownames(age.corr.coeff) <- dat$site.names

# correlations between sites within a locus
cor.by.loc <- lapply(loci, function(l){
  temp <- cor(meth.dat[-1, grep(l, names(meth.dat))]) %>% data.frame() %>%
    rownames_to_column(var = "loc1")
  pivot_longer(temp, cols = 2:ncol(temp), names_to = "loc2")
})
names(cor.by.loc) <- loci

plots <- lapply(1:length(cor.by.loc), function(l){
  ggplot(cor.by.loc[[l]]) + 
  geom_tile(aes(x = loc1, y = loc2, fill = value))+
  scale_fill_gradientn(colours=brewer.pal(9,"Greys"))+
#  scale_fill_viridis(option = "magma") +
  theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0))+
  labs(title=loci[l], x="",y="")+
  coord_fixed()
})

pdf(file = "results-raw/correlations.by.locus.pdf")
plots
dev.off()
