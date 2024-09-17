rm(list = ls())
library(tidyverse)
library(sn)
library(swfscMisc)
library(runjags)
library(parallel)
library(abind)
library(gridExtra)
source("misc_funcs.R")
set.seed(1)


# Load and format age data ------------------------------------------------

crc.data <- read.csv(
  'data/Pseudorca_AgeEstimates_Simplified_2024FEBv2.csv',
  na.strings = c(NA, 'NA', '')
) 


# split out and replicate ids
crc.data <- do.call(
  rbind,
  lapply(1:nrow(crc.data), function(i) {
    id <- strsplit(crc.data$swfsc.id[i], ",")[[1]] |> 
      stringi::stri_trim()
    df <- crc.data[rep(i, length(id)), ]
    df$swfsc.id <- id
    df
  })
)

# format labid
crc.data <- crc.data |> 
  mutate(
    swfsc.id = zero.pad(as.numeric(swfsc.id)),
    swfsc.id = ifelse(is.na(swfsc.id), NA, paste0("z0", swfsc.id))
  )


# format age class
age.class.map <- data.frame(
  class = c("calf", "juvenile", "sub.adult", "adult.f.sm", "adult.pm"),
  f.ages = c(2, 5, 9, max(crc.data$age.best), NA),
  m.ages = c(2, 8, 14, 24, max(crc.data$age.best))
)
crc.data <- crc.data |>  
  mutate(
    age.class = ifelse(
      sex == "Female", 
      as.numeric(cut(age.best, c(0, age.class.map$f), include = T, right = T)),
      as.numeric(cut(age.best, c(0, age.class.map$m), include = T, right = T))
    ),
    age.class = factor(
      age.class.map$class[age.class], 
      levels = age.class.map$class
    )
  )

crc.cols <- c(
  "crc.id", "swfsc.id", "biopsy.id", "date.biopsy", "sex",
  "age.best", "age.confidence", "age.min", "age.max", "age.class", "pair.id"
)
conf.map <- read_csv(
  'data/CRC_confidence_rating_mappings.csv',
  name_repair = 'minimal',
  show_col_types = FALSE
) |> 
  column_to_rownames('person') |> 
  as.matrix()
age.df <- crc.data[, crc.cols] |> 
  mutate( 
    confidence.wt = unname(conf.map["mk", age.confidence]),
    date.biopsy = as.POSIXct(date.biopsy, format = "%m/%d/%y"),
    sex.num = as.numeric(factor(sex))
  )  |> 
  filter(!is.na(age.min)) |> 
  arrange(swfsc.id)
attributes(age.df$date.biopsy)$tzone <- "HST"


# change minimum age of z0178364 to 10 since KMR GLG age was 12.5 and CRC best
# age is 10
age.df$age.min[age.df$swfsc.id == 'z0178364'] <- 10


# fit Skew Normal parameters
age.sn.optim <- mclapply(1:nrow(age.df), function(i) {
  sn.params(
    age.df$age.best[i],
    age.df$age.min[i],
    age.df$age.max[i],
    p = 0.975,
    shape.const = 4
  )
}, mc.cores = 14)
age.sn.params <- t(sapply(age.sn.optim, function(x) x$par))
colnames(age.sn.params) <- c("sn.location", "sn.scale", "sn.shape")
age.df <- cbind(age.df, age.sn.params) |> 
  mutate(
    age.range = age.max - age.min,
    age.var = sn.scale ^ 2 * (1 - ((2 * swfscMisc::sn.delta(sn.shape) ^ 2) / pi)),
    dens.unif = 1 / age.range,
    wtd.unif = (1 - confidence.wt) / age.range,
    sn.age.offset = age.best - (sn.scale * swfscMisc::sn.m0(sn.shape))
  ) 

# compute weighted density scaling constant
# so that weighted density integral = SN integral between min and max
# age.df$wtd.dens.const <- sapply(1:nrow(age.df), function(i) {
#   x <- age.df[i, ]
# 
#   # integral under raw Uniform-weighted Skew Normal
#   wtd.int <- integrate(
#     crcDensity, 
#     lower = -Inf, 
#     upper = Inf,
#     age.min = x$age.min,
#     age.max = x$age.max,
#     dp = c(x$sn.location, x$sn.scale, x$sn.shape),
#     wt = x$confidence.wt,
#     const = 1
#   )$value
#   
#   # constant to multiply by raw Uniform-weighted Skew Normal density 
#   # so integral from min to max integrates to 1
#   1 / wtd.int
# })

age.df <- arrange(age.df, swfsc.id)


# CRC confidence colors ---------------------------------------------------

conf.colors <- setNames(
  colorRampPalette(viridis::magma(30, direction = -1)[-(1:5)])(5),
  5:1
)
conf.colors[5] <- 'gray40'


# Load and format methylation data ----------------------------------------

load("data/Pcra.epi.data.for.Eric.Rdata")
load("data/res.sum.Rdata")

epi.df <- do.call(
  rbind,
  lapply(names(res.sum), function(locus) {
    do.call(
      rbind,
      lapply(names(res.sum[[locus]]), function(type) {
        x <- res.sum[[locus]][[type]]
        if(is.null(x)) return(NULL)
        x$coverage |> 
          select(-avg) |> 
          pivot_longer(-id, names_to = "site", values_to = "coverage") |> 
          left_join(
            pivot_longer(x$freq.meth, -id, names_to = "site", values_to = "freq.meth"),
            by = c("id", "site")
          ) |> 
          left_join(
            pivot_longer(x$errors, -id, names_to = "site", values_to = "errors"),
            by = c("id", "site")
          ) |>
          mutate(locus = locus, type = type)
      })
    )
  })
) |>   
  mutate(
    site = as.numeric(site),
    loc.site = paste0(locus, "_", zero.pad(site)),
    id.site = paste0(id, "_", loc.site),
    corrected.cov = coverage - errors,
    pct.meth = freq.meth / corrected.cov,
    type = gsub(".sum", "", type)
  ) |> 
  arrange(type, id, locus, site) 

# create map data frame associating loci and sites
locus.map <- epi.df |> 
  select(loc.site, locus, site) |> 
  unique() |> 
  mutate(locus.num = as.numeric(factor(locus))) |> 
  arrange(locus, site) |> 
  group_by(locus) |> 
  mutate(locus.site.num = as.numeric(factor(site))) |> 
  ungroup()

# create map vector associating loci and sites
locus <- locus.map |> 
  select(loc.site, locus.num) |> 
  deframe()


# Calculate conversion rate -----------------------------------------------

non.cpg <- epi.df |>
  filter(type == "non.CpG") |>
  group_by(id) |>
  summarize(
    coverage = sum(corrected.cov, na.rm = TRUE),
    freq.meth = sum(freq.meth, na.rm = TRUE),
    .groups = "drop"
  ) |> 
  mutate(
    converted = coverage - freq.meth,
    pct.conv = converted / coverage,
    conv.shape1 = converted + 1,
    conv.shape2 = (coverage - converted) + 1
  )

cpg <- epi.df |> 
  filter(type == "CpG") |> 
  mutate(
    meth.shape1 = ifelse(
      is.na(freq.meth) | is.na(corrected.cov), 
      1, freq.meth + 1
    ),
    meth.shape2 = ifelse(
      is.na(freq.meth) | is.na(corrected.cov), 
      1, (corrected.cov - freq.meth) + 1
    )
  ) |> 
  left_join(
    select(non.cpg, id, pct.conv, conv.shape1, conv.shape2),
    by = "id"
  )


# Plot age distributions --------------------------------------------------

graphics.off()
pdf(paste0("results/summary - CRC age priors.pdf"), height = 10, width = 10)

p <- age.df %>%
  mutate(
    id = factor(swfsc.id, levels = .$swfsc.id[order(-age.best, -age.min, -age.max)]),
    Confidence = factor(age.confidence)
  ) |> 
  ggplot(aes(y = id)) +
  geom_segment(
    aes(x = age.min, xend = age.max, yend = id, color = Confidence), 
    linewidth = 2, alpha = 0.8
  ) +
  geom_point(
    aes(x = age.best, color = Confidence), 
    shape = 21, size = 3.5, fill = "white"
  ) + 
  scale_color_manual(values = conf.colors) +
  labs(x = "Age", y = NULL) + 
  scale_x_continuous(breaks = seq(0, max(age.df$age.max), 5)) +
  theme(
    text = element_text(size = 10),
    legend.position = "top"
  )
print(p)

p <- age.df |> 
  mutate(Confidence = factor(age.confidence)) |> 
  ggplot(aes(age.best, age.range)) +
  geom_point(aes(fill = Confidence), color = "white", shape = 21, size = 4) +
  scale_fill_manual(values = conf.colors) +
  labs(x = "CRC best age", y = "Maximum - minimum age") +
  theme(legend.position = "top")
print(p)


for(cr.df in split(age.df, age.df$age.confidence)) {
  conf <- unique(cr.df$age.confidence)
  p <- lapply(1:nrow(cr.df), function(i) {
    crcAgeDist(cr.df$swfsc.id[i], cr.df, type = 'density')
  }) |>
    bind_rows() |>
    ggplot() +
    geom_vline(aes(xintercept = age.min), linetype = 'dashed', data = cr.df) +
    geom_vline(aes(xintercept = age.max), linetype = 'dashed', data = cr.df) +
    geom_vline(aes(xintercept = age.best), data = cr.df) +
    geom_line(
      aes(x = age, y = density),
      color = conf.colors[as.character(conf)],
      linewidth = 1.5
    ) +
    labs(x = 'Age', y = 'Density', title = paste0("Confidence rating = ", conf)) +
    facet_wrap(~ swfsc.id, scales = 'free_y') +
    theme_minimal() +
    theme(
      axis.line.x = element_line(),
      axis.line.y = element_line()
    )
  print(p)
}

dev.off()


# Extract conversion, methylation, and coverage matrices ----------------

meth <- list(
  freq = cpg |> 
    select(id, loc.site, freq.meth) |> 
    pivot_wider(names_from = loc.site, values_from = freq.meth) |> 
    column_to_rownames("id") |> 
    as.matrix(),
  cov = cpg |> 
    select(id, loc.site, corrected.cov) |> 
    pivot_wider(names_from = loc.site, values_from = corrected.cov) |> 
    column_to_rownames("id") |> 
    as.matrix(),
  conversion = non.cpg |> 
    select(id, converted, coverage) |> 
    arrange(id) |> 
    column_to_rownames("id") |> 
    as.matrix()
)
meth$obs.pct <- meth$freq / meth$cov
meth$est.pct <- 1 - ((1 - meth$obs.pct) / (meth$conversion[, 1] / meth$conversion[, 2]))
meth$est.pct[meth$est.pct < 0] <- 0


# IDs with good coverage and no missing sites -----------------------------

ids.to.keep <- cpg |> 
  filter(id %in% age.df$swfsc.id & type == "CpG") |> 
  group_by(id) |> 
  summarize(
    is.complete = !any(is.na(freq.meth)),
    median.cov = median(coverage, na.rm = TRUE), 
    .groups = "drop"
  ) |> 
  filter(median.cov >= 1000 & is.complete) |> 
  pull("id")


# Summarize and select sites ----------------------------------------------

site.smry <- cpg |> 
  filter(id %in% ids.to.keep & type == "CpG") |> 
  mutate(
    lo.meth = swfscMisc::logOdds(freq.meth / corrected.cov),
    lo.meth = ifelse(is.infinite(lo.meth), NA, lo.meth)
  ) |> 
  left_join(age.df, by = c(id = "swfsc.id")) |> 
  group_by(loc.site) |> 
  summarize(
    age.cor = cor(lo.meth, age.best, use = "complete.obs"),
    lo.var = var(lo.meth, na.rm = TRUE),
    cov.var = var(corrected.cov, na.rm = TRUE),
    cov.median = median(corrected.cov, na.rm = TRUE),
    .groups = "drop"
  ) |> 
  arrange("loc.site")

sites.to.keep <- site.smry |> 
  filter(cov.median >= 1000) |> 
  pull("loc.site")



# Paired sample age difference --------------------------------------------

pair.df <- age.df |>
  filter(!is.na(age.df$pair.id) & swfsc.id %in% ids.to.keep)

num.samples <- 1000000
pair.df <- split(pair.df, pair.df$pair.id) |>
  lapply(function(df) {
    if(n_distinct(df$swfsc.id) < 2) NULL else {
      combn(df$swfsc.id, 2, function(id) {
        pairs <- data.frame(swfsc.id = id, pair.id = unique(df$pair.id)) |>
          left_join(select(age.df, swfsc.id, date.biopsy), by = 'swfsc.id') |>
          arrange(desc(date.biopsy))
        
        older <- pairs$swfsc.id[1]
        younger <- pairs$swfsc.id[2]
        
        older.sn.dp <- age.df |>
          filter(swfsc.id == older) |>
          select(sn.location, sn.scale, sn.shape) |>
          unlist()
        younger.sn.dp <- age.df |>
          filter(swfsc.id == younger) |>
          select(sn.location, sn.scale, sn.shape) |>
          unlist()
        pair.diff <- sn::rsn(num.samples, dp = older.sn.dp) - sn::rsn(num.samples, dp = younger.sn.dp)
        
        data.frame(
          older = older,
          younger = younger,
          pair.id = unique(df$pair.id),
          year.diff = as.numeric(difftime(pairs$date.biopsy[1], pairs$date.biopsy[2], 'days')) / 365.25,
          mean = mean(pair.diff),
          sd = sd(pair.diff)
        )
      }, simplify = FALSE) |>
        bind_rows()
    }
  }) |>
  bind_rows()


# Estimate true methylation -----------------------------------------------

post <- run.jags(
  model = "model {
    for(i in 1:num.ind) {
      # Probability of conversion
      p.conv[i] ~ dunif(0, 1) 
      conversion[i, 1] ~ dbinom(p.conv[i], conversion[i, 2])
      
      for(s in 1:num.sites) {
        # Probability of minimum amount of methylation (0 reads)
        p.zero.meth[i, s] ~ dunif(0, 1)
        zeroes[i, s] ~ dbinom(p.zero.meth[i, s], meth.cov[i, s])
      
        # Probability of observed methylation
        obs.p.meth[i, s] ~ dunif(0, 1)
        freq.meth[i, s] ~ dbinom(obs.p.meth[i, s], meth.cov[i, s])
        
        # Probability of actual methylation
        p.meth[i, s] <- ifelse(
          (1 - obs.p.meth[i, s]) >= p.conv[i], 
          p.zero.meth[i, s], 
          1 - ((1 - obs.p.meth[i, s]) / p.conv[i])
        )
      }
    }
  }",
  monitor = c("deviance", "p.meth"), 
  data = list(
    num.ind = length(ids.to.keep),
    num.sites = length(sites.to.keep),
    freq.meth = meth$freq[ids.to.keep, sites.to.keep],
    meth.cov = meth$cov[ids.to.keep, sites.to.keep],
    conversion = meth$conversion[ids.to.keep, ],
    zeroes = matrix(0, nrow = length(ids.to.keep), ncol = length(sites.to.keep))
  ),
  inits = function() list(
    .RNG.name = "lecuyer::RngStream",
    .RNG.seed = sample(1:9999, 1)
  ),
  modules = c("glm", "lecuyer"),
  summarise = FALSE,
  method = "parallel",
  n.chains = 10,
  adapt = 100,
  burnin = 1e3,
  sample = ceiling(1e4 / 10),
  thin = 10
)
post$timetaken <- swfscMisc::autoUnits(post$timetaken)
save(post, file = "results/pr.meth.posterior.rdata")

p <- swfscMisc::runjags2list(post)
dimnames(p$p.meth) <- list(
  swfsc.id = ids.to.keep, 
  loc.site = sites.to.keep, 
  rep = 1:dim(p$p.meth)[3]
)
logit.meth.normal.params <- qlogis(p$p.meth) |> 
  as.data.frame.table(responseName = 'logit.meth', stringsAsFactors = FALSE) |> 
  group_by(swfsc.id, loc.site) |> 
  summarize(
    mean.logit = mean(logit.meth),
    sd.logit = sd(logit.meth),
    .groups = 'drop'
  )


# Save data ---------------------------------------------------------------

save(
  age.df, conf.colors, epi.df, cpg, non.cpg, meth,
  locus.map, locus, site.smry,
  ids.to.keep, sites.to.keep, pair.df,
  logit.meth.normal.params,
  file = "age_and_methylation_data.rdata"
)



df <- logit.meth.normal.params |> 
  mutate(
    meth.lci = qnorm(0.025, mean.logit, sd.logit),
    meth.uci = qnorm(0.975, mean.logit, sd.logit)
  ) |> 
  left_join(
    select(age.df, age.best, age.min, age.max, age.confidence, swfsc.id),
    by = 'swfsc.id'
  ) |> 
  mutate(age.confidence = factor(age.confidence))

graphics.off()
pdf('results/summary - logit meth by age.pdf')
for(i in sort(sites.to.keep)) {
  g <- df |> 
    filter(loc.site == i) |> 
    ggplot() +
    geom_vline(
      aes(xintercept = median(mean.logit)), 
      color = 'gray10', alpha = 0.5, linetype = 'dashed'
    ) +
    geom_hline(
      aes(yintercept = median(age.df$age.best)),
      color = 'gray10', alpha = 0.5, linetype = 'dashed'
    ) +
    geom_segment(
      aes(x = meth.lci, xend = meth.uci, y = age.best, yend = age.best, color = age.confidence),
      alpha = 0.6, linewidth = 0.2
    ) +
    geom_segment(
      aes(x = mean.logit, xend = mean.logit, y = age.min, yend = age.max, color = age.confidence),
      alpha = 0.6, linewidth = 0.2
    ) +
    geom_point(
      aes(mean.logit, age.best, fill = age.confidence), 
      color = 'white', shape = 21, size = 3
    ) +
    scale_fill_manual(values = conf.colors) +
    scale_color_manual(values = conf.colors) +
    labs(x = 'logit(Pr(meth))', y = 'CRC age', title = i) + 
    theme(
      legend.position = 'top',
      legend.title = element_blank()
    )
  print(g)
}
dev.off()
