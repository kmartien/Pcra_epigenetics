rm(list = ls())
library(tidyverse)
library(mgcv)

load("../age_and_methylation_data.rdata")

model.ages <- age.df |>
  filter(swfsc.id %in% ids.to.keep) |>   
  left_join(
    logit.meth.normal.params |> 
      select(swfsc.id, loc.site, mean.logit) |>
      pivot_wider(names_from = 'loc.site', values_from = 'mean.logit'),
    by = 'swfsc.id'
  )


site.gam <- sapply(sites.to.keep, function(site) {
  fit <- gam(
    as.formula(paste0('age.best ~ s(', site, ', bs = "ts")')),
    data = model.ages,
    weights = model.ages$confidence.wt
  )
  
  site.x <- seq(min(model.ages[[site]]), max(model.ages[[site]]), length.out = 100)
  pred <- predict(
    fit,
    site.x |> 
      data.frame() |> 
      setNames(site),
    type = 'response',
    se.fit = TRUE
  ) |> 
    bind_cols() |> 
    setNames(c('fit', 'se')) |> 
    mutate(
      lower = fit - (2 * se),
      upper = fit + (2 * se),
      logit.meth = site.x
    )
  
  list(site = site, gam = fit, pred = pred)
}, simplify = FALSE)


sapply(site.gam, function(x) summary(x$gam)$s.table[1, 'p-value']) |> 
  enframe(name = 'loc.site', value = 'p.val') |> 
  filter(p.val <= 0.05) |> 
  pull('loc.site') |> 
  saveRDS('gam significant sites.rds')


full.gam <- gam(
  as.formula(
    paste0(
      'age.best ~ ',
      paste0('s(', sites.to.keep, ', bs = "ts")', collapse = ' + ')
    )
  ),
  data = model.ages,
  control = gam.control(nthreads = 10, ncv.threads = 10)
)
full.p.values <- summary(full.gam)$s.table[, 'p-value']
save(full.p.values, full.gam, file = 'full_gam.rdata')


graphics.off()
pdf('gam fit by site.pdf')
for(site in names(site.gam)) {
  g <- logit.meth.normal.params |> 
    mutate(
      meth.lci = qnorm(0.025, mean.logit, sd.logit),
      meth.uci = qnorm(0.975, mean.logit, sd.logit)
    ) |> 
    left_join(
      select(age.df, age.best, age.min, age.max, age.confidence, swfsc.id),
      by = 'swfsc.id'
    ) |> 
    mutate(age.confidence = factor(age.confidence)) |> 
    filter(loc.site == site) |> 
    ggplot() +
    geom_vline(
      aes(xintercept = median(mean.logit)),
      color = 'gray10', alpha = 0.3, linetype = 'dashed'
    ) +
    geom_hline(
      aes(yintercept = median(model.ages$age.best)),
      color = 'gray10', alpha = 0.3, linetype = 'dashed'
    ) +
    geom_ribbon(
      aes(x = logit.meth, ymin = lower, ymax = upper),
      alpha = 0.3,
      data = site.gam[[site]]$pred
    ) +
    geom_line(
      aes(x = logit.meth, y = fit),
      linewidth = 1.5,
      data = site.gam[[site]]$pred
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
    labs(x = 'logit(Pr(meth))', y = 'CRC age', title = site) + 
    theme(
      legend.position = 'top',
      legend.title = element_blank()
    )
  print(g)
}
dev.off()
