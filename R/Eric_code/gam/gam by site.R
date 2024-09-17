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
  fit <- NULL
  for(k in 2:20) {
    tryCatch({
    site.fit <- gam(
      as.formula(paste0('age.best ~ te(', site, ', k = ', k, ', bs = "ts")')),
      data = model.ages,
      weights = 1 / model.ages$age.var,
      optimizer = c('outer', 'bfgs')
    )}, error = function(e) NULL)
    if(is.null(site.fit)) break
    if(k.check(site.fit)[1, 'p-value'] > 0.05) {
      fit <- site.fit
      break
    }
  }
  fit
}, simplify = FALSE)

x <- lapply(site.gam, function(z) as.data.frame(k.check(z))) |> 
  bind_rows()
x


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




sig.sites <- sapply(
  site.gam, 
  function(x) {
    smry <- summary(x$gam)
    c(
      p.val = smry$s.table[1, 'p-value'], 
      r.sq = smry$r.sq,
      dev.expl = smry$dev.expl,
      edf = smry$edf,
      dispersion = smry$dispersion,
      residual.df = smry$residual.df
    )
  }
) |> 
  t() |> 
  as.data.frame() |> 
  rownames_to_column('loc.site') |> 
  filter(p.val <= 0.05) |> 
  arrange(desc(r.sq), loc.site)

saveRDS(sig.sites, 'gam significant sites.rds')


graphics.off()
pdf('gam fit by site.pdf')
for(site in sig.sites$loc.site) {
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


full.gam.df <- filter(model.ages, age.confidence >= 4)
full.gam <- gam(
  as.formula(
    paste0(
      'age.best ~ ',
      paste0('s(', sig.sites$loc.site, ', bs = "ts")', collapse = ' + ')
    )
  ),
  data = full.gam.df,
  weights = 1 / full.gam.df$age.var,
  control = gam.control(nthreads = 10, ncv.threads = 10)
)
full.p.values <- summary(full.gam)$s.table[, 'p-value']
names(full.p.values) <- stringr::str_remove(names(full.p.values), "s\\(")
names(full.p.values) <- stringr::str_remove(names(full.p.values), "\\)")

save(full.p.values, full.gam, full.gam.2, file = 'full_gam.rdata')

sites.2 <- names(full.p.values)[full.p.values <= 0.05]
full.gam.2 <- gam(
  as.formula(
    paste0(
      'age.best ~ ',
      paste0('s(', sites.2, ', bs = "ts")', collapse = ' + ')
    )
  ),
  data = full.gam.df,
  weights = 1 / full.gam.df$age.var,
  control = gam.control(nthreads = 10, ncv.threads = 10)
)




