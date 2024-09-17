unlink('~/temp', recursive = TRUE, force = TRUE)
dir.create('~/temp')
Sys.setenv(TEMP = '~/temp', TMP = '~/temp')

library(tidyverse)

# estimate skew-normal parameters from best, min, and max
sn.params <- function(
    age.mode, age.min, age.max, p, shape.const = 5, scale.const = 0.5, 
    maxit = 5000
) {
  age.mean <- mean(c(age.min, age.max))
  shape <- shape.const * (age.mean - age.mode) / (age.mean - age.min)
  scale <- scale.const * (age.max - age.min)
  optim(
    par = c(age.mode, scale, shape), 
    fn = function(dp, obs, p) {
      if(dp[2] <= 0) dp[2] <- .Machine$double.eps 
      abs(obs[1] - swfscMisc::sn.mode(dp)) +
        abs(obs[2] - sn::qsn(1 - p, dp = dp, solver = 'RFB')) +
        abs(obs[3] - sn::qsn(p, dp = dp, solver = 'RFB'))
    }, 
    obs = c(age.mode, age.min, age.max), 
    p = p,
    control = list(
      maxit = maxit, 
      abstol = .Machine$double.eps,
      reltol = 1e-12
    )
  )
}

cwsnDensity <- function(age, age.range, dp, wt) {
  dens.sn <- sn::dsn(age, dp = dp)
  wtd.unif <- (1 - wt) / age.range
  (wtd.unif + (dens.sn * wt))
}

maxCWSNdensity <- function(df) {
  max(sapply(1:nrow(df), function(i) {
    x <- df[i, ]
    cwsnDensity(
      x$age.best, 
      x$age.range,
      dp = c(x$sn.location, x$sn.scale, x$sn.shape),
      wt = x$confidence.wt
    )
  }))
}

rcwsn <- function(n, age.min, age.max, dp, wt) {
  accepted <- numeric(0)
  max.lik <- cwsnDensity(swfscMisc::sn.mode(dp), age.max - age.min, dp, wt)
  while(length(accepted) < n) {
    unif.age <- runif(1, age.min, age.max) 
    if(runif(1, 0, 1) < cwsnDensity(unif.age, age.max - age.min, dp, wt) / max.lik) {
      accepted <- c(accepted, unif.age)
    }
  }
  accepted
}

ranAges <- function(df) {
  sapply(1:nrow(df), function(i) {
    rcwsn(
      n = 1, 
      age.min = df$age.min[i], 
      age.max = df$age.max[i], 
      dp = unlist(df[i, c('sn.location', 'sn.scale', 'sn.shape')]),
      wt = df$confidence.wt[i]
    )
  })
}

crcAgeDist <- function(
    id, age.df,
    n.samples = if(match.arg(type) == 'samples') 10000 else 1000,
    type = c('samples', 'density'),
    from = 0, to = 80
) {
  type <- match.arg(type)
  x <- age.df |> 
    filter(swfsc.id == id) |> 
    as.list()
  
  age <- if(type == 'samples') {
    runif(n.samples * 10, x$age.min, x$age.max)
  } else {
    seq(from, to, length.out = n.samples)
  }
  
  cwsn.dens <- cwsnDensity(
    age = age, 
    age.range = x$age.range, 
    dp = c(x$sn.location, x$sn.scale, x$sn.shape), 
    wt = x$confidence.wt
  )
  
  if(type == 'samples') {
    sample(age[cwsn.dens > 0], n.samples, p = cwsn.dens[cwsn.dens > 0])
  } else {
    data.frame(swfsc.id = id, age = age, density = cwsn.dens) 
  }
}

sampleAgeMeth <- function(ages, meth, num.reps, num.cores) {
  parallel::mclapply(
    1:num.reps, 
    function(rep) {
      meth$ran.meth <- rnorm(nrow(meth), meth$mean.logit, meth$sd.logit)
      ran.df <- ages |>
        filter(swfsc.id %in% ids.to.keep) |>   
        left_join(
          meth |> 
            select(swfsc.id, loc.site, ran.meth) |>
            pivot_wider(names_from = 'loc.site', values_from = 'ran.meth'),
          by = 'swfsc.id'
        ) 
      ran.df$ran.age <- sapply(1:nrow(ran.df), function(i) {
        rcwsn(
          1,
          ran.df$age.min[i],
          ran.df$age.max[i],
          ran.df[i, c('sn.location', 'sn.scale', 'sn.shape')],
          ran.df$confidence.wt[i]
        )
      })
      ran.df |> 
        select(swfsc.id, confidence.wt, ran.age, all_of(sort(unique(meth$loc.site)))) |> 
        as.data.frame()
    },
    mc.cores = num.cores
  )
}


# predictAge <- function(
    #     post, age.df, age.sample, pr.meth.sample, ids, sites, 
#     max.age = 75, 
#     num.samples = min(length(post$deviance), dim(age.sample)[2]), 
#     num.cores = 8) {
#   sex <- setNames(age.df$sex, age.df$swfsc.id)
#   ran.meth <- sample(1:dim(pr.meth.sample)[3], ceiling(sqrt(num.samples) * 2))
#   ran.post <- sample(1:dim(post$exp.age)[2], ceiling(sqrt(num.samples) * 2))
#   
#   age.func <- function(x, post, post.i, sex) {
#     if(is.null(post$w.sex.b)) {
#       post$intercept[post.i] + colSums(post$b.prime[, post.i] * x)
#     } else {
#       post$intercept.prime[sex, post.i] + colSums(post$b.prime[, sex, post.i] * x)
#     }
#   }
#   
#   parallel::mclapply(ids, function(i) {
#     pred.age <- as.vector(apply(
#       qlogis(pr.meth.sample[i, sites, ran.meth]), 2,
#       age.func, post = post, post.i = ran.post, sex = sex[i]
#     ))
#     is.train <- i %in% dimnames(post$exp.age)[[1]]
#     data.frame(
#       swfsc.id = i,
#       type = if(is.train) 'Training' else 'Testing',
#       pred.age = sample(pred.age[dplyr::between(pred.age, 0, max.age)], num.samples),
#       crc.age = sample(age.sample[i, ], num.samples),
#       exp.age = if(is.train) {
#         post$exp.age[i, sample(dim(post$exp.age)[2], num.samples)]
#       } else {
#         rep(NA, num.samples)
#       }
#     )
#   }, mc.cores = num.cores) |> 
#     bind_rows() |> 
#     left_join(select(age.df, swfsc.id, age.confidence), by = 'swfsc.id') |>  
#     mutate(
#       age.confidence = factor(age.confidence),
#       type = factor(type, levels = c('Training', 'Testing'))
#     ) 
# }


postDiag <- function(p, vars, num.cores = 4, sites = NULL, num.samples = 10000) {
  gc(verbose = FALSE)
  library(runjags)
  load('age_and_methylation_data.rdata')
  load('age.pr.meth.samples.rdata')
  sites <- if(is.null(sites)) sites.to.keep else sites
  
  # Plot posterior samples
  tryCatch(plot(
    p$rj, vars = vars,
    file = paste0(p$label, '.plots.', p$end.time, '.pdf'),
    plot.type = c('trace', 'ecdf', 'histogram', 'autocorr')
  ), error = function(e) e)
  
  # # Log-likelihood for loo analysis
  loglik <- swfscMisc::runjags2list(p$rj, collapse.chains = FALSE)$weighted.dens |>
    aperm(c(2, 3, 1)) |>
    log()
  
  smry <- list(
    post.smry = tryCatch(summary(p$rj, vars = vars), error = function(e) NA),
    age.loo = loo::loo(
      loglik,
      r_eff = loo::relative_eff(exp(loglik)),
      cores = num.cores
    ),
    age.pred = predictAge(
      p$post, age.df, age.sample, pr.meth.sample, ids.to.keep, sites, 
      num.samples = num.samples, num.cores = num.cores
    )
  )
  
  save(smry, file = paste0(p$label, '.diagnostics.', p$end.time, '.rdata'))
  invisible(smry)
}


modelData <- function(model.ages, meth, sites, with.sex) {
  model.data <- list(
    num.sites = length(sites),
    freq.meth = meth$freq[model.ages$swfsc.id, sites, drop = FALSE],
    meth.cov = meth$cov[model.ages$swfsc.id, sites, drop = FALSE],
    conversion = meth$conversion[model.ages$swfsc.id, , drop = FALSE],
    zeroes = matrix(0, nrow = nrow(model.ages), ncol = length(sites)),
    num.ind = nrow(model.ages),
    age.min = model.ages$age.min,
    age.max = model.ages$age.max,
    max.age = max(model.ages$age.max),
    age.mean = mean(model.ages$age.best),
    age.tau = 1 / var(model.ages$age.best),
    scale = model.ages$sn.scale,
    shape = model.ages$sn.shape,
    age.offset = model.ages$sn.age.offset,
    dens.unif.1m = model.ages$dens.unif.1m,
    confidence.wt = model.ages$confidence.wt,
    wtd.dens.const = model.ages$wtd.dens.const,
    const = 1.1 * max(sapply(1:nrow(model.ages), function(i) {
      x <- model.ages[i, ]
      crcDensity(
        age = x$age.best,
        age.min = x$age.min,
        age.max = x$age.max,
        dp = c(x$sn.location, x$sn.scale, x$sn.shape),
        wt = x$confidence.wt,
        const = x$wtd.dens.const
      )
    })),
    min.dens = .Machine$double.xmin,
    ones = rep(1, nrow(model.ages))
  )
  if(with.sex) model.data$sex <- model.ages$sex.num
  model.data
}


jags.model.with.sex <- 'model {
  # Prior for intercept
  intercept ~ dnorm(age.mean, age.tau) T(0, max.age)

  # Prior for sex effect on intercept
  w.sex.intercept ~ dbern(0.5)
  
  for(x in 1:2) {
    # Sex effect on intercept
    sex.intercept[x] ~ dnorm(0, 1)
    
    # Weighted intercept
    intercept.prime[x] <- intercept + (sex.intercept[x] * w.sex.intercept)
  }
  
  for(s in 1:num.sites) {
    # Site coefficient
    b[s] ~ dnorm(0, 5e-4)
    
    # Prior for site inclusion
    w.b[s] ~ dbern(0.5)
    
    # Prior for sex effect on coefficient
    w.sex.b[s] ~ dbern(0.5)
    
    for(x in 1:2) {
      # Sex effect on site coefficient
      sex.b[s, x] ~ dnorm(0, 1e-2)
      
      # Weighted coefficient
      b.prime[s, x] <- (b[s] + (sex.b[s, x] * w.sex.b[s])) * w.b[s]
    }
    
    for(i in 1:num.ind) {
      # Probability of methylation for 0 reads
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
      
      b.meth[i, s] <- b.prime[s, sex[i]] * logit(p.meth[i, s])
    }
  }
    
  for(i in 1:num.ind) {
    # Probability of conversion
    p.conv[i] ~ dunif(0, 1)
    conversion[i, 1] ~ dbinom(p.conv[i], conversion[i, 2])
    
    # Linear model predicting age
    exp.age[i] <- intercept.prime[sex[i]] + sum(b.meth[i, ])
    
    # Normal age using predicted age as mode of skew normal
    norm.age[i] <- (exp.age[i] - age.offset[i]) / scale[i]
      
    # Skew normal density, restricted to min and max ages
    dsn[i] <- (2 / scale[i]) * dnorm(norm.age[i], 0, 1) * pnorm(norm.age[i] * shape[i], 0, 1)
    
    # Uniform-weighted density (corrected for confidence) or minimum density
    weighted.dens[i] <- ifelse(
      exp.age[i] >= age.min[i] && exp.age[i] <= age.max[i],
      ((dens.unif.1m[i] + (dsn[i] * confidence.wt[i])) * wtd.dens.const[i]) / const,
      min.dens
    )
    
    # "Ones-trick" likelihood of density 
    # (Kruschke et al 2010, DBDA2E Section 8.6.1 pp. 214-215)
    ones[i] ~ dbern(weighted.dens[i])
    1 ~ dbern(0.3)
    
    exp.age <- ....
    
    norm.dens <- dnorm(exp.age, mu, tau)
    
    ones ~ dbern(norm.dens / const)
  }
}'

monitor.with.sex <- c(
  'deviance', 'weighted.dens', 'exp.age', 'p.meth',
  'intercept.prime', 'w.sex.intercept', 
  'b.prime', 'w.sex.b', 'w.b'
)


jags.model.no.sex <- 'model {
  # Prior for intercept
  intercept ~ dnorm(age.mean, age.tau) T(0, max.age)
  
  for(s in 1:num.sites) {
    # Site coefficient
    b[s] ~ dnorm(0, 5e-4)
    
    # Prior for site inclusion
    w.b[s] ~ dbern(0.5)
      
    # Weighted coefficient
    b.prime[s] <- b[s] * w.b[s]
    
    for(i in 1:num.ind) {
      # Probability of methylation for 0 reads
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
      
      b.meth[i, s] <- b.prime[s] * logit(p.meth[i, s])
    }
  }
    
  for(i in 1:num.ind) {
    # Probability of conversion
    p.conv[i] ~ dunif(0, 1)
    conversion[i, 1] ~ dbinom(p.conv[i], conversion[i, 2])
    
    # Linear model predicting age
    exp.age[i] <- intercept + sum(b.meth[i, ])
    
    # Normal age using predicted age as mode of skew normal
    norm.age[i] <- (exp.age[i] - age.offset[i]) / scale[i]
      
    # Skew normal density, restricted to min and max ages
    dsn[i] <- (2 / scale[i]) * dnorm(norm.age[i], 0, 1) * pnorm(norm.age[i] * shape[i], 0, 1)
    
    # Uniform-weighted density (corrected for confidence) or minimum density
    weighted.dens[i] <- ifelse(
      exp.age[i] >= age.min[i] && exp.age[i] <= age.max[i],
      ((dens.unif.1m[i] + (dsn[i] * confidence.wt[i])) * wtd.dens.const[i]) / const,
      min.dens
    )
    
    # "Ones-trick" likelihood of density 
    # (Kruschke et al 2010, DBDA2E Section 8.6.1 pp. 214-215)
    ones[i] ~ dbern(weighted.dens[i])
  }
}'

monitor.no.sex <- c(
  'deviance', 'weighted.dens', 'exp.age', 'p.meth',
  'intercept', 'b.prime', 'w.b'
)


runModel <- function(
    label, model.ages, meth, sites, with.sex,
    chains, adapt, burnin, total.samples, thin) {
  library(runjags)
  
  p <- list(
    label = label,
    nodename = Sys.info()['nodename'],
    ids = model.ages$swfsc.id,
    sites = sites.to.keep,
    rj = run.jags(
      data = modelData(model.ages, meth, sites, with.sex),
      model = if(with.sex) jags.model.with.sex else jags.model.no.sex,
      monitor = if(with.sex) monitor.with.sex else monitor.no.sex, 
      inits = function() list(
        .RNG.name = 'lecuyer::RngStream',
        .RNG.seed = sample(1:9999, 1)
      ),
      modules = c('glm', 'lecuyer'),
      summarise = FALSE,
      jags.refresh = 30,
      keep.jags.files = FALSE,
      tempdir = FALSE,
      method = 'parallel',
      n.chains = chains,
      adapt = adapt,
      burnin = burnin,
      sample = ceiling(total.samples / chains),
      thin = thin
    )
  )
  p$end.time <- format(Sys.time(), '%y%m%d_%H%M')
  p$elapsed <- swfscMisc::autoUnits(p$rj$timetaken)
  p$post <- swfscMisc::runjags2list(p$rj)
  for(i in 1:length(p$post)) { 
    d <- dim(p$post[[i]])
    if(is.null(d)) next
    j <- which(d == length(p$ids))
    if(length(j) > 0) dimnames(p$post[[i]])[[j]] <- p$ids
    j <- which(d == length(p$sites))
    if(length(j) > 0) dimnames(p$post[[i]])[[j]] <- p$sites
    j <- which(d == 2)
    if(length(j) > 0) dimnames(p$post[[i]])[[j]] <- c('Female', 'Male')
  }
  save(p, file = paste0(p$label, '.posterior.', p$end.time, '.rdata'))
  
  on.gcp <- p$nodename == 'r-instance1'
  postDiag(
    p, c('deviance', 'exp.age'), 
    num.cores = if(on.gcp) 10 else 4, 
    num.samples = total.samples
  )
  
  print(p$elapsed)
  
  # upload data and shutdown instance
  if(on.gcp) {
    system(
      paste0('gsutil cp *.', p$end.time, '.* gs://nmfs-swfscbeam-static/'), 
      timeout = 30
    )
    system('gsutil cp *.Rout gs://nmfs-swfscbeam-static/')
    system('gcloud compute instances stop r-instance1 --zone us-west2-b')
  }
}




# 
# summaryReports <- function(label, folder = 'results') {
#   rmarkdown::render(
#     "summary_MCMC_diagnostics.Rmd", 
#     params = list(folder = folder, timestamp = ts),
#     output_file = file.path(folder, paste0(ts, '_summary_MCMC_diagnostics.pdf'))
#   )
#   
#   rmarkdown::render(
#     "summary_methylation.Rmd", 
#     params = list(folder = folder, timestamp = ts),
#     output_file = file.path(folder, paste0(ts, '_summary_methylation.pdf'))
#   )
#   
#   rmarkdown::render(
#     "summary_site_coefficients.Rmd", 
#     params = list(folder = folder, timestamp = ts),
#     output_file = file.path(folder, paste0(ts, '_summary_site_coefficients.pdf'))
#   )
#   
#   rmarkdown::render(
#     "summary_predict_training.Rmd", 
#     params = list(folder = folder, timestamp = ts),
#     output_file = file.path(folder, paste0(ts, '_summary_predict_training.pdf'))
#   )
#   
#   rmarkdown::render(
#     "summary_predict_testing.Rmd", 
#     params = list(folder = folder, timestamp = ts),
#     output_file = file.path(folder, paste0(ts, '_summary_predict_testing.pdf'))
#   )
# }


predictAgeID <- function(meth.df, id, post, meth.const = 20, max.age = 80, 
                         num.samples = 1e3) {
  id.pred <- NULL
  sites <- dimnames(post$b)[[1]]
  
  while(length(id.pred) < num.samples) {
    lo.meth <- t(sapply(rownames(meth.df), function(i) {
      rnorm(num.samples, meth.df[i, 'logit.meth.mean'], meth.df[i, 'logit.meth.sd'])
    })) + meth.const
    i <- rep(1:length(post$deviance), length.out = num.samples)
    age <- post$intercept[i] + colSums(post$w[sites, i] * post$b[sites, i] * lo.meth[sites, ]) 
    to.keep <- age > 0 & age <= max.age
    age <- age[to.keep]
    num.left <- num.samples - length(id.pred)
    if(length(age) > num.left) age <- sample(age, num.left)
    id.pred <- c(id.pred, age)
  }
  
  id.pred
}


predictAge <- function(age.df, post, meth.params, 
                       meth.const = 20, max.age = 80, num.samples = 1e3,
                       mc.cores = 10) {
  parallel::mclapply(
    ids.to.keep,
    function(id) {
      meth.params |>
        filter(swfsc.id == id) |>
        column_to_rownames('loc.site') |>
        predictAgeID(id, post, meth.const, max.age, num.samples)
    },
    mc.cores = mc.cores
  ) |>
    setNames(ids.to.keep) |>
    bind_cols() |>
    pivot_longer(everything(), names_to = 'swfsc.id', values_to = 'pred.age') |>
    mutate(type = ifelse(swfsc.id %in% p$ids, 'Training', 'Testing')) 
}

plotPredAgeID <- function(pred.age, age.df) {
  g <- pred.age |>
    left_join(age.df, by = 'swfsc.id') |> 
    ggplot() +
    geom_histogram(aes(pred.age, fill = type)) +
    geom_vline(aes(xintercept = age.min), linetype = 'dashed') +
    geom_vline(aes(xintercept = age.max), linetype = 'dashed') +
    geom_vline(aes(xintercept = age.best)) +
    facet_wrap(~ swfsc.id)
  print(g)
}

plotPredCRCAge <- function(pred.age, age.df) {
  g <- pred.age |>
    group_by(swfsc.id, type) |>
    summarize(
      age.mode = modeest::venter(pred.age),
      hdi.lower = unname(HDInterval::hdi(pred.age))[1],
      hdi.upper = unname(HDInterval::hdi(pred.age))[2],
      .groups = 'drop'
    ) |>
    left_join(
      select(age.df, swfsc.id, age.best, age.min, age.max, age.confidence),
      by = 'swfsc.id'
    ) |>
    mutate(age.confidence = factor(age.confidence)) |>
    ggplot() +
    geom_abline(slope = 1, color = "black", linetype = "dashed") +
    geom_segment(
      aes(
        x = age.best,
        xend = age.best,
        y = hdi.lower,
        yend = hdi.upper,
        color = age.confidence
      ),
      linewidth = 0.2,
      show.legend = FALSE
    ) +
    geom_segment(
      aes(
        x = age.min,
        xend = age.max,
        y = age.mode, yend = age.mode,
        color = age.confidence
      ),
      linewidth = 0.2,
      show.legend = FALSE
    ) +
    geom_point(
      aes(x = age.best, y = age.mode, fill = age.confidence),
      shape = 21, size = 3.5, color = 'white'
    ) +
    labs(
      x = "CRC age",
      y = paste("Posterior predicted age"),
      fill = "Confidence",
      title = 'A)'
    ) +
    scale_color_manual(values = conf.colors) +
    scale_fill_manual(values = conf.colors) +
    coord_cartesian(xlim = c(0, 75), ylim = c(0, 75), expand = FALSE) +
    facet_wrap(~ type, ncol = 1)
  
  print(g)
}