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

sampleAgeMeth <- function(ages, meth, ran.age = FALSE, ran.meth = FALSE) {
  ages |>  
    mutate(age.ran = ranAges(ages)) |> 
    left_join(
      meth |> 
        mutate(ran.meth = rnorm(n(), mean.logit, sd.logit)) |> 
        select(swfsc.id, loc.site, ran.meth) |>
        pivot_wider(names_from = 'loc.site', values_from = 'ran.meth'),
      by = 'swfsc.id'
    ) 
}


# fit GAM model to training data
fitTrainGAM <- function(df, sites, resp) {
  # run GAM starting at k = 3 until convergence
  fit <- NULL
  for(k in 3:10) {
    fit <- tryCatch(
      gam(
        formula = as.formula(paste(
          resp, '~',
          paste0('te(', sites, ', k = ', k, ', bs = "ts")', collapse = ' + ')
        )),
        data = df,
        select = TRUE
      ),
      error = function(e) NULL
    )
    if(!is.null(fit)) break
  }
  if(is.null(fit)) NULL else list(gam = fit, k = k)
}

# predict testing data from fitted model
predictTestGAM <- function(fit, cv.df, resp){
  if(is.null(fit)) return(NULL)
  
  pred <- predict(
    fit$gam,
    newdata = cv.df, 
    type = 'response',
    se.fit = TRUE
  )
  
  tibble(
    swfsc.id = cv.df$swfsc.id,
    k = fit$k,
    age.resp = cv.df[[resp]],
    age.pred = unname(ifelse(pred$fit < 0, 0, pred$fit)),
    lci = unname(pred$fit + (pred$se.fit * qnorm(0.025))),
    uci = unname(pred$fit + (pred$se.fit * qnorm(0.975)))
  ) |> 
    mutate(
      lci = ifelse(lci < 0, 0, lci),
      uci = ifelse(uci < 0, 0, uci),
      # I think resid, dev, and sq.err should be calculated based on age.best, not age.resp
      resid = age.pred - age.resp,
      dev = abs(resid),
      sq.err = resid ^ 2
    )
}

predictAllIDsGAM <- function(train.df, model.df, sites, resp) {
  rbind( 
    # cross-validation model for CR 4 & 5
    lapply(train.df$swfsc.id, function(cv.id) {
      fitTrainGAM(filter(train.df, swfsc.id != cv.id), sites, resp) |> 
        predictTestGAM(filter(model.df, swfsc.id == cv.id), resp)
    }) |> 
      bind_rows(),
    # full CR 4 & 5 model to predict CR 2 & 3
    fitTrainGAM(train.df, sites, resp) |> 
      predictTestGAM(filter(model.df, age.confidence %in% 2:3), resp)
  ) 
}

fitTrainSVM <- function(df, sites, resp, tune.obj) {
  fit <- svm(
    formula = as.formula(paste0(resp, ' ~ .')), 
    data = select(df, c(resp, all_of(sites))),
    cost = tune.obj$best.parameters$cost,
    gamma = tune.obj$best.parameters$gamma)
}


predictTestSVM <- function(fit, cv.df, sites, resp){
  if(is.null(fit)) return(NULL)
  
  pred <- predict(
    fit, 
    select(cv.df, all_of(sites))
  )
  
  tibble(
    swfsc.id = cv.df$swfsc.id,
    age.resp = cv.df[[resp]],
    age.pred = unname(ifelse(pred < 0, 0, pred)),
  ) |> 
    mutate(
      resid = age.pred - age.resp,
      dev = abs(resid),
      sq.err = resid ^ 2
    )
}

predictAllIDsSVM <- function(train.df, model.df, sites, resp, tune.obj) {
  rbind( 
    # cross-validation model for CR 4 & 5
    lapply(train.df$swfsc.id, function(cv.id) {
      fitTrainSVM(filter(train.df, swfsc.id != cv.id), sites, resp, tune.obj) |> 
        predictTestSVM(filter(model.df, swfsc.id == cv.id), sites, resp)
    }) |> 
      bind_rows(),
    # full CR 4 & 5 model to predict CR 2 & 3
    fitTrainSVM(train.df, sites, resp, tune.obj) |> 
      predictTestSVM(filter(model.df, age.confidence %in% 2:3), sites, resp)
  ) 
}

fitTrainRF <- function(df, sites, resp, rf.params) {
  fit <- randomForest(
    y = df$resp,
    x = select(df, all_of(sites)),
    mtry = rf.params$mtry,
    ntree = 10000,
    sampsize = rf.params$sampsize,
    replace = FALSE
  )
  
  tibble(
    swfsc.id = df$swfsc.id,
    age.resp = cv.df[[resp]],
    age.pred = fit$predicted,
  ) |> 
    mutate(
      resid = age.pred - age.resp,
      dev = abs(resid),
      sq.err = resid ^ 2
    )
}


predictTestRF <- function(fit, cv.df, sites, resp){
  if(is.null(fit)) return(NULL)
  
  pred <- predict(
    fit, 
    select(cv.df, all_of(sites))
  )
  
  tibble(
    swfsc.id = cv.df$swfsc.id,
    age.resp = cv.df[[resp]],
    age.pred = unname(ifelse(pred < 0, 0, pred)),
  ) |> 
    mutate(
      resid = age.pred - age.resp,
      dev = abs(resid),
      sq.err = resid ^ 2
    )
}

predictAllIDsRF <- function(train.df, model.df, sites, resp, rf.params) {
  fit <- randomForest(
    formula = as.formula(paste0(resp, ' ~ .')), 
    data = select(train.df, c(resp, all_of(sites))),
    #y = train.df$resp,
    #x = select(train.df, all_of(sites)),
    mtry = rf.params$mtry,
    ntree = 10000,
    sampsize = rf.params$sampsize,
    replace = FALSE
  )
  
  test.df <- filter(model.df, age.confidence %in% 2:3)
  pred <- predict(
    fit, 
    select(test.df, all_of(sites))
  )
  
  rbind( 
    # OOB for CR 4 & 5
    tibble(
      swfsc.id = train.df$swfsc.id,
      age.resp = train.df[[resp]],
      age.pred = fit$predicted,
    ),
    # full CR 4 & 5 model to predict CR 2 & 3
    tibble(
      swfsc.id = test.df$swfsc.id,
      age.resp = test.df[[resp]],
      age.pred = pred,
    )) |>
    mutate(
      resid = age.pred - age.resp,
      dev = abs(resid),
      sq.err = resid ^ 2
    )
}
