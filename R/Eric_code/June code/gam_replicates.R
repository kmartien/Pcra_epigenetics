rm(list = ls())
library(tidyverse)
library(mgcv)

load("../age_and_methylation_data.rdata")
ran.df <- readRDS('../age_meth_sample.rds')
sites <- readRDS('../random forest/rf_imp_sites.rds')

ids.to.use <- age.df |> 
  filter(swfsc.id %in% ids.to.keep & age.confidence %in% 4:5) |> 
  pull('swfsc.id')

gam.rep <- lapply(ids.to.use, function(cv.id) {
  cat(format(Sys.time()), ': ', cv.id, ', ', match(cv.id, ids.to.use), '/', length(ids.to.use), '\n', sep = '')
  gam.fit <- parallel::mclapply(1:500, function(i) {
    model.df <- ran.df[[i]] |> 
      filter(swfsc.id %in% ids.to.use & swfsc.id != cv.id) |> 
      column_to_rownames('swfsc.id')
    tryCatch(
      gam(
        formula = as.formula(paste('ran.age ~', paste0('s(', sites, ', bs = "ts")', collapse = ' + '))),
        data = model.df,
        select = TRUE
      ), 
      error = function(e) NULL
    )
  }, mc.cores = 10)
  list(cv.id = cv.id, gam.fit = gam.fit)
})

save(gam.rep, file = format(Sys.time(), 'gam_replicates_%Y%m%d_%H%M.rdata'))

# if(Sys.info()['nodename'] == 'r-instance1') {
#   system('gsutil cp gam_replicates*.rdata gs://nmfs-swfscbeam-static/')
#   system('gcloud compute instances stop r-instance1 --zone us-west2-b')
# }