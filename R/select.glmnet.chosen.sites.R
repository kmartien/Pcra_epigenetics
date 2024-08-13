select.glmnet.chosen.sites <- function(params, incl.prob.threshold){
  load(file = paste0("results/glmnet_stability_results-sites", params$site.select.params, "_logit_mincov100.rda"))
  filter(site.incl.counts, prob >= incl.prob.threshold)
}