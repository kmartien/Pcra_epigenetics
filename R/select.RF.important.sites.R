library(rfPermute)
select.RF.important.sites <- function(params, incl.prob.threshold){

  load(file = paste0("results/rf_stability_results-sites", params$site.select.params, "_logit_mincov100.rda"))
  temp <- filter(importance.sum, pval <= 0.05) %>% group_by(site) %>% summarise(prob = n()) %>% filter(prob > incl.prob.threshold)
  rm(imp.pvals, importance.sum, rf.list)
  return(temp)
}