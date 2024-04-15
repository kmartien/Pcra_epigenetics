# imp should have two columns:
#  "site" is the name of the CpG site in the form locus.position
#  the second column is the measure of importance (e.g., inclusion, %IncMSE, etc.)
#      the name of the second column does not matter
plot.site.importance <- function(imp){
  
  imp <- bind_cols(do.call(rbind, strsplit(imp$site, split = "[.]")), imp)
  names(imp) <- c("locus", "pos", "site", "importance")
  imp$pos <- as.numeric(imp$pos)

  ggplot(imp, aes(x = pos, y = importance)) +
    geom_point() +
    facet_wrap(~ locus)
}
