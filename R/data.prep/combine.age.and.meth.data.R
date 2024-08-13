combine.age.and.meth.data <- function(description){
  load(paste0("data/corrected.meth_",description,".Rdata"))
  load("data/age.data.rda")
  load("data/sites.and.inds.from.Eric.rda")
  wts <- data.frame(age.confidence = 1:5, wt = c(0.05, 0.2, 0.55, 0.75, 1))
  age.data <- left_join(age.data, wts)
  
  age.data <- filter(age.data, swfsc.labid %in% ids.to.keep)
  sites.2.keep <- gsub("_00", replacement = ".", sites.2.keep)
  sites.2.keep <- gsub("_0", replacement = ".", sites.2.keep)
  sites.2.keep <- gsub("_", replacement = ".", sites.2.keep)
  corrected.meth <- select(corrected.meth, sites.2.keep)
  
  site.names <- names(corrected.meth)
  corrected.meth <- cbind(id=rownames(corrected.meth),corrected.meth)
  names(age.data)[which(names(age.data) %in% c("swfsc.labid"))] <- c("id")
  age.data$numeric.sex <- 0
  age.data$numeric.sex[which(age.data$sex=="Female")] <- 1
  age.data$numeric.sex[which(age.data$sex=="Male")] <- 0
  age.data$decade <- floor(age.data$age.best/10)
  age.data$decade[which(age.data$decade>3)] <- 3
  first.meth.col <- dim(age.data)[2] + 1
  
  all.samples <- left_join(age.data,corrected.meth)
  return(list(site.names = site.names, first.meth.col = first.meth.col, all.samples = all.samples))
}
