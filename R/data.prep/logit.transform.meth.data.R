#load("data/uncorrected.pct.meth.min.cov.100.Rdata")
#load("data/pct.conversion.rda")

logit.transform.meth.data <- function(meth.dat, pct.conversion){
  
  corrected.logit.meth <- do.call('cbind',lapply(1:length(meth.dat), function(amp){
    s <- meth.dat[[amp]]
    last.pos <- dim(s$CpG.sum$freq.meth)[2]
    meth <- s$CpG.sum$freq.meth[,2:last.pos]
    cov <- s$CpG.sum$coverage[,2:last.pos]
    err <- s$CpG.sum$errors[,2:last.pos]
    x <- meth/(cov-err)
    x <- x/pct.conversion
    rownames(x) <- s$CpG.sum$freq.meth[,1]
    colnames(x) <- paste(names(meth.dat)[amp],colnames(meth),sep=".")
    for (j in 1:ncol(x)){
      col.min.nonzero <- min(x[which(x[,j] > 0),j])
      x[which(x[,j] == 0),j] <- col.min.nonzero - (max(x[,j], na.rm = TRUE) - col.min.nonzero) * 0.01
    }
    return(log(x/(1-x)))
  }))
}
