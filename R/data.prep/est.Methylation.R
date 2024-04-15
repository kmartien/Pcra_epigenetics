library(dplyr)
library(tidyr)
library(readr)
library(seqinr)
library(swfscGenetics)
source("R/calc_methylation.R")
source("R/logit.transform.meth.data.R")

min.cov <- 100
description <- paste0("min.cov.",min.cov)

bs.refs <- read.fasta("data-raw/Pcra.bisulfite.amplicon.wPrimers.fasta")
orig.refs <- read.fasta("data-raw/Pcra.unconverted.amplicon.wPrimer.fasta")
primers <- read.fasta("data-raw/Pcra.epi.primers.fasta")
amplicons <- list("TET2","FLJ0945","DIRAS3","VGF","HMG20B","KCNC4","PRDM12","GRIA2")

#determine the locations within the references that aren't part of the primers
target.positions <- sapply(names(bs.refs), function(p){
  target.start <- length(primers[[which(names(primers)==paste(p,"_F1", sep=""))]]) +1
  target.end <- length(bs.refs[[which(names(bs.refs)==p)]]) - length(primers[[which(names(primers)==paste(p,"_R1", sep=""))]])
  return(c(start=target.start, end=target.end))
})

pileup.path <- "data-raw/MS15.galaxy.alignments/Aligned.to.targets.w.primers/pileups"
fnames <- list.files(path = pileup.path, 
                     pattern="\\.fastq.tabular.csv$")
num.samps <- length(fnames)

plp.all <- do.call("rbind",lapply(fnames, function(f){
  read.csv(paste0(pileup.path,"/",f)) %>% as.data.frame
}))

id <- do.call("c",lapply(plp.all$fname, function(f){
  strsplit(f,"_",fixed=TRUE)[[1]][1]
}))
plp.all$fname <- id
names(plp.all)[c(1,3)] <- c("id","position")
#plp.all <- plp.all[-which(is.na(plp.all$A)),] #this takes care of the error in pileupCallRun with insertions

meth.dat <- lapply(amplicons, function(amp){
  sub.plp <- subset(plp.all,chrom==amp)
  bs.ref.seq <- bs.refs[[which(names(bs.refs)==amp)]]
  orig.ref.seq <- orig.refs[[which(names(orig.refs)==amp)]]
  target.pos <- target.positions[,which(colnames(target.positions)==amp)]
  calc_methylation(sub.plp,bs.ref.seq, orig.ref.seq, target.pos, min.cov)
})
names(meth.dat) <- amplicons

save(meth.dat, file = paste0("data/uncorrected.pct.meth.",description,".Rdata"))

mean.cov <- do.call("cbind",lapply(meth.dat,function(s){s$CpG.sum$coverage$avg}))
rownames(mean.cov) <- meth.dat$TET2$CpG.sum$coverage$id
colnames(mean.cov) <- amplicons
write.csv(mean.cov, "results/mean.coverage.Aligned.to.targets.with.primers.csv")

non.CpG.coverage <- colSums(do.call('rbind', lapply(meth.dat, function(s){
  last.pos <- dim(s$non.CpG.sum$freq.meth)[2]
  if (!is.null(last.pos)){
    tot.cov <- rowSums(s$non.CpG.sum$coverage[,2:last.pos],na.rm = TRUE) - rowSums(s$non.CpG.sum$errors[,2:last.pos],na.rm = TRUE)
  }
})))

non.CpG.Ts <- colSums(do.call('rbind',lapply(meth.dat, function(s){
  last.pos <- dim(s$non.CpG.sum$freq.meth)[2]
  if (!is.null(last.pos)){
#    tot.cov <- rowSums(s$non.CpG.sum$coverage[,2:last.pos],na.rm = TRUE)
    tot.Ts <- (rowSums(s$non.CpG.sum$coverage[,2:last.pos],na.rm = TRUE)-rowSums(s$non.CpG.sum$freq.meth[,2:last.pos],na.rm = TRUE)-rowSums(s$non.CpG.sum$errors[,2:last.pos],na.rm = TRUE))
  }
})))

pct.conversion <- non.CpG.Ts/non.CpG.coverage
save(pct.conversion, file = "data/pct.conversion.rda")

corrected.pct.meth <- do.call('cbind',lapply(1:length(meth.dat), function(amp){
  s <- meth.dat[[amp]]
  last.pos <- dim(s$CpG.sum$freq.meth)[2]
  x <- (s$CpG.sum$freq.meth[,2:last.pos]/(s$CpG.sum$coverage[,2:last.pos]-s$CpG.sum$errors[,2:last.pos]))
  x <- x/pct.conversion
  rownames(x) <- s$CpG.sum$freq.meth[,1]
  colnames(x) <- paste(amplicons[amp],colnames(s$CpG.sum$freq.meth)[2:last.pos],sep=".")
  return(x)
}))
#write.csv(corrected.pct.meth, file=paste("Percent.methylated.",description,".csv",sep=""))
save(corrected.pct.meth,file=paste0("data/corrected.meth.pct.",description,".Rdata"))

corrected.logit.meth <- logit.transform.meth.data(meth.dat, pct.conversion)
save(corrected.logit.meth, file = paste0("data/corrected.meth.logit.",description,".Rdata"))
