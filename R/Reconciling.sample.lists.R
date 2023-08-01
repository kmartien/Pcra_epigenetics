library(dplyr)

load("results/From.Eric/age_and_methylation_data.rdata")

#median coverage per individual without na.rm = TRUE
med.cov.ind <- sapply(1:nrow(meth.cov), function(i){
  median(meth.cov[i,])
})
print(rownames(meth.cov)[which(med.cov.ind < 1000)])

meth.cov.no.na.rm <- meth.cov[-which(med.cov.ind < 1000),]

missing.data.ind <- sapply(1:nrow(meth.cov.no.na.rm), function(i){
  length(which(is.na(meth.cov.no.na.rm[i,])))
})
print(rownames(meth.cov)[missing.data.ind >0])

#median coverage per individual with na.rm = TRUE
med.cov.ind <- sapply(1:nrow(meth.cov), function(i){
  median(meth.cov[i,], na.rm = TRUE)
})
print(rownames(meth.cov)[which(med.cov.ind < 1000)])

meth.cov <- meth.cov[-which(med.cov.ind < 1000),]

missing.data.ind <- sapply(1:nrow(meth.cov), function(i){
  length(which(is.na(meth.cov[i,])))
})
print(rownames(meth.cov)[missing.data.ind >0])
