library(dplyr)
library(tidyr)
library(swfscMisc)

load("data/meth.dat.Rdata")
meth.samps <- meth.dat[[1]]$CpG.sum$coverage$id
full.CRC.data <- read.csv("data-raw/CRC.age.estimates.2022SEP.csv")

names(full.CRC.data)[c(1,2,4,5,6,7,21,22,25,27,28,33)] <- 
  c("crc.id",
    "social.cluster",
    "sex",
    "swfsc.labid",
    "biopsy.id",
    "date.biopsy",
    "age.best",
    "age.confidence",
    "age.min",
    "age.max.CRR",
    "ADULT.max.best",
    "pair.id")

# Eric wants:
#  [1] crc.id (first column in sheet)
#  [5] swfsc.labid (so it matches the id in the methylation data, eg., z00#####)
#  [42] age.class (Calf, Juvenile, Subadult, Adult) - this should only be if the age class (at time of biopsy) is certain
#  [4] sex (Male, Female)
#  [6] biopsy.id (to match with the field ID in the minimum tooth age spreadsheet from Kelly)
#  [7] date.biopsy (YY-MM-DD)
#  [23] age.best (CRC best estimate of age at biopsy)
#  [24] age.confidence (confidence rating 1:5)
#  [27] age.min (minimum age at biopsy)
#  [41] age.max (best maximum age at biopsy)
#  [35] pair.id (number grouping same individual biopsied multiple times)


full.CRC.data$age.max <- sapply(1:nrow(full.CRC.data), function(i){
  if(is.na(full.CRC.data$ADULT.max.best[i]) && is.na(full.CRC.data$age.max.CRR[i])) {
    return(NA)} else{
  min(full.CRC.data$ADULT.max.best[i], full.CRC.data$age.max.CRR[i], na.rm = TRUE) }
})

full.CRC.data$swfsc.labid <- as.numeric(full.CRC.data$swfsc.labid)
full.CRC.data$swfsc.labid <- paste0("z0",zero.pad(full.CRC.data$swfsc.labid))

full.CRC.data$age.class <- "adult"
full.CRC.data$age.class[which(full.CRC.data$sex=="Female" & full.CRC.data$age.best < 10)] <- "sub-adult"
full.CRC.data$age.class[which(full.CRC.data$sex=="Female" & full.CRC.data$age.best < 6)] <- "juvenile"
full.CRC.data$age.class[which(full.CRC.data$sex=="Female" & full.CRC.data$age.best < 3)] <- "calf"

full.CRC.data$age.class[which(full.CRC.data$sex=="Male")] <- "adult-PM"
full.CRC.data$age.class[which(full.CRC.data$sex=="Male" & full.CRC.data$age.best < 25)] <- "adult-SM"
full.CRC.data$age.class[which(full.CRC.data$sex=="Male" & full.CRC.data$age.best < 15)] <- "sub-adult"
full.CRC.data$age.class[which(full.CRC.data$sex=="Male" & full.CRC.data$age.best < 9)] <- "juvenile"
full.CRC.data$age.class[which(full.CRC.data$sex=="Male" & full.CRC.data$age.best < 3)] <- "calf"

age.data <- full.CRC.data[which(full.CRC.data$swfsc.labid %in% meth.samps),c(1,2,5,40,4,6,7,21,22,25,39,33)]

save(age.data,file = "data/age.data.rda")
