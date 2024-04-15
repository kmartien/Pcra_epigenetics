library(dplyr)
library(tidyr)
library(swfscMisc)

load("data/meth.dat.Rdata")
meth.samps <- meth.dat[[1]]$CpG.sum$coverage$id
full.CRC.data <- read.csv("data-raw/CRC.age.estimates.2024FEB.csv")

names(full.CRC.data)[c(1,2,4,5,6,7,25,26,29,34,38)] <- 
  c("crc.id",
    "social.cluster",
    "sex",
    "swfsc.labid",
    "biopsy.id",
    "date.biopsy",
    "age.best",
    "age.confidence",
    "age.min",
    "age.max",
    "pair.id")

# Eric wants:
#  [1] crc.id (first column in sheet)
#  [5] swfsc.labid (so it matches the id in the methylation data, eg., z00#####)
#  [44] age.class (Calf, Juvenile, Subadult, Adult) - this should only be if the age class (at time of biopsy) is certain
#  [4] sex (Male, Female)
#  [6] biopsy.id (to match with the field ID in the minimum tooth age spreadsheet from Kelly)
#  [7] date.biopsy (YY-MM-DD)
#  [25] age.best (CRC best estimate of age at biopsy)
#  [26] age.confidence (confidence rating 1:5)
#  [29] age.min (minimum age at biopsy)
#  [34] age.max (best maximum age at biopsy)
#  [38] pair.id (number grouping same individual biopsied multiple times)


full.CRC.data$swfsc.labid <- as.numeric(full.CRC.data$swfsc.labid)
full.CRC.data <- full.CRC.data[-which(is.na(full.CRC.data$swfsc.labid)),]
full.CRC.data$swfsc.labid <- paste0("z0",zero.pad(full.CRC.data$swfsc.labid))

full.CRC.data$age.class <- "adult"
full.CRC.data$age.class[which(full.CRC.data$sex=="Female" & full.CRC.data$age.best < 10)] <- "sub-adult"
full.CRC.data$age.class[which(full.CRC.data$sex=="Female" & full.CRC.data$age.best < 6)] <- "juvenile"
full.CRC.data$age.class[which(full.CRC.data$sex=="Female" & full.CRC.data$age.best < 3)] <- "calf"

full.CRC.data$age.class[which(full.CRC.data$sex=="Male")] <- "adult-PM"
full.CRC.data$age.class[which(full.CRC.data$sex=="Male" & full.CRC.data$age.best < 25)] <- "adult-DPM"
full.CRC.data$age.class[which(full.CRC.data$sex=="Male" & full.CRC.data$age.best < 15)] <- "adult"
full.CRC.data$age.class[which(full.CRC.data$sex=="Male" & full.CRC.data$age.best < 10)] <- "sub-adult"
full.CRC.data$age.class[which(full.CRC.data$sex=="Male" & full.CRC.data$age.best < 6)] <- "juvenile"
full.CRC.data$age.class[which(full.CRC.data$sex=="Male" & full.CRC.data$age.best < 3)] <- "calf"

age.data <- full.CRC.data[which(full.CRC.data$swfsc.labid %in% meth.samps),c(1,2,5,44,4,6,7,25,26,29,34,38)]

save(age.data,file = "data/age.data.rda")
