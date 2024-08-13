library(viridis)

sex.palette <- viridis(n = 5)[c(2,4)]
conf.palette <- colorRampPalette(viridis::magma(30, direction = -1)[-(1:5)])(5)[5:1]
conf.palette[1] <- 'gray40'
names(conf.palette) <- c(1:5)
#training.set.palette <- mako(n = 5)[c(1,3,5)] # also looks good, but the blue might be too similar to CR = 5 color
training.set.palette <- c("black","grey65", "grey95")

save(sex.palette, conf.palette, training.set.palette, file = "data/color.palettes.rda")
