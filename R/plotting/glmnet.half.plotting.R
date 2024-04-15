library(dplyr)
library(ggplot2)
library(gridExtra)
library(viridis)

label <- "model.6"
load("data/age.data.rda")
load("data/bayesian.models.smry.rdata")

sex.palette <- viridis(n = 5)[c(2,4)]
conf.palette <- colorRampPalette(viridis::magma(30, direction = -1)[-(1:5)])(5)[5:1]
conf.palette[1] <- 'gray40'
names(conf.palette) <- c(1:5)
# mako(n = 5)[c(2,4)] also looks good, but the blue might be too similar to CR = 5 color

# Regression plot of glmnet Testing results
load(paste0("results/glmnet.absolute.age.", label, ".results.rda"))
age.sum$error <- abs(age.sum$predicted.age - age.sum$age.point)
test.dat <- filter(age.sum, type == "Testing")

med.error <- group_by(test.dat, confidence) %>% summarise(round(median(error),2))
names(med.error) <- c("CR", "MAE")
glmnet.regress <- do.call(rbind,lapply(2:5, function(cr){
  dat <- filter(test.dat, confidence == cr)
  regr <- lm(predicted.age~age.point, data = dat)$coefficients
  names(regr) <- c("intercept", "slope")
  
  cor.coeff <- round(cor.test(dat$age.point, dat$predicted.age, method = "pearson")$estimate,2)
  return(c(cor.coeff, regr))
}))
colnames(glmnet.regress)[1] <- "Corr"
fit.sum <- bind_cols(med.error, glmnet.regress)

p.cv.testing <- ggplot(test.dat) +
  geom_abline(slope = 1, color = "black", linewidth = 0.5, linetype = 2) +
  geom_point(aes(x = age.point, y = predicted.age, col = as.character(confidence)), size = 3) +
  annotation_custom(tableGrob(fit.sum[,1:3], theme = ttheme_minimal(), rows = NULL), xmin = 30, xmax = 40, ymin = 0, ymax = 10) +
  scale_color_manual(values = conf.palette, name = "Confidence") +
  ggtitle(paste0("Testing data half of samples CR = ", label)) +
  labs(x = "Age.best", y = "Predicted age") +
  xlim(0,40) + ylim(0,43) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = c(0.2, 0.8)
  )
for(i in 1:nrow(fit.sum)){
  cr <- i+1
  p.cv.testing <- p.cv.testing + 
    geom_abline(slope = filter(fit.sum, CR == cr)$slope, intercept = filter(fit.sum, CR == cr)$intercept, color = conf.palette[cr])
}

# Regression plot of Bayes Testing results
dat <- model.smry[[which(names(model.smry) == label)]]
names(dat)[1] <- "swfsc.labid"
test.dat <- left_join(dat, age.data) %>% mutate(error = abs(pred.mode - age.best)) %>%
  filter(type == "Testing")

med.error <- group_by(test.dat, age.confidence) %>% summarise(round(median(error),2))
names(med.error) <- c("CR", "MAE")
bayes.regress <- do.call(rbind,lapply(2:5, function(cr){
  dat <- filter(test.dat, age.confidence == cr)
  regr <- lm(pred.mode~age.best, data = dat)$coefficients
  names(regr) <- c("intercept", "slope")
  
  cor.coeff <- round(cor.test(dat$age.best, dat$pred.mode, method = "pearson")$estimate,2)
  return(c(cor.coeff, regr))
}))
colnames(bayes.regress)[1] <- "Corr"
fit.sum <- bind_cols(med.error, bayes.regress)

p.bayes <- ggplot(test.dat) +
  geom_abline(slope = 1, color = "black", linewidth = 0.5, linetype = 2) +
  geom_point(aes(x = age.best, y = pred.mode, col = as.character(age.confidence)), size = 3) +
  scale_color_manual(values = conf.palette, name = "Confidence") +
  annotation_custom(tableGrob(fit.sum[,1:3], theme = ttheme_minimal(), rows = NULL), xmin = 30, xmax = 40, ymin = 0, ymax = 10) +
  ggtitle(paste0("Model 3 Testing samples CR = ", label)) +
  #  scale_fill_discrete(c("red","green")) +
  labs(x = "Age.best", y = "Predicted age") +
  #  facet_wrap(~ sex, ncol = 1) +
  xlim(0,40) + ylim(0,43) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = c(0.2, 0.8)
  )
for(i in 1:nrow(fit.sum)){
  cr <- i+1
  p.bayes <- p.bayes + 
    geom_abline(slope = filter(fit.sum, CR == cr)$slope, intercept = filter(fit.sum, CR == cr)$intercept, color = conf.palette[cr])
}

#plots <- unlist(plots, recursive = FALSE)
plots <- list(p.cv.testing, p.bayes)
plots$nrow <- 1
jpeg(file = paste0("results/glmnet.regression.plots.half.",label,".jpg"), width = 1000, height = 500)
do.call(grid.arrange,plots)
dev.off()
