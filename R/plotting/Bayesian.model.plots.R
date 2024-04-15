library(dplyr)
library(ggplot2)
library(gridExtra)
library(viridis)

training.min.CR <- 4
load(paste0("results/glmnet.absolute.age.minCR", training.min.CR, ".results.rda"))
load(paste0("results/glmnet.LOOV.minCR", training.min.CR, ".results.rda"))
load("results/glmnet.sex.cluster.rda")
load("data/age.data.rda")
load("data/bayesian.models.smry.rdata")

sex.palette <- viridis(n = 5)[c(2,4)]
conf.palette <- colorRampPalette(viridis::magma(30, direction = -1)[-(1:5)])(5)[5:1]
conf.palette[1] <- 'gray40'
names(conf.palette) <- c(1:5)
# mako(n = 5)[c(2,4)] also looks good, but the blue might be too similar to CR = 5 color

age.data$age.confidence <- as.character(age.data$age.confidence)

plots <- lapply(1:4, function(m){
  mdl <- model.smry[[m]]
  names(mdl)[1] <- "swfsc.labid"
  mdl <- left_join(mdl, age.data) %>% mutate(error = pred.mode - age.best) %>%
    filter(age.confidence %in% c(4,5))
  p <- lapply(c("Training", "Testing"), function(t){
    dat <- filter(mdl, type == t)
    if (nrow(dat) == 0) return(NULL) else {
    regress <- lm(pred.mode~age.best, data = dat)
    cor.coeff <- cor.test(dat$age.best, dat$pred.mode, method = "pearson")$estimate
    cor.coeff.label <- data.frame(txt = paste0("Corr. coeff.: ", round(cor.coeff, 2)), 
                     xpos = 40, ypos = 0, h = 1, v = 0)
    regression.plot <- ggplot(dat) +
      geom_abline(slope = 1, color = "black", linewidth = 0.5, linetype = 2) +
      geom_abline(slope = regress$coefficients[2], intercept = regress$coefficients[1]) +
      geom_point(aes(x = age.best, y = pred.mode, col = age.confidence), size = 3) +
      scale_color_manual(values = conf.palette) +
      geom_text(data = cor.coeff.label, aes(x = xpos, y = ypos, hjust = h, vjust = v, label = txt)) +
      labs(x = "Age.best", y = "Predicted age") +
      ggtitle(paste0(names(model.smry)[m], " ", t)) +
      #  facet_wrap(~ sex, ncol = 1) +
      xlim(0,40) + ylim(0,40) +
      theme_minimal() +
      theme(
        text = element_text(size = 15),
        legend.position = c(0.2, 0.8)
      )
    
    error.hist <- ggplot(dat, aes(x = error, fill = age.confidence)) +
      geom_histogram(bins = 35, binwidth = 1, color = "black") +
      scale_fill_manual(values = conf.palette) +
      labs(x = "Predicted age - Age.best", y = "Count") +
#      xlim(-30, 30) +
      xlim(-8, 8) +
      ylim(0, 15) +
      ggtitle(paste0(names(model.smry)[m], " ", t)) +
      theme_minimal() +
      theme(
        text = element_text(size = 15),
        legend.position = c(0.2, 0.8)
      )
    list(regr.plot = regression.plot, error.hist = error.hist)
    }
  })
  names(p) <- c("Training", "Testing")
  return(p)
})
names(plots) <- names(model.smry)

regr.plots <- unlist(lapply(plots, function(p){list(p[[1]]$regr.plot, p[[2]]$regr.plot)}), recursive = FALSE)
regr.plots <- regr.plots[-which(sapply(regr.plots, is.null))]
hist.plots <- unlist(lapply(plots, function(p){list(p[[1]]$error.hist, p[[2]]$error.hist)}), recursive = FALSE)
hist.plots <- hist.plots[-which(sapply(hist.plots, is.null))]

regr.plots$nrow <- 3
jpeg(file = "results/Bayes.regression.plots.CR4-5.jpg", width = 600, height = 900)
do.call(grid.arrange,regr.plots)
dev.off()

hist.plots$nrow <- 3
jpeg(file = "results/Bayes.error.histograms.CR4-5.jpg", width = 600, height = 900)
do.call(grid.arrange,hist.plots)
dev.off()

jpeg(file = "results/model1b.error.histogram.jpg")
hist.plots[[1]]
dev.off()

# Regression plot of model 1b testing data (CR < 4)
mdl <- model.smry[[1]]
names(mdl)[1] <- "swfsc.labid"
dat <- left_join(mdl, age.data) %>% mutate(error = pred.mode - age.best) %>%
  filter(type == "Testing")
regress <- lm(pred.mode~age.best, data = dat)
R2 <- data.frame(txt = paste0("R-squared: ", round(summary(regress)$adj.r.squared, 2)), 
                 xpos = 40, ypos = 0, h = 1, v = 0)
regression.plot <- ggplot(dat) +
  geom_abline(slope = 1, color = "black", linewidth = 0.5, linetype = 2) +
  geom_abline(slope = regress$coefficients[2], intercept = regress$coefficients[1]) +
  geom_point(aes(x = age.best, y = pred.mode, col = sex), size = 3) +
  scale_color_manual(values = conf.palette) +
  geom_text(data = R2, aes(x = xpos, y = ypos, hjust = h, vjust = v, label = txt)) +
  labs(x = "Age.best", y = "Predicted age") +
  ggtitle(paste0(names(model.smry)[1], " Testing data (CR < 4)")) +
  #  facet_wrap(~ sex, ncol = 1) +
  xlim(0,40) + ylim(0,40) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = c(0.2, 0.8)
  )

jpeg(file = "results/Bayes.error.regression.model1bTesting.jpg", width = 600, height = 300)
regression.plot
dev.off()

# Distribution of deviations for model 1b testing data (CR < 4)
mdl <- model.smry[[1]]
names(mdl)[1] <- "swfsc.labid"
dat <- left_join(mdl, age.data) %>% mutate(error = pred.mode - age.best) %>%
  filter(type == "Training")
dat$age.range <- dat$age.max - dat$age.min
deviation.plot <- ggplot(dat) +
  geom_point(aes(x = error, y = age.range, col = age.confidence)) +
  scale_colour_manual(values = conf.palette) +
  xlim(-6,6) +
  labs(x = "Modal predicted age - CRC age.best", y = "Age.max - age.min") +
  ggtitle("Model 1b Training (CR >=4)") +
  theme_minimal() +
  theme(
    text = element_text(size = 15)
  )

jpeg(file = "results/Deviation.plot.bayes.jpg", width = 800, height = 600)
deviation.plot
dev.off()
