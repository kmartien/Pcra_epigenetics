library(dplyr)
library(ggplot2)
library(gridExtra)

training.min.CR <- 3
load(paste0("results/glmnet.absolute.age.minCR", training.min.CR, ".results.rda"))
load(paste0("results/glmnet.LOOV.minCR", training.min.CR, ".results.rda"))
load("data/age.data.rda")

# Age distribution of training data set
p.age.distribution <- ggplot(age.sum, aes(x = age.point, fill = sex)) +
  geom_histogram(binwidth = 1, color = "black") +
  labs(x = "Age.best", y = "Count") +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = c(0.6, 0.8)
  )

jpeg(file = paste0("results/age.distribution.histogram.minCR", training.min.CR, ".jpg"))
p.age.distribution
dev.off()

# Regression plots of glmnet cv.fit results and LOOCV results
cv.fit.regression <- lm(predicted.age~age.point, data = age.sum)
R2.cv.git <- data.frame(txt = paste0("R-squared: ", round(summary(cv.fit.regression)$adj.r.squared, 2)), 
                 xpos = 40, ypos = 0, h = 1, v = 0)
loov.regression <- lm(predicted.age~age.point, data = loov.res)
R2.loov <- data.frame(txt = paste0("R-squared: ", round(summary(loov.regression)$adj.r.squared, 2)), 
                      xpos = 40, ypos = 0, h = 1, v = 0)

R2.loov
mean(abs(loov.res$error))
median(abs(loov.res$error))

p.cv.fit <- ggplot(age.sum) +
  geom_abline(slope = 1, color = "black", linewidth = 0.5, linetype = 2) +
  geom_abline(slope = cv.fit.regression$coefficients[2], intercept = cv.fit.regression$coefficients[1]) +
  geom_point(aes(x = age.point, y = predicted.age, col = sex), size = 3) +
  geom_text(data = R2.cv.git, aes(x = xpos, y = ypos, hjust = h, vjust = v, label = txt)) +
  #  scale_fill_discrete(c("red","green")) +
  labs(x = "Age.best", y = "Predicted age") +
  #  facet_wrap(~ sex, ncol = 1) +
  xlim(0,40) + ylim(0,40) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = c(0.2, 0.8)
  )

p.loov <- ggplot(loov.res) +
  geom_abline(slope = 1, color = "black", linewidth = 0.5, linetype = 2) +
  geom_abline(slope = loov.regression$coefficients[2], intercept = loov.regression$coefficients[1]) +
  geom_point(aes(x = age.point, y = predicted.age, col = sex), size = 3) +
  geom_text(data = R2.loov, aes(x = xpos, y = ypos, hjust = h, vjust = v, label = txt)) +
  #  scale_fill_discrete(c("red","green")) +
  labs(x = "Age.best", y = "LOOCV predicted age") +
  #  facet_wrap(~ sex, ncol = 1) +
  xlim(0,40) + ylim(0,40) +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = "none"
  )

jpeg(file = paste0("results/glmnet.regression.plots.minCR", training.min.CR, ".jpg"), width = 900, height = 450)
plots <- list(p.cv.fit, p.loov)
plots$nrow <- 1
do.call(grid.arrange,plots)
dev.off()

# Histogram of age errors from LOOCV results
loov.hist <- ggplot(loov.res, aes(x = error, fill = sex)) +
  geom_histogram(bins = 35, binwidth = 1, color = "black") +
  labs(x = "Predicted age - Age.best", y = "Count") +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = c(0.2, 0.8)
  )
  
jpeg(file = paste0("results/loov.error.histogram.minCR", training.min.CR, ".jpg"))
loov.hist
dev.off()

