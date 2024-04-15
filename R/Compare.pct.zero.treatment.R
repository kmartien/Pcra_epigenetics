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
load(paste0("results/glmnet.absolute.age.", label, ".pct.results.rda"))
age.sum.combined <- age.sum
load(paste0("results/glmnet.absolute.age.", label, ".pct.no.zero.results.rda"))
age.sum.combined <- left_join(age.sum.combined, select(age.sum, c(id, predicted.age)), by = "id")
names(age.sum.combined)[c(14,15)] <- c("pred.age.zeros", "pred.age.no.zeros")
age.sum.combined$pct.diff <- (age.sum.combined$pred.age.no.zeros - age.sum.combined$pred.age.zeros)/ age.sum.combined$age.point

conf.palette <- colorRampPalette(viridis::magma(30, direction = -1)[-(1:5)])(5)[5:1]
conf.palette[1] <- 'gray40'
names(conf.palette) <- c(1:5)

cor.coeff.zeros <- round(cor.test(age.sum.combined$age.point, age.sum.combined$pred.age.zeros, method = "pearson")$estimate,2)
cor.coeff.no.zeros <- round(cor.test(age.sum.combined$age.point, age.sum.combined$pred.age.no.zeros, method = "pearson")$estimate,2)

g.pct.diff <- ggplot(filter(age.sum.combined, type == "Testing")) +
  geom_point(aes(x = age.point, y = pct.diff, shape = type, col = as.character(confidence)), size = 3) +
  scale_color_manual(values = conf.palette, name = "Confidence") +
  ggtitle("GLMNET prediction with and without zeros replaced") +
  labs(x = "Age.best", y = "Proportion age increase when zeros replaced") +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = c(0.8, 0.8)
  )
g.compare.zeros <- ggplot(filter(age.sum.combined, type == "Testing")) +
  geom_abline(slope = 1, color = "black", linewidth = 0.5, linetype = 2) +
  geom_point(aes(x = age.point, y = pred.age.zeros, col = as.character(confidence)), size = 3) +
  scale_color_manual(values = conf.palette, name = "Confidence") +
  ggtitle(paste0("GLMNET prediction with zeros; correlation = ", cor.coeff.zeros)) +
  labs(x = "Age.best", y = "Pred age with zeros") +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = c(0.9, 0.2)
  )
g.compare.no.zeros <- ggplot(filter(age.sum.combined, type == "Testing")) +
  geom_abline(slope = 1, color = "black", linewidth = 0.5, linetype = 2) +
  geom_point(aes(x = age.point, y = pred.age.no.zeros, col = as.character(confidence)), size = 3) +
  scale_color_manual(values = conf.palette, name = "Confidence") +
  ggtitle(paste0("GLMNET prediction zeros replaced; correlation = ", cor.coeff.no.zeros)) +
  labs(x = "Age.best", y = "Pred age zeros replaced") +
  theme_minimal() +
  theme(
    text = element_text(size = 15),
    legend.position = c(0.9, 0.2)
  )
#jpeg(file = "results/compare.glmnet.zero.treatment.jpg", width = 1000, height = 2500)
pdf(file = "results/compare.glmnet.zero.treatment.pdf")
#grid.arrange(g.pct.diff, g.compare.zeros, g.compare.no.zeros, nrow = 3)
g.pct.diff
g.compare.zeros
g.compare.no.zeros
dev.off()
