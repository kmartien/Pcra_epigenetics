library(ggplot2)
library(dplyr)
library(gridExtra)
library(tidyverse)
load("data/color.palettes.rda")
load("data/model.params.rda")

description <- "logit.min.cov.100"

mae <- list()

# load glmnet MAE summaries
mae$glmnet <- read.csv("results-raw/glmnet.age.mae.summary.csv")

# load RF MAE summaries
mae$RF <- read.csv("results-raw/rf.age.mae.summary.csv")

# load SVM MAE summaries
mae$SVM <- read.csv("results-raw/svm.age.mae.summary.csv")

mae.df <- do.call(bind_rows, lapply(1:length(mae), function(x){
    bind_cols(method = names(mae)[x], select(mae[[x]], c(i, CR, MAE)))
}))
#mae.df$model <- as.factor(mae.df$i)
names(mae.df)[2] <- "model"
mae.df$CR <- fct_inorder(mae.df$CR)
mae.df <- left_join(mae.df, model.params)


plots <- 
  ggplot(filter(mae.df, CR %in% c("4","5","4&5")), aes(fill = CR, y = MAE, x = model, col = "black")) +
  geom_bar(colour = "black", position="dodge", stat = "identity") +
  labs(title = "Median Age Error for samples with CR = 4 or 5", x = "model",
       y = "Median age error (yrs)") +
  #    ggtitle(paste0("Median Age Error for samples with CR = ", cr)) +
  scale_fill_manual(values = training.set.palette, name = "CR", 
                    labels = c("CR \u2265 4", "CR \u2265 5", "CR \u2265 4&5")) +
  theme(axis.text.x = element_text(angle = 30, hjust=1)) +
  facet_wrap(~method, nrow = 3)
jpeg(file = "results-raw/method.and.wt.comparison2.jpg", width = 1480, height = 720)
plots
dev.off()

age.error <- list()
age.error$glmnet <- read.csv("results-raw/glmnet.age.loov.error.summary.csv")
age.error$rf <- read.csv("results-raw/rf.age.loov.error.summary.csv")
age.error$svm <- read.csv("results-raw/svm.age.loov.error.summary.csv")
age.error.df <- do.call(bind_rows, lapply(1:length(age.error), function(x){
  bind_cols(method = names(age.error)[x], select(age.error[[x]], c(model, id, age.best, error, sex, age.confidence)))
})) %>%
  left_join(model.params)
age.error.df$abs.error <- abs(age.error.df$error)


all.sites <- filter(age.error.df, age.confidence > 3) %>% filter(site.select.regr.meth == "none")
all.sites$training.cr <- as.factor(all.sites$training.cr)
all.sites$wt.label <- paste0(all.sites$method, "-", all.sites$weight)
all.sites <- filter(all.sites, !wt.label %in% c("svm-CRC", "svm-linear"))
all.sites$training.cr.label <- paste0(all.sites$method, "-", all.sites$training.cr)

training.cr.labs <- c("ALL", "CR3+", "CR4+")
names(training.cr.labs) <- c("2", "3", "4")
plots <- 
  ggplot(all.sites) +
  geom_boxplot(aes(x = training.cr.label, y = abs.error)) +
  labs(title = "Models trained with all CpG sites", x = "Regression method and training sample set",
       y = "Absolute age error (yrs)") +  
  theme(text = element_text(size = 24), axis.text.x = element_text(angle = 30, hjust=1)) +
#  facet_wrap(~training.cr, nrow = 3, labeller = labeller(training.cr = training.cr.labs))  
  facet_wrap(~weight, nrow = 3)  
jpeg(file = "results-raw/age.error.boxplot.allCpGs.jpg", width = 720, height = 800)
plots
dev.off()


# plots <- 
#   ggplot(filter(mae.df, CR %in% c("4","5","4&5")), aes(fill = min.cr, y = MAE, x = m.wt, col = "black")) +
#   geom_bar(colour = "black", position="dodge", stat = "identity") +
#   labs(title = "Median Age Error for samples with CR = 4 or 5", x = "Regression method",
#        y = "Median age error (yrs)") +
#   #    ggtitle(paste0("Median Age Error for samples with CR = ", cr)) +
#   scale_fill_manual(values = training.set.palette, name = "Training samples:", 
#                     labels = c("CR \u2265 2", "CR \u2265 3", "CR \u2265 4")) +
#   theme(axis.text.x = element_text(angle = 30, hjust=1)) +
#   facet_wrap(~CR, nrow = 3, labeller = label_both)
# jpeg(file = "results-raw/method.and.wt.comparison.jpg", width = 480, height = 720)
# plots
# dev.off()
# 
# # same basic plot, but only showing CR = 4&5
# plot.4.5 <- 
#   ggplot(filter(mae.df, CR=="4&5"), aes(fill = min.cr, y = MAE, x = m.wt, col = "black")) +
#   geom_bar(colour = "black", position="dodge", stat = "identity") +
#   #    labs(title = "Median Age Error for samples with CR = 4 and 5", 
#   labs(x = "Regression method", y = "Median age error") +
#   scale_fill_manual(values = training.set.palette, name = "Training\nsamples:", 
#                     labels = c("CR \u2265 2", "CR \u2265 3", "CR \u2265 4")) +
#   theme(text = element_text(size = 20), axis.text.x = element_text(size = 15, angle = 30, hjust=1))
# jpeg(file = "results-raw/method.and.wt.comparison.4and5.only.jpg", width = 800, height = 480)
# plot.4.5
# dev.off()
# 
