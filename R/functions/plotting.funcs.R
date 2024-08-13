library(ggplot2)
library(gridExtra)
library(viridis)

# Age distribution of training data set
p.age.distribution <- function(dat, min.CR){
  ggplot(filter(dat, age.confidence >= min.CR), aes(x = age.best, fill = as.character(age.confidence))) +
    geom_histogram(bins = 35, binwidth = 1, color = "black") +
    scale_fill_manual(values = conf.palette, name = "Confidence") +
    labs(x = "Age.best", y = "Count") +
    theme_minimal() +
    theme(
      text = element_text(size = 15),
      legend.position = c(0.6, 0.8)
    )
}

# Bar chart of social cluster membership by sex
p.cluster.distribution <- function(dat) {
  dat$cluster.louvain <- as.numeric(dat$cluster.louvain)
  ggplot(dat, aes(x = cluster.louvain, col = I("grey"), fill = sex)) +
    geom_histogram(binwidth = 1) +
    scale_fill_manual(values = sex.palette) +
    labs(x = "Social cluster", y = "Count") +
    theme_minimal() +
    theme(
      text = element_text(size = 15),
      legend.position = c(0.9, 0.8)
    )
}

# plot loov results by min.CR
plot.loov.res <- function(loov.res, min.CR) {# expected argument is loov.res
  dat <- filter(loov.res, age.confidence >= min.CR)
  med.error <- group_by(dat, age.confidence) %>% summarise(round(median(abs(error)),2))
  names(med.error) <- c("CR", "MAE")
  loov.regress <- do.call(rbind,lapply(min.CR:5, function(cr){
    dat.filtered <- filter(dat, age.confidence == cr)
    regr <- lm(predicted.age~age.best, data = dat.filtered)$coefficients
    names(regr) <- c("intercept", "slope")
    
    cor.coeff <- round(cor.test(dat.filtered$age.best, dat.filtered$predicted.age, method = "pearson")$estimate,2)
    return(c(cor.coeff, regr))
  }))
  colnames(loov.regress)[1] <- "Corr"
  fit.sum <- bind_cols(med.error, loov.regress)
  
  p.loov <- ggplot(dat) +
    geom_abline(slope = 1, color = "black", linewidth = 0.5, linetype = 2) +
    geom_point(aes(x = age.best, y = predicted.age, col = as.character(age.confidence)), size = 3) +
    annotation_custom(tableGrob(fit.sum[,1:3], theme = ttheme_minimal(base_size = 10), rows = NULL), xmin = 30, xmax = 40, ymin = 0, ymax = 10) +
    scale_color_manual(values = conf.palette, name = "Confidence") +
    ggtitle(paste0("LOOV min CR = ", min.CR)) +
    labs(x = "Age.best", y = "LOOCV predicted age") +
    #  facet_wrap(~ sex, ncol = 1) +
    xlim(0,40) + ylim(0,52) +
    theme_minimal() +
    theme(
      text = element_text(size = 15),
      legend.position = c(0.2, 0.8)
    )
  for(i in 1:nrow(fit.sum)){
    cr <- i+(min.CR-1)
    p.loov <- p.loov + 
      geom_abline(slope = filter(fit.sum, CR == cr)$slope, intercept = filter(fit.sum, CR == cr)$intercept, color = conf.palette[cr])
  }
  return(list(p.loov = p.loov, fit.sum = fit.sum))
}

# Distribution of deviations for glmnet
plot.deviation <- function(dat, min.CR){
  dat$error <- dat$predicted.age - dat$age.best
  dat$age.range <- dat$age.max - dat$age.min
  dat$age.confidence <- as.character(dat$age.confidence)
  ggplot(filter(dat, age.confidence >= min.CR)) +
    geom_point(aes(x = error, y = age.range, col = age.confidence)) +
    scale_colour_manual(values = conf.palette) +
    labs(x = "Predicted age - CRC age.best", y = "Age.max - age.min") +
    ggtitle(paste0("glmnet Training (CR >=", min.CR, ")")) +
    theme_minimal() +
    theme(
      text = element_text(size = 15)
    )
}

# Histogram of age errors from LOOCV results
loov.hist <- function(dat, min.cr){ # expected argument is loov.res
  dat$error <- dat$predicted.age - dat$age.best
  ggplot(filter(dat, age.confidence >= min.cr), aes(x = error, fill = as.character(age.confidence))) +
    geom_histogram(bins = 35, binwidth = 1, color = "black") +
    scale_fill_manual(values = conf.palette, name = "Confidence") +
    labs(x = "Predicted age - Age.best", y = "Count") +
    theme_minimal() +
    theme(
      text = element_text(size = 15),
      legend.position = c(0.2, 0.8)
    )
}

