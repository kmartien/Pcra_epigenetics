
# nrep <- 1000
# nfolds <- 10
# alpha.vec <- seq(0.0005, 0.25, by = 0.0005)
# 
# # cv.glmnet replicates with different values of alpha 
# glmnet.tune <- parallel::mclapply(1:nrep, function(i) {
#   foldid <- sample(1:nfolds, size = nrow(model.ages), replace = TRUE)
#   lapply(alpha.vec, function(alpha) {
#     cv.fit <- cv.glmnet(
#       x = as.matrix(model.ages[, sites.to.keep]),
#       y = model.ages$age.best,
#       weights = model.ages$confidence.wt,
#       alpha = alpha,
#       nfolds = nfolds,
#       foldid = foldid
#     ) 
#     list(
#       loc.site = cv.fit |> 
#         coef() |> 
#         as.matrix() |> 
#         as.data.frame() |> 
#         rownames_to_column('loc.site') |> 
#         filter(s1 > 0 & loc.site != '(Intercept)') |> 
#         pull(loc.site),
#       smry = data.frame(
#         alpha = alpha,
#         cvm.1se = cv.fit$cvm[cv.fit$lambda == cv.fit$lambda.1se], 
#         lambda.1se = cv.fit$lambda.1se,
#         cvm.min = cv.fit$cvm[cv.fit$lambda == cv.fit$lambda.min], 
#         lambda.min = cv.fit$lambda.min
#       )
#     )
#   })
# }, mc.cores = 5)
# 
# # summarize cvm and lambda across replicates
# glmnet.tune.smry <- glmnet.tune |> 
#   lapply(function(tune.rep) {
#     tune.rep |> 
#       lapply(function(alpha.rep) alpha.rep$smry) |> 
#       bind_rows()
#   }) |> 
#   bind_rows() |> 
#   group_by(alpha) |> 
#   summarize(
#     median.cvm = median(cvm.1se),
#     lower.cvm = unname(quantile(cvm.1se, 0.025)),
#     upper.cvm = unname(quantile(cvm.1se, 0.975)),
#     median.lambda.1se = median(lambda.1se),
#     lower.lambda.1se = unname(quantile(lambda.1se, 0.025)),
#     upper.lambda.1se = unname(quantile(lambda.1se, 0.975)),
#     .groups = 'drop'
#   )
# 
# # plot alpha vs. median lambda
# min.median.cvm <- glmnet.tune.smry |> 
#   filter(median.cvm == min(median.cvm)) |> 
#   as.data.frame()
# min.median.cvm
# 
# ggplot(glmnet.tune.smry, aes(alpha, median.cvm)) +
#   geom_line() +
#   geom_vline(xintercept = min.median.cvm$alpha) +
#   geom_point(aes(fill = median.lambda.1se), shape = 21) +
#   xlim(c(0.001, 0.04)) +
#   scale_fill_distiller(palette = 'Spectral') +
#   labs(x = 'alpha', y = 'median.cvm')

# number and % of sites retained for each value of alpha

# glmnet.site.by.alpha <- lapply(glmnet.tune, function(tune.rep) {
#   tune.rep |> 
#     lapply(function(alpha.rep) {
#       data.frame(alpha = alpha.rep$smry$alpha, loc.site = alpha.rep$loc.site)
#     }) |> 
#     bind_rows() 
# }) |> 
#   bind_rows() |> 
#   group_by(alpha, loc.site) |> 
#   summarize(prop = 100 * (n() / nrep), .groups = 'drop')
# 
# gridExtra::grid.arrange(
#   glmnet.site.by.alpha |> 
#     filter(prop > 0.95) |> 
#     group_by(alpha) |> 
#     summarize(num.sites = n_distinct(loc.site), .groups = 'drop') |> 
#     ggplot() +
#     geom_line(aes(alpha, num.sites)) + 
#     geom_vline(xintercept = min.median.cvm$alpha),
#   
#   glmnet.site.by.alpha |> 
#     left_join(locus.map, by = 'loc.site') |> 
#     ggplot() +
#     geom_tile(aes(loc.site, alpha, fill = prop)) +
#     geom_hline(yintercept = min.median.cvm$alpha) + 
#     scale_fill_distiller(palette = 'YlOrRd', direction = 1) +
#     facet_wrap(~ locus, scales = 'free') +
#     theme_minimal() +
#     theme(
#       axis.text.x = element_text(angle = 90, hjust = 0),
#       legend.position = 'top'
#     ),
#   heights = c(1, 4)
# )

# glmnet.site.by.alpha |> 
#   filter(alpha == min.median.cvm$alpha) |> 
#   pull('loc.site') |> 
#   saveRDS('glmnet sites.rds')





# use caret::train

# caret.tune <- caret::train(
#   age.best ~ .,
#   data = model.ages[, c('age.best', sites.to.keep)],
#   method = 'glmnet',
#   weights = model.ages$confidence.wt,
#   trControl = trainControl(method = 'cv'),
#   tuneGrid = expand.grid(alpha = alpha.vec, lambda = seq(30, 120, by = 0.5)),
#   tuneLength = 100
# )
# 
# caret.tuned.best <- cv.glmnet(
#   x = as.matrix(model.ages[, sites.to.keep]),
#   y = model.ages$age.best,
#   weights = model.ages$confidence.wt,
#   alpha = caret.tune$results |> 
#     rownames_to_column('rep') |> 
#     filter(rep == rownames(caret.tune$bestTune)) |> 
#     pull('alpha')
# ) 
# print(caret.tuned.best)
# plot(caret.tuned.best)