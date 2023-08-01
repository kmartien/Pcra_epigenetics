rm(list = ls())
library(tidyverse)
library(abind)

load("Pcra.epi.data.for.Eric.Rdata")

# combineList <- function(x, name) {
#   x <- x[[name]]
#   if(is.null(x)) return(NULL)
#   abind(
#     coverage = x$coverage %>% 
#       column_to_rownames("id") %>% 
#       select(-avg) %>% 
#       as.matrix(),
#     freq.meth = x$freq.meth %>% 
#       column_to_rownames("id") %>% 
#       as.matrix(),
#     errors = x$errors %>% 
#       column_to_rownames("id") %>% 
#       as.matrix(),
#     along = 3
#   )
# }
# cpg <- lapply(meth.dat, combineList, name = "CpG.sum")
# non.cpg <- lapply(meth.dat, combineList, name = "non.CpG.sum")

epi.df <- do.call(
  rbind,
  lapply(names(meth.dat), function(locus) {
    do.call(
      rbind,
      lapply(names(meth.dat[[locus]]), function(type) {
        x <- meth.dat[[locus]][[type]]
        if(is.null(x)) return(NULL)
        left_join(
          x$coverage %>% 
            select(-avg) %>% 
            pivot_longer(-id, names_to = "site", values_to = "coverage"),
          pivot_longer(x$freq.meth, -id, names_to = "site", values_to = "freq.meth"),
          by = c("id", "site")
        ) %>% 
          left_join(
            pivot_longer(x$errors, -id, names_to = "site", values_to = "errors"),
            by = c("id", "site")
          ) %>% 
          mutate(
            locus = locus,
            type = type
          )
      })
    )
  })
) %>% 
  arrange(type, id, locus, site) %>% 
  left_join(
    enframe(
      setNames(pct.conversion, meth.dat[[1]][[1]]$coverage$id),
      name = "id",
      value = "pct.conversion"
    ),
    by = "id"
  ) %>% 
  mutate(site = as.numeric(site))

saveRDS(epi.df, file = "methylation df.rdata")