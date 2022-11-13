### Here we run all the analysis for clustering
library(tidyverse)
library(klaR)
library(glue)

# Import multimorbidity data
mm <- readRDS(file = "data/03-mm.RDS")

# Remove IDs (and CKD) for clustering
cldata <- map(.x = mm, .f = function(x) dplyr::select(x, -lopnr, -ckd) %>% as.matrix())
names(cldata) <- names(mm)

# Define possible number of clusters
.n_clusters <- 2:10

# Over all eGFR values, fit all possible clusters combinations
fit_clusters <- function(data, possible_n) {
  res <- map(.x = possible_n, .f = function(n) {
    ccc <- NULL
    while (is.null(ccc)) {
      ccc <- tryCatch(
        {
          kmodes(data = data, modes = n, iter.max = 20)
        },
        error = function(e) NULL
      )
    }
    ccc
  })
  names(res) <- possible_n
  res
}
set.seed(398769832) # for reproducibility
pb <- txtProgressBar(min = 0, max = length(cldata), style = 3)
full_res <- map(.x = seq_along(cldata), .f = function(i) {
  y <- fit_clusters(data = cldata[[i]], possible_n = .n_clusters)
  setTxtProgressBar(pb, value = i)
  y
})
names(full_res) <- names(cldata)
close(pb)

# Make elbow plots
ep_res <- map_dfr(.x = seq_along(full_res), .f = function(i) {
  tibble(eGFR = names(full_res)[i], n = as.numeric(names(full_res[[i]])), mwd = map_dbl(.x = full_res[[i]], .f = function(x) mean(x$withindiff)))
})
ep <- ggplot(ep_res, aes(x = n, y = mwd)) +
  geom_line(linetype = "dotted") +
  geom_point() +
  facet_wrap(~eGFR, labeller = label_both, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(.n_clusters), ceiling(max(.n_clusters) / 2) * 2, by = 2)) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  labs(x = "Number of Clusters", y = "Average Within Cluster Distance")
ggsave(plot = ep, filename = "figures/04-elbow-plots.pdf", height = 5, width = 5 * 29.7 / 21, device = cairo_pdf)
ggsave(plot = ep, filename = "figures/04-elbow-plots.png", height = 5, width = 5 * 29.7 / 21, device = ragg::agg_png, res = 600, units = "in")
# Here I derive the gradients at each number to identify when the elbow plot flattens
grep_res <- map_dfr(.x = names(full_res), .f = function(e) {
  spfun <- splinefun(x = ep_res$n[ep_res$eGFR == e], y = ep_res$mwd[ep_res$eGFR == e])
  tibble(eGFR = e, n = .n_clusters, grad = spfun(.n_clusters, deriv = 1))
})
dep <- ggplot(grep_res, aes(x = n, y = grad)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_line(linetype = "dotted") +
  geom_point() +
  facet_wrap(~eGFR, labeller = label_both, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(.n_clusters), ceiling(max(.n_clusters) / 2) * 2, by = 2)) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  labs(x = "Number of Clusters", y = "Gradient of Elbow Plot")
ggsave(plot = dep, filename = "figures/04-elbow-plots-grad.pdf", height = 5, width = 5 * 29.7 / 21, device = cairo_pdf)
ggsave(plot = dep, filename = "figures/04-elbow-plots-grad.png", height = 5, width = 5 * 29.7 / 21, device = ragg::agg_png, res = 600, units = "in")

# Export all clustering results too
saveRDS(object = full_res, file = "data/04-fitted-clusters.RDS")

# As a sensitivity analysis, re-running the algorithm a bunch of times
# to test the choice of clusters. To do this, we run each clustering
# algorithm 'times' more times and create new elbow plots
epfun <- function(data, possible_n, times, eGFR, .pb) {
  out <- map_dfr(.x = seq(times), .f = function(i) {
    y <- fit_clusters(data, possible_n)
    y <- map_dbl(.x = y, .f = function(x) {
      setTxtProgressBar(pb = .pb, value = .pb$getVal() + 1)
      mean(x$withindiff)
    })
    tibble::tibble(i = i, wd = y, k = possible_n)
  })
  out[["eGFR"]] <- eGFR
  return(out)
}
.times <- 20
pb <- txtProgressBar(max = length(.n_clusters) * length(cldata) * .times, style = 3)
test <- map(.x = seq_along(cldata), .f = function(j) epfun(cldata[[j]], possible_n = .n_clusters, times = .times, eGFR = names(cldata)[j], .pb = pb))
close(pb)
testb <- map_dfr(.x = seq_along(test), .f = function(i) {
  ff <- splinefun(x = test[[i]]$k, y = test[[i]]$wd, method = "natural")
  tibble(
    k = sort(unique(test[[i]]$k)),
    wd_hat = ff(x = k),
    wd_grad = ff(x = k, deriv = 1L),
    eGFR = names(cldata)[i]
  )
})
tep <- ggplot(testb, aes(x = k, y = wd_hat)) +
  geom_line(linetype = "dotted") +
  geom_point() +
  facet_wrap(~eGFR, labeller = label_both, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(.n_clusters), ceiling(max(.n_clusters) / 2) * 2, by = 2)) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  labs(x = "Number of Clusters", y = "Average Within Cluster Distance", caption = glue("Plots smoothed over {.times} repetitions of the clustering algorithm."))
ggsave(plot = tep, filename = "figures/04-elbow-plots-test.pdf", height = 5, width = 5 * 29.7 / 21, device = cairo_pdf)
ggsave(plot = tep, filename = "figures/04-elbow-plots-test.png", height = 5, width = 5 * 29.7 / 21, device = ragg::agg_png, res = 600, units = "in")
tdep <- ggplot(testb, aes(x = k, y = wd_grad)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
  geom_line(linetype = "dotted") +
  geom_point() +
  facet_wrap(~eGFR, labeller = label_both, scales = "free_y") +
  scale_x_continuous(breaks = seq(min(.n_clusters), ceiling(max(.n_clusters) / 2) * 2, by = 2)) +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  labs(x = "Number of Clusters", y = "Gradient of Elbow Plot", caption = glue("Plots smoothed over {.times} repetitions of the clustering algorithm."))
ggsave(plot = tdep, filename = "figures/04-elbow-plots-grad-test.pdf", height = 5, width = 5 * 29.7 / 21, device = cairo_pdf)
ggsave(plot = tdep, filename = "figures/04-elbow-plots-grad-test.png", height = 5, width = 5 * 29.7 / 21, device = ragg::agg_png, res = 600, units = "in")
