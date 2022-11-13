# Age-specific analysis for age >= 65
# At each interpolated date, we calculate age and select only subjects with age >= 65
# Then, we run the clustering algorithm and create some output
library(klaR)
library(glue)
library(formattable)
library(knitr)
library(patchwork)
library(tidyverse)
library(ggrepel)

# Import previous data
data <- readRDS(file = "path/to/file")
data <- distinct(data, lopnr, dob)
mm <- readRDS(file = "data/03-mm.RDS")
ip <- readRDS(file = "data/02-interpolated-index-dates.RDS")

# Age at index date
ip <- map(.x = ip, .f = function(x) {
  left_join(x, data, by = "lopnr") %>%
    mutate(age = as.numeric(index_date - dob) / 365.242) %>%
    filter(age >= 65) %>%
    pull(lopnr)
})

# Reduce size of mm
mm <- map(.x = seq_along(mm), .f = function(i) filter(mm[[i]], lopnr %in% ip[[i]]))
names(mm) <- names(ip)

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
  labs(x = "Number of Clusters", y = "Average Within Cluster Distance", caption = "Analysis for age >= 65")
ggsave(plot = ep, filename = "figures/41-age-specific-above-elbow-plots.pdf", height = 5, width = 5 * 29.7 / 21, device = cairo_pdf)
ggsave(plot = ep, filename = "figures/41-age-specific-above-elbow-plots.png", height = 5, width = 5 * 29.7 / 21, device = ragg::agg_png, res = 600, units = "in")
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
  labs(x = "Number of Clusters", y = "Gradient of Elbow Plot", caption = "Analysis for age >= 65")
ggsave(plot = dep, filename = "figures/41-age-specific-above-elbow-plots-grad.pdf", height = 5, width = 5 * 29.7 / 21, device = cairo_pdf)
ggsave(plot = dep, filename = "figures/41-age-specific-above-elbow-plots-grad.png", height = 5, width = 5 * 29.7 / 21, device = ragg::agg_png, res = 600, units = "in")

# Export all clustering results too
saveRDS(object = full_res, file = "data/41-age-specific-above-fitted-clusters.RDS")

# Summarise the optimal clustering
# We pick 'n_to_pick' clusters for each eGFR value
n_to_pick <- c("4", "4", "4", "5", "6", "7")
full_res <- map(.x = seq_along(full_res), .f = function(i) full_res[[i]][[n_to_pick[i]]])
names(full_res) <- names(mm)

# Now, we add clustering information to each multimorbidity dataset
summ_mm <- map_dfr(.x = seq_along(mm), .f = function(x) {
  df <- mm[[x]]
  df[["lopnr"]] <- NULL
  df[["cluster"]] <- full_res[[x]]$cluster
  df[["eGFR"]] <- names(full_res)[x]
  df
}) %>%
  group_by(eGFR, cluster) %>%
  summarise(across(.fns = mean)) %>%
  ungroup() %>%
  select(-ckd) %>%
  pivot_longer(cols = 3:29) %>%
  rename(Condition = name) %>%
  mutate(Condition = factor(Condition, levels = c("alcohol_misuse", "asthma", "afib", "cancer", "ckd", "cpain", "cvhepatitis", "cirrhosis", "dementia", "depression", "diabetes", "epilepsy", "chf", "hypertension", "hypothyroidism", "ibd", "ibs", "multiple_sclerosis", "mi", "parkinson", "cpd", "pud", "pvd", "psoriasis", "rheum_arthritis", "schizofrenia", "severe_constipation", "stroke"), labels = c("Alcohol Misuse", "Asthma", "Atrial Fibrillation", "Cancer", "CKD", "Chronic Pain", "Chronic Viral Hepatitis B", "Cirrhosis", "Dementia", "Depression", "Diabetes", "Epilepsy", "Heart Failure", "Hypertension", "Hypothyroidism", "IBD", "IBS", "Multiple Sclerosis", "MI", "Parkinson's Disease", "Pulmonary Disease", "PUD", "PVD", "Psoriasis", "Rheumatoid Arthritis", "Schizophrenia", "Severe Constipation", "Stroke"))) %>%
  arrange(eGFR, Condition)
# Export this data
summ_mm %>%
  pivot_wider(names_from = cluster, values_from = value) %>%
  arrange(eGFR, Condition) %>%
  write_csv(file = "data/41-age-specific-above-proportion-conditions-by-clusters.csv")
# Heatmap
p_heat <- summ_mm %>%
  mutate(Condition = forcats::fct_rev(Condition)) %>%
  ggplot(aes(x = factor(cluster), y = Condition, fill = value)) +
  geom_tile() +
  facet_wrap(. ~ eGFR, labeller = label_both, scales = "free_x") +
  scale_fill_viridis_c() +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(legend.position = "top") +
  labs(x = "Cluster", fill = "Proportion with Condition", caption = "Analysis for age >= 65")
ggsave(plot = p_heat, filename = "figures/41-age-specific-above-heatmap-by-clusters.pdf", width = 6, height = 6 * sqrt(2), device = cairo_pdf)
ggsave(plot = p_heat, filename = "figures/41-age-specific-above-heatmap-by-clusters.png", width = 6, height = 6 * sqrt(2), device = ragg::agg_png, res = 600, units = "in")

# Calculate clusters that are "over-expressed"
overall_prevalences <- map_dfr(.x = seq_along(mm), .f = function(i) {
  mm[[i]] %>%
    select(-lopnr, -ckd) %>%
    summarise(across(.fns = mean)) %>%
    mutate(eGFR = names(mm)[i]) %>%
    pivot_longer(cols = -eGFR)
}) %>%
  mutate(name = factor(name, levels = c("alcohol_misuse", "asthma", "afib", "cancer", "cpain", "cvhepatitis", "cirrhosis", "dementia", "depression", "diabetes", "epilepsy", "chf", "hypertension", "hypothyroidism", "ibd", "ibs", "multiple_sclerosis", "mi", "parkinson", "cpd", "pud", "pvd", "psoriasis", "rheum_arthritis", "schizofrenia", "severe_constipation", "stroke"), labels = c("Alcohol Misuse", "Asthma", "Atrial Fibrillation", "Cancer", "Chronic Pain", "Chronic Viral Hepatitis B", "Cirrhosis", "Dementia", "Depression", "Diabetes", "Epilepsy", "Heart Failure", "Hypertension", "Hypothyroidism", "IBD", "IBS", "Multiple Sclerosis", "MI", "Parkinson's Disease", "Pulmonary Disease", "PUD", "PVD", "Psoriasis", "Rheumatoid Arthritis", "Schizophrenia", "Severe Constipation", "Stroke"))) %>%
  arrange(eGFR, name)
oedf <- pivot_wider(data = summ_mm, names_from = "cluster", values_from = "value")
oedf <- split(x = oedf, f = oedf$eGFR)
oedf <- map(.x = seq_along(oedf), .f = function(i) {
  n <- oedf[[i]]$Condition
  x <- select(oedf[[i]], -eGFR, -Condition) %>%
    as.matrix()
  rownames(x) <- n
  x
})
oedf <- map(.x = seq_along(oedf), .f = function(i) {
  x <- oedf[[i]]
  prev <- (x >= 0.2)
  avgs <- filter(overall_prevalences, eGFR == names(mm)[i]) %>%
    select(-eGFR, -name) %>%
    pull()
  avgs <- matrix(data = avgs, nrow = nrow(x), ncol = ncol(x), byrow = FALSE)
  over <- x / avgs
  over[is.nan(over)] <- 0
  over <- (over >= 2)
  out <- as.data.frame(over * prev)
  out$Condition <- rownames(x)
  rownames(out) <- NULL
  out <- out %>%
    mutate(Condition = factor(Condition)) %>%
    mutate(Condition = forcats::fct_rev(Condition))
  out
})
names(oedf) <- names(full_res)
oedf <- map_dfr(.x = seq_along(oedf), .f = function(i) {
  oedf[[i]] %>%
    mutate(eGFR = names(oedf)[i])
}) %>%
  pivot_longer(cols = 1:max(as.numeric(n_to_pick))) %>%
  filter(!is.na(value))
p_def <- ggplot(oedf, aes(x = name, y = Condition, fill = factor(value))) +
  geom_tile(show.legend = FALSE, alpha = 2 / 3) +
  facet_wrap(~eGFR, as.table = TRUE, labeller = label_both, scales = "free_x") +
  scale_discrete_manual(aesthetics = "fill", values = c("grey90", "red")) +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(legend.position = "top") +
  labs(x = "Cluster", fill = "Proportion with Condition", caption = "Analysis for age >= 65")
ggsave(plot = p_def, filename = "figures/41-age-specific-above-defining-conditions-by-clusters.pdf", width = 6, height = 6 * sqrt(2), device = cairo_pdf)
ggsave(plot = p_def, filename = "figures/41-age-specific-above-defining-conditions-by-clusters.png", width = 6, height = 6 * sqrt(2), device = ragg::agg_png, res = 600, units = "in")

# Cluster-wise proportions
.palette <- c(
  vermillion = "#D55E00",
  sky_blue = "#56B4E9",
  bluish_green = "#009E73",
  reddish_purple = "#CC79A7",
  orange = "#E69F00",
  blue = "#0072B2",
  yellow = "#F0E442",
  black = "#000000"
)
names(.palette) <- NULL
cdesc <- readxl::read_excel("data/41-age-specific-above-SCREAM-cluster-descriptions-v2.xlsx")
cdesc$eGFR <- as.character(cdesc$eGFR)
cwp <- map_dfr(.x = seq_along(mm), .f = function(i) {
  out <- mm[[i]]
  out$cluster <- full_res[[i]]$cluster
  out <- out %>%
    mutate(eGFR = names(mm)[i]) %>%
    left_join(cdesc, by = c("eGFR", "cluster")) %>%
    select(lopnr, eGFR, cluster, description, system)
}) %>%
  mutate(
    eGFR = factor(eGFR, levels = c("90", "75", "60", "45", "30", "15")),
    system = factor(system, levels = c("Cardiovascular", "Endocrine", "Cancer", "Mental Health, Pain", "Respiratory", "Dermatological", "Non-Specific"))
  ) %>%
  group_by(eGFR) %>%
  mutate(N = n()) %>%
  ungroup() %>%
  group_by(eGFR, system, description, N) %>%
  summarise(n = n()) %>%
  mutate(
    p = n / N,
    label = formattable::percent(p)
  ) %>%
  arrange(eGFR, system, description) %>%
  group_by(eGFR) %>%
  mutate(ord = row_number()) %>%
  ungroup()
p_prop <- ggplot(cwp, aes(x = eGFR, y = p, fill = system, color = system, group = ord)) +
  geom_col(alpha = 1 / 2, width = 9 / 10) +
  geom_label_repel(aes(label = str_wrap(description, width = 20)), color = "black", position = position_stack(vjust = 0.5), size = 2, alpha = 4 / 5, show.legend = FALSE, box.padding = 0.1, label.padding = 0.1, force_pull = 1, min.segment.length = 0) +
  scale_fill_manual(values = .palette) +
  scale_color_manual(values = .palette) +
  scale_y_continuous(labels = function(x) formattable::percent(x, 0)) +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(legend.position = "bottom") +
  labs(x = "eGFR", y = "Cluster-Wise Proportion", fill = "", color = "", caption = "Analysis for age >= 65")
ggsave(plot = p_prop, filename = "figures/41-age-specific-above-cluster-wise-proportion.pdf", height = 5, width = 5 * sqrt(2), device = cairo_pdf)
ggsave(plot = p_prop, filename = "figures/41-age-specific-above-cluster-wise-proportion.png", height = 5, width = 5 * sqrt(2), device = ragg::agg_png, res = 600, units = "in")
