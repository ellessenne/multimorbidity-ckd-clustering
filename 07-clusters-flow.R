### Script to subset subjects that are in common within cluster
library(tidyverse)
library(ggalluvial)
library(matrixStats)
library(ggrepel)
library(readxl)

# Custom, colourblind-friendly palette
# This is tuned for 6 clusters...
.palette <- c(
  vermillion = "#D55E00",
  sky_blue = "#56B4E9",
  bluish_green = "#009E73",
  reddish_purple = "#CC79A7",
  yellow = "#F0E442",
  orange = "#E69F00",
  blue = "#0072B2",
  black = "#000000"
)
names(.palette) <- NULL
# Import multimorbidity IDs
mm <- readRDS(file = "data/03-mm.RDS")
subdf <- map(.x = seq_along(mm), .f = function(i) {
  mm[[i]] %>%
    select(lopnr)
})

# Reduce
subdf <- reduce(.x = subdf, .f = inner_join, by = "lopnr")

# Recover multimorbidity data
mm <- readRDS(file = "data/03-mm.RDS")
mm_full <- map_dfr(.x = seq_along(mm), .f = function(i) {
  out <- mm[[i]] %>%
    mutate(eGFR = names(mm)[i])
  out$mm <- out %>%
    select(-lopnr, -eGFR, -ckd) %>%
    as.matrix() %>%
    rowSums2()
  out %>%
    mutate(mmc = case_when(
      mm >= 8 ~ "8+",
      TRUE ~ as.character(mm)
    )) %>%
    mutate(mmc = factor(mmc, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8+")))
}) %>%
  mutate(eGFR = factor(eGFR, levels = c("90", "75", "60", "45", "30", "15")))
mm_common <- map_dfr(.x = seq_along(mm), .f = function(i) {
  out <- mm[[i]] %>%
    filter(lopnr %in% subdf$lopnr) %>%
    mutate(eGFR = names(mm)[i])
  out$mm <- out %>%
    select(-lopnr, -eGFR, -ckd) %>%
    as.matrix() %>%
    rowSums2()
  out %>%
    mutate(mmc = case_when(
      mm >= 8 ~ "8+",
      TRUE ~ as.character(mm)
    )) %>%
    mutate(mmc = factor(mmc, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8+")))
}) %>%
  mutate(eGFR = factor(eGFR, levels = c("90", "75", "60", "45", "30", "15")))

# Make alluvial plot (?)
pa_full <- mm_full %>%
  ggplot(aes(x = eGFR, alluvium = lopnr, stratum = mmc, fill = mmc, color = mmc)) +
  geom_flow() +
  geom_stratum(alpha = 0.9) +
  scale_fill_viridis_d(drop = FALSE) +
  scale_color_viridis_d(drop = FALSE) +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  labs(y = "Number of Subjects", fill = "", color = "", caption = "Number of concurrent morbidities across categories of eGFR")
ggsave(plot = pa_full, filename = "figures/07-alluvial-full.pdf", height = 5, width = 5 * sqrt(2), device = cairo_pdf)
ggsave(plot = pa_full, filename = "figures/07-alluvial-full.png", height = 5, width = 5 * sqrt(2), device = ragg::agg_png, res = 600, units = "in")

pa_common <- mm_common %>%
  ggplot(aes(x = eGFR, alluvium = lopnr, stratum = mmc, fill = mmc, color = mmc, label = mmc)) +
  geom_flow() +
  geom_stratum(alpha = 0.9) +
  geom_text(aes(family = "Arial Narrow"), stat = "stratum", color = "white") +
  scale_fill_viridis_d(drop = FALSE) +
  scale_color_viridis_d(drop = FALSE) +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(legend.position = "none") +
  labs(y = "Number of Subjects", fill = "", caption = "Number of concurrent morbidities across categories of eGFR")
ggsave(plot = pa_common, filename = "figures/07-alluvial-common.pdf", height = 5, width = 5 * sqrt(2), device = cairo_pdf)
ggsave(plot = pa_common, filename = "figures/07-alluvial-common.png", height = 5, width = 5 * sqrt(2), device = ragg::agg_png, res = 600, units = "in")

# Here we make alluvial plots for actual cluster descriptions
fc <- readRDS("data/04-fitted-clusters.RDS")
cdesc <- read_excel("data/05-SCREAM-cluster-descriptions-v5.xlsx") %>%
  mutate(
    description = tools::toTitleCase(description),
    system = tools::toTitleCase(system),
    eGFR = as.character(eGFR)
  )
n_to_pick <- c("7", "9", "7", "6", "6", "5")

# All subjects
fc_full <- map_dfr(.x = seq_along(fc), .f = function(i) {
  out <- mm[[i]]
  out$cluster <- fc[[i]][[n_to_pick[i]]]$cluster
  out <- out %>%
    mutate(eGFR = names(mm)[i]) %>%
    left_join(cdesc, by = c("eGFR", "cluster")) %>%
    select(lopnr, eGFR, cluster, description, system)
}) %>%
  mutate(
    eGFR = factor(eGFR, levels = c("90", "75", "60", "45", "30", "15")),
    system = factor(system, levels = c("Cardiovascular", "Endocrine", "Cancer", "Mental Health, Pain", "Respiratory", "Non-Specific"))
  )
fc_full <- fc_full %>%
  ggplot(aes(x = eGFR, alluvium = lopnr, stratum = description, fill = system, color = system, label = description)) +
  geom_flow() +
  geom_stratum(alpha = 2 / 3) +
  geom_text(aes(family = "Arial Narrow"), stat = "stratum", color = "black", size = 3) +
  scale_fill_manual(values = .palette) +
  scale_color_manual(values = .palette) +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(legend.position = "none") +
  labs(y = "Number of Subjects", fill = "")
ggsave(plot = fc_full, filename = "figures/07-alluvial-fc-full.pdf", height = 6, width = 6 * sqrt(2), device = cairo_pdf)
ggsave(plot = fc_full, filename = "figures/07-alluvial-fc-full.png", height = 6, width = 6 * sqrt(2), device = ragg::agg_png, res = 600, units = "in")
# Then, common subjects
fc_common <- map_dfr(.x = seq_along(fc), .f = function(i) {
  out <- mm[[i]]
  out$cluster <- fc[[i]][[n_to_pick[i]]]$cluster
  out <- out %>%
    mutate(eGFR = names(mm)[i]) %>%
    left_join(cdesc, by = c("eGFR", "cluster")) %>%
    filter(lopnr %in% subdf$lopnr) %>%
    select(lopnr, eGFR, cluster, description, system)
}) %>%
  mutate(
    eGFR = factor(eGFR, levels = c("90", "75", "60", "45", "30", "15")),
    system = factor(system, levels = c("Cardiovascular", "Endocrine", "Cancer", "Mental Health, Pain", "Respiratory", "Non-Specific"))
  )
fc_common <- fc_common %>%
  ggplot(aes(x = eGFR, alluvium = lopnr, stratum = description, fill = system, color = system, label = description)) +
  geom_flow() +
  geom_stratum(alpha = 2 / 3) +
  geom_text(aes(family = "Arial Narrow"), stat = "stratum", color = "black", size = 3) +
  scale_fill_manual(values = .palette) +
  scale_color_manual(values = .palette) +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(legend.position = "none") +
  labs(y = "Number of Subjects", fill = "")
ggsave(plot = fc_common, filename = "figures/07-alluvial-fc-common.pdf", height = 6, width = 6 * sqrt(2), device = cairo_pdf)
ggsave(plot = fc_common, filename = "figures/07-alluvial-fc-common.png", height = 6, width = 6 * sqrt(2), device = ragg::agg_png, res = 600, units = "in")

# Some other options for plotting...
p_prop <- map_dfr(.x = seq_along(fc), .f = function(i) {
  out <- mm[[i]]
  out$cluster <- fc[[i]][[n_to_pick[i]]]$cluster
  out <- out %>%
    mutate(eGFR = names(mm)[i]) %>%
    left_join(cdesc, by = c("eGFR", "cluster")) %>%
    select(lopnr, eGFR, cluster, description, system)
}) %>%
  mutate(
    eGFR = factor(eGFR, levels = c("90", "75", "60", "45", "30", "15")),
    system = factor(system, levels = c("Cardiovascular", "Endocrine", "Cancer", "Mental Health, Pain", "Respiratory", "Non-Specific"))
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
  ungroup() %>%
  ggplot(aes(x = eGFR, y = p, fill = system, color = system, group = ord)) +
  geom_col(alpha = 1 / 2, width = 9 / 10) +
  geom_label_repel(aes(label = str_wrap(description, width = 20)), color = "black", position = position_stack(vjust = 0.5), size = 2, alpha = 4 / 5, show.legend = FALSE, box.padding = 0.1, label.padding = 0.1, force_pull = 1, min.segment.length = 0) +
  scale_fill_manual(values = .palette) +
  scale_color_manual(values = .palette) +
  scale_y_continuous(labels = function(x) formattable::percent(x, 0)) +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(legend.position = "bottom") +
  labs(x = "eGFR", y = "Cluster-Wise Proportion", fill = "", color = "")
ggsave(plot = p_prop, filename = "figures/07-cluster-wise-proportion.pdf", height = 5, width = 5 * sqrt(2), device = cairo_pdf)
ggsave(plot = p_prop, filename = "figures/07-cluster-wise-proportion.png", height = 5, width = 5 * sqrt(2), device = ragg::agg_png, res = 600, units = "in")
