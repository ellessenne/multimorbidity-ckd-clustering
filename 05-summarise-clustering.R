### Here we summarise the optimal clustering
library(tidyverse)
library(glue)
library(formattable)
library(knitr)
library(patchwork)

# Import data
mm <- readRDS(file = "data/03-mm.RDS")
.eGFR <- names(mm)
ipdata <- readRDS(file = "data/02-interpolated-index-dates.RDS")
demo <- readRDS(file = "path/to/file")
demo <- select(demo, lopnr, female, dob) %>%
  distinct()
# New columns:
# - 'female' = 1 if females, 0 otherwise
# - 'dob' is date of birth

# Import clustering data
fitted <- readRDS(file = "data/04-fitted-clusters.RDS")

# We pick 'n_to_pick' clusters for each eGFR value
# This was chosen subjectively by the analysts
# Make sure you choose this properly!
n_to_pick <- c("7", "9", "7", "6", "6", "5")
fitted <- map(.x = seq_along(fitted), .f = function(i) fitted[[i]][[n_to_pick[i]]])
names(fitted) <- names(mm)

# Now, we add clustering information to each multimorbidity dataset
mm <- map(.x = seq_along(mm), .f = function(x) {
  df <- mm[[x]]
  df[["cluster"]] <- fitted[[x]]$cluster
  df <- left_join(df, ipdata[[x]], by = "lopnr") %>%
    left_join(demo, by = "lopnr") %>%
    mutate(age = as.numeric(index_date - dob) / 365.242) %>%
    select(-index_date, -dob)
})
names(mm) <- names(fitted)

# Summarise by cluster
# Proportions, this is more relevant
summ_mm <- map_dfr(.x = seq_along(mm), .f = function(i) {
  mm[[i]] %>%
    select(-lopnr) %>%
    group_by(cluster) %>%
    summarise(across(.fns = mean), n = n()) %>%
    ungroup() %>%
    pivot_longer(cols = -cluster) %>%
    mutate(value = case_when(
      name == "age" ~ paste0(formattable::comma(value)),
      name == "n" ~ paste0(formattable::comma(value, 0)),
      TRUE ~ paste0(formattable::percent(value))
    )) %>%
    pivot_wider(names_from = "cluster", values_from = "value", names_prefix = "Cluster_") %>%
    mutate(eGFR = .eGFR[i]) %>%
    rename(Condition = name)
}) %>%
  mutate(Condition = factor(Condition, levels = c("n", "age", "female", "alcohol_misuse", "asthma", "afib", "cancer", "ckd", "cpain", "cvhepatitis", "cirrhosis", "dementia", "depression", "diabetes", "epilepsy", "chf", "hypertension", "hypothyroidism", "ibd", "ibs", "multiple_sclerosis", "mi", "parkinson", "cpd", "pud", "pvd", "psoriasis", "rheum_arthritis", "schizofrenia", "severe_constipation", "stroke"), labels = c("n", "Age", "Gender = Female", "Alcohol Misuse", "Asthma", "Atrial Fibrillation", "Cancer", "CKD", "Chronic Pain", "Chronic Viral Hepatitis B", "Cirrhosis", "Dementia", "Depression", "Diabetes", "Epilepsy", "Heart Failure", "Hypertension", "Hypothyroidism", "IBD", "IBS", "Multiple Sclerosis", "MI", "Parkinson's Disease", "Pulmonary Disease", "PUD", "PVD", "Psoriasis", "Rheumatoid Arthritis", "Schizophrenia", "Severe Constipation", "Stroke"))) %>%
  arrange(eGFR, Condition) %>%
  select(eGFR, Condition, starts_with("Cluster"))
write_excel_csv(x = summ_mm, file = "data/05-proportion-conditions-by-clusters.csv")
# Numbers, this is useful later
summ_mm_n <- map_dfr(.x = seq_along(mm), .f = function(i) {
  mm[[i]] %>%
    select(-lopnr, -age, -female) %>%
    group_by(cluster) %>%
    summarise(across(.fns = sum), n = n()) %>%
    ungroup() %>%
    pivot_longer(cols = -cluster) %>%
    mutate(value = paste0(formattable::comma(value, 0))) %>%
    pivot_wider(names_from = "cluster", values_from = "value", names_prefix = "Cluster_") %>%
    mutate(eGFR = .eGFR[i]) %>%
    rename(Condition = name)
}) %>%
  mutate(Condition = factor(Condition, levels = c("n", "age", "female", "alcohol_misuse", "asthma", "afib", "cancer", "ckd", "cpain", "cvhepatitis", "cirrhosis", "dementia", "depression", "diabetes", "epilepsy", "chf", "hypertension", "hypothyroidism", "ibd", "ibs", "multiple_sclerosis", "mi", "parkinson", "cpd", "pud", "pvd", "psoriasis", "rheum_arthritis", "schizofrenia", "severe_constipation", "stroke"), labels = c("n", "Age", "Gender = Female", "Alcohol Misuse", "Asthma", "Atrial Fibrillation", "Cancer", "CKD", "Chronic Pain", "Chronic Viral Hepatitis B", "Cirrhosis", "Dementia", "Depression", "Diabetes", "Epilepsy", "Heart Failure", "Hypertension", "Hypothyroidism", "IBD", "IBS", "Multiple Sclerosis", "MI", "Parkinson's Disease", "Pulmonary Disease", "PUD", "PVD", "Psoriasis", "Rheumatoid Arthritis", "Schizophrenia", "Severe Constipation", "Stroke"))) %>%
  arrange(eGFR, Condition) %>%
  select(eGFR, Condition, starts_with("Cluster"))
write_excel_csv(x = summ_mm_n, file = "data/05-number-conditions-by-clusters.csv")

sink(file = "data/05-tables-by-egfr.md")
for (e in .eGFR) {
  cat("eGFR:", e)
  kable(
    x = filter(summ_mm, eGFR == e) %>% select(-eGFR),
    format = "markdown"
  ) %>%
    print()
  cat("\n\n")
}
sink()

# Make a plot with the information above
plot_mm <- map_dfr(.x = seq_along(mm), .f = function(i) {
  mm[[i]] %>%
    select(-lopnr, -age, -female) %>%
    group_by(cluster) %>%
    summarise(across(.fns = mean)) %>%
    ungroup() %>%
    pivot_longer(cols = -cluster) %>%
    mutate(eGFR = .eGFR[i]) %>%
    # left_join(cdesc[[i]], by = "cluster") %>%
    rename(Cluster = cluster)
}) %>%
  mutate(name = factor(name, levels = c("n", "age", "female", "alcohol_misuse", "asthma", "afib", "cancer", "ckd", "cpain", "cvhepatitis", "cirrhosis", "dementia", "depression", "diabetes", "epilepsy", "chf", "hypertension", "hypothyroidism", "ibd", "ibs", "multiple_sclerosis", "mi", "parkinson", "cpd", "pud", "pvd", "psoriasis", "rheum_arthritis", "schizofrenia", "severe_constipation", "stroke"), labels = c("n", "Age", "Gender = Female", "Alcohol Misuse", "Asthma", "Atrial Fibrillation", "Cancer", "CKD", "Chronic Pain", "Chronic Viral Hepatitis B", "Cirrhosis", "Dementia", "Depression", "Diabetes", "Epilepsy", "Heart Failure", "Hypertension", "Hypothyroidism", "IBD", "IBS", "Multiple Sclerosis", "MI", "Parkinson's Disease", "Pulmonary Disease", "PUD", "PVD", "Psoriasis", "Rheumatoid Arthritis", "Schizophrenia", "Severe Constipation", "Stroke")))
p <- ggplot(plot_mm, aes(x = name, y = value)) +
  geom_col() +
  facet_grid(Cluster ~ eGFR, labeller = label_both) +
  scale_y_continuous(labels = formattable::percent) +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(x = "Condition", y = "Proportion")
ggsave(plot = p, filename = "figures/05-conditions-by-clusters.pdf", height = 13, width = 13 * sqrt(2), device = cairo_pdf)
ggsave(plot = p, filename = "figures/05-conditions-by-clusters.png", height = 13, width = 13 * sqrt(2), device = ragg::agg_png, res = 600, units = "in")

# Another plot with a 'heatmap':
hmdf <- readr::read_csv("data/05-proportion-conditions-by-clusters.csv") %>%
  filter(!(Condition %in% c("n", "Age", "Gender = Female", "CKD"))) %>%
  pivot_longer(cols = starts_with("Cluster_")) %>%
  mutate(value = stringr::str_remove_all(string = value, pattern = "%")) %>%
  mutate(name = stringr::str_remove_all(string = name, pattern = "Cluster_")) %>%
  mutate(value = readr::parse_double(value)) %>%
  mutate(Condition = factor(Condition)) %>%
  mutate(Condition = forcats::fct_rev(Condition)) %>%
  filter(!is.na(value))
p_heat <- ggplot(hmdf, aes(x = name, y = Condition, fill = value)) +
  geom_tile() +
  facet_wrap(. ~ eGFR, labeller = label_both, scales = "free_x") +
  scale_fill_viridis_c() +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(legend.position = "top") +
  labs(x = "Cluster", fill = "Proportion with Condition")
ggsave(plot = p_heat, filename = "figures/05-heatmap-by-clusters.pdf", width = 6, height = 6 * sqrt(2), device = cairo_pdf)
ggsave(plot = p_heat, filename = "figures/05-heatmap-by-clusters.png", width = 6, height = 6 * sqrt(2), device = ragg::agg_png, res = 600, units = "in")

# Calculate clusters that are "over-expressed"
overall_prevalences <- map_dfr(.x = seq_along(mm), .f = function(i) {
  mm[[i]] %>%
    select(-lopnr, -female, -age, -cluster, -ckd) %>%
    summarise(across(.fns = mean)) %>%
    mutate(eGFR = .eGFR[i]) %>%
    pivot_longer(cols = -eGFR)
}) %>%
  mutate(name = factor(name, levels = c("n", "age", "female", "alcohol_misuse", "asthma", "afib", "cancer", "ckd", "cpain", "cvhepatitis", "cirrhosis", "dementia", "depression", "diabetes", "epilepsy", "chf", "hypertension", "hypothyroidism", "ibd", "ibs", "multiple_sclerosis", "mi", "parkinson", "cpd", "pud", "pvd", "psoriasis", "rheum_arthritis", "schizofrenia", "severe_constipation", "stroke"), labels = c("n", "Age", "Gender = Female", "Alcohol Misuse", "Asthma", "Atrial Fibrillation", "Cancer", "CKD", "Chronic Pain", "Chronic Viral Hepatitis B", "Cirrhosis", "Dementia", "Depression", "Diabetes", "Epilepsy", "Heart Failure", "Hypertension", "Hypothyroidism", "IBD", "IBS", "Multiple Sclerosis", "MI", "Parkinson's Disease", "Pulmonary Disease", "PUD", "PVD", "Psoriasis", "Rheumatoid Arthritis", "Schizophrenia", "Severe Constipation", "Stroke"))) %>%
  arrange(eGFR, name)

oedf <- pivot_wider(data = hmdf, names_from = "name", values_from = "value")
oedf <- split(x = oedf, f = oedf$eGFR)
oedf <- map(.x = oedf, .f = function(x) {
  n <- x$Condition
  x <- select(x, -eGFR, -Condition) %>%
    as.matrix()
  rownames(x) <- n
  x
})
oedf <- map(.x = seq_along(oedf), .f = function(i) {
  x <- oedf[[i]]
  prev <- (x >= 20)
  avgs <- filter(overall_prevalences, eGFR == .eGFR[i]) %>%
    select(-eGFR, -name) %>%
    pull()
  avgs <- matrix(data = avgs * 100, nrow = nrow(x), ncol = ncol(x), byrow = FALSE)
  over <- x / avgs
  over <- (over >= 2)
  out <- as.data.frame(over * prev)
  out$Condition <- rownames(x)
  rownames(out) <- NULL
  out <- out %>%
    mutate(Condition = factor(Condition)) %>%
    mutate(Condition = forcats::fct_rev(Condition))
  out
})
names(oedf) <- names(fitted)
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
  labs(x = "Cluster", fill = "Proportion with Condition")
ggsave(plot = p_def, filename = "figures/05-defining-conditions-by-clusters.pdf", width = 6, height = 6 * sqrt(2), device = cairo_pdf)
ggsave(plot = p_def, filename = "figures/05-defining-conditions-by-clusters.png", width = 6, height = 6 * sqrt(2), device = ragg::agg_png, res = 600, units = "in")

# Merge p_heat and p_def
p_heat_def <- (p_heat + labs(y = "")) +
  (p_def + labs(y = "") + theme(axis.text.y = element_blank())) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", family = "Arial Narrow"))
ggsave(plot = p_heat_def, filename = "figures/05-heat-def.pdf", width = 9, height = 9, device = cairo_pdf)
ggsave(plot = p_heat_def, filename = "figures/05-heat-def.png", width = 9, height = 9, device = ragg::agg_png, res = 600, units = "in")
