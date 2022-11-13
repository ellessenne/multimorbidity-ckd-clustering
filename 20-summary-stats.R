# This script computes descriptive quantities and summary statistics for the ms
library(tidyverse)
library(formattable)
library(matrixStats)
library(survival)
library(lme4)
library(mvtnorm)

# Unique patients
base_data <- readRDS(file = "data/00-base-data.RDS")
base_data_d <- distinct(base_data, lopnr, female, dob)
sink(file = "data/20-overall.md")
cat("Unique patients from SCREAM:", paste0(comma(nrow(base_data), 0)), "\n\n")
sink()

# By eGFR
iid <- readRDS(file = "data/02-interpolated-index-dates.RDS")
iid <- map(.x = iid, .f = function(x) pull(x, lopnr))
iid <- c(iid$`15`, iid$`30`, iid$`45`, iid$`60`, iid$`75`, iid$`90`)
iid <- unique(iid)
iid <- left_join(data.frame(lopnr = iid), base_data_d, by = "lopnr")
sink(file = "data/20-overall.md", append = TRUE)
cat("Unique patients used with the clustering analysis:", paste0(comma(nrow(iid), 0)), "\n")
cat("Of which, females:", paste0(comma(sum(iid$female == 1), 0)), ",", paste0(percent(mean(iid$female == 1), 2)), "\n")
sink()

# Add cholesterol
chol <- readRDS(file = "data/00-chol.RDS") %>%
  filter(vtype != "IP") %>%
  select(-vtype) %>%
  group_by(lopnr, datum, test) %>%
  summarise(result = median(result)) %>%
  ungroup()
hdl <- chol %>%
  filter(test == "hdl") %>%
  select(-test) %>%
  rename(hdl = result)
ldl <- chol %>%
  filter(test == "ldl") %>%
  select(-test) %>%
  rename(ldl = result)
iid <- readRDS(file = "data/02-interpolated-index-dates.RDS")
hdl_ldl <- map(.x = iid, .f = function(x) {
  A <- hdl %>%
    left_join(x, by = "lopnr") %>%
    mutate(dd = abs(as.numeric(datum - index_date))) %>%
    filter(dd <= 365.242) %>%
    group_by(lopnr, index_date) %>%
    mutate(dd_min = min(dd)) %>%
    ungroup() %>%
    filter(dd == dd_min) %>%
    select(lopnr, hdl)
  B <- ldl %>%
    left_join(x, by = "lopnr") %>%
    mutate(dd = abs(as.numeric(datum - index_date))) %>%
    filter(dd <= 365.242) %>%
    group_by(lopnr, index_date) %>%
    mutate(dd_min = min(dd)) %>%
    ungroup() %>%
    filter(dd == dd_min) %>%
    select(lopnr, ldl)
  AB <- full_join(A, B, by = "lopnr")
  AB <- left_join(x, AB, by = "lopnr")
  return(AB)
})

# Table 1A
mm <- readRDS(file = "data/03-mm.RDS")
mmdf <- map(.x = mm, .f = function(x) {
  data.frame(
    lopnr = x$lopnr,
    mm = rowSums2(x = as.matrix(x[, -1]))
  )
})
table1a <- map_dfr(.x = seq_along(iid), .f = function(i) {
  indf <- left_join(iid[[i]], base_data_d, by = "lopnr") %>%
    left_join(mmdf[[i]], by = "lopnr") %>%
    mutate(age = as.numeric(index_date - dob) / 365.242) %>%
    left_join(hdl_ldl[[i]], by = c("lopnr", "index_date"))
  tibble(
    eGFR = names(iid)[i],
    `Number of patients` = paste0(comma(nrow(indf), 0)),
    `Age (years), Median (IQI)` = paste0(comma(median(indf$age), 2), " (", comma(fivenum(indf$age)[2], 2), " to ", comma(fivenum(indf$age)[4], 2), ")"),
    `Sex, Female (%)` = paste0(comma(sum(indf$female == 1), 0), " (", percent(mean(indf$female == 1), 2), ")"),
    `Condition  count (%), 0` = paste0(comma(sum(indf$mm == 0), 0), " (", percent(mean(indf$mm == 0), 2), ")"),
    `Condition  count (%), 1` = paste0(comma(sum(indf$mm == 1), 0), " (", percent(mean(indf$mm == 1), 2), ")"),
    `Condition  count (%), 2` = paste0(comma(sum(indf$mm == 2), 0), " (", percent(mean(indf$mm == 2), 2), ")"),
    `Condition  count (%), 3` = paste0(comma(sum(indf$mm == 3), 0), " (", percent(mean(indf$mm == 3), 2), ")"),
    `Condition  count (%), 4+` = paste0(comma(sum(indf$mm >= 4), 0), " (", percent(mean(indf$mm >= 4), 2), ")"),
    # `HDL, Median (IQI)` = paste0(comma(fivenum(indf$hdl)[3], 2), " (", comma(fivenum(indf$hdl)[2], 2), " to ", comma(fivenum(indf$hdl)[4], 2), ")"),
    # `HDL, Missing` = paste0(comma(sum(is.na(indf$hdl)), 0), " (", percent(mean(is.na(indf$hdl)), 2), ")"),
    # `LDL, Median (IQI)` = paste0(comma(fivenum(indf$ldl)[3], 2), " (", comma(fivenum(indf$ldl)[2], 2), " to ", comma(fivenum(indf$ldl)[4], 2), ")"),
    # `LDL, Missing` = paste0(comma(sum(is.na(indf$ldl)), 0), " (", percent(mean(is.na(indf$ldl)), 2), ")")
  )
})
table1a <- table1a %>%
  pivot_longer(cols = -eGFR) %>%
  pivot_wider(names_from = eGFR, values_from = value) %>%
  select(name, `90`, `75`, `60`, `45`, `30`, `15`) %>%
  separate(col = "name", into = c("name", "name2"), sep = ", ")
write_excel_csv(x = table1a, file = "data/20-table-1a.csv")

# Table with size by cluster (2A)
oo <- readRDS(file = "data/06-outcome-dfs.RDS")
table2a <- map_dfr(.x = seq_along(oo), .f = function(i) {
  with(
    oo[[i]],
    table(cluster, description)
  ) %>%
    as.data.frame() %>%
    filter(Freq > 0) %>%
    mutate(eGFR = names(oo)[i]) %>%
    mutate(cluster = as.character(cluster))
}) %>%
  select(eGFR, cluster, Freq, description) %>%
  rename(`eGFR category` = eGFR, Cluster = cluster, n = Freq, `Key condition(s)` = description) %>%
  arrange(`eGFR category`, Cluster) %>%
  mutate(n = paste0(comma(n, 0)))
table2a
write_excel_csv(x = table2a, file = "data/20-table-2a.csv")

# Follow-up
fwup <- map_dfr(.x = seq_along(oo), .f = function(i) {
  #
  f_acm <- survfit(Surv(t_acm, acm == 0) ~ 1, data = oo[[i]])
  f_acm <- tibble(
    median = paste0(comma(quantile(f_acm)$quantile["50"], 2)),
    lower = paste0(comma(quantile(f_acm)$lower["50"], 2)),
    upper = paste0(comma(quantile(f_acm)$upper["50"], 2)),
    outcome = "ACM"
  )
  #
  f_mace3 <- survfit(Surv(t_mace3, mace3 == 0) ~ 1, data = oo[[i]])
  f_mace3 <- tibble(
    median = paste0(comma(quantile(f_mace3)$quantile["50"], 2)),
    lower = paste0(comma(quantile(f_mace3)$lower["50"], 2)),
    upper = paste0(comma(quantile(f_mace3)$upper["50"], 2)),
    outcome = "MACE3"
  )
  #
  f_mace4 <- survfit(Surv(t_mace4, mace4 == 0) ~ 1, data = oo[[i]])
  f_mace4 <- tibble(
    median = paste0(comma(quantile(f_mace4)$quantile["50"], 2)),
    lower = paste0(comma(quantile(f_mace4)$lower["50"], 2)),
    upper = paste0(comma(quantile(f_mace4)$upper["50"], 2)),
    outcome = "MACE4"
  )
  out <- bind_rows(
    f_acm,
    f_mace3,
    f_mace4
  )
  out[["eGFR"]] <- names(oo)[i]
  return(out)
})
fwup
write_excel_csv(x = fwup, file = "data/20-table-fwup.csv")

# Adding example to describe eGFR interpolation procedure
model <- readRDS(file = "data/01-fitted-lmm.RDS")
beta <- fixef(model)
Sigma <- as.matrix(VarCorr(model)$lopnr)
residual_sigma <- attr(x = VarCorr(model), "sc")
set.seed(643132422)
bs <- rmvnorm(n = 2, sigma = Sigma)
fakedf <- bind_rows(
  tibble(
    id = 1,
    b0 = bs[1, 1],
    b1 = bs[1, 2],
    time = runif(n = 10, min = 0, max = 12)
  ),
  tibble(
    id = 2,
    b0 = bs[2, 1],
    b1 = bs[2, 2],
    time = runif(n = 10, min = 0, max = 12)
  )
) %>%
  arrange(id, time) %>%
  mutate(e = rnorm(n = nrow(.), sd = residual_sigma)) %>%
  mutate(m = (b0 + beta[1]) + (b1 + beta[2]) * time) %>%
  mutate(eGFR = m + e) %>%
  mutate(date = as.Date("2006-01-01") + time * 365.242) %>%
  mutate(id = factor(id, levels = 1:2, labels = c("ID = 1", "ID = 2")))

p_interpol <- ggplot(fakedf, aes(x = date, y = eGFR, color = id)) +
  geom_point(aes(shape = id)) +
  geom_line(aes(group = id, y = m), linetype = "dashed", show.legend = FALSE) +
  scale_y_continuous(breaks = seq(0, 150, by = 15)) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  scale_color_brewer(type = "qual", palette = "Dark2") +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(legend.position = "top") +
  labs(x = "Calendar Time", color = "", shape = "") +
  # ID = 1, eGFR = 75
  geom_segment(aes(color = "ID = 1", x = as.Date("2006-01-01"), xend = as.Date("2006-01-01") + (75 - beta[1] - bs[1, 1]) / (beta[2] + bs[1, 2]) * 365.242, y = 75, yend = 75), linetype = "dotted") +
  geom_segment(aes(color = "ID = 1", x = as.Date("2006-01-01") + (75 - beta[1] - bs[1, 1]) / (beta[2] + bs[1, 2]) * 365.242, xend = as.Date("2006-01-01") + (75 - beta[1] - bs[1, 1]) / (beta[2] + bs[1, 2]) * 365.242, y = 75, yend = min(eGFR)), linetype = "dotted") +
  geom_text(data = summarise(fakedf, min_eGFR = min(eGFR)), aes(color = "ID = 1", x = as.Date("2006-01-01") + (75 - beta[1] - bs[1, 1]) / (beta[2] + bs[1, 2]) * 365.242, y = min_eGFR, label = "Index Date (ID=1, eGFR=75)", family = "Arial Narrow"), show.legend = FALSE, vjust = 1) +
  # ID = 1, eGFR = 60
  geom_segment(aes(color = "ID = 1", x = as.Date("2006-01-01"), xend = as.Date("2006-01-01") + (60 - beta[1] - bs[1, 1]) / (beta[2] + bs[1, 2]) * 365.242, y = 60, yend = 60), linetype = "dotted") +
  geom_segment(aes(color = "ID = 1", x = as.Date("2006-01-01") + (60 - beta[1] - bs[1, 1]) / (beta[2] + bs[1, 2]) * 365.242, xend = as.Date("2006-01-01") + (60 - beta[1] - bs[1, 1]) / (beta[2] + bs[1, 2]) * 365.242, y = 60, yend = min(eGFR)), linetype = "dotted") +
  geom_text(data = summarise(fakedf, min_eGFR = min(eGFR)), aes(color = "ID = 1", x = as.Date("2006-01-01") + (60 - beta[1] - bs[1, 1]) / (beta[2] + bs[1, 2]) * 365.242, y = min_eGFR, label = "Index Date (ID=1, eGFR=60)", family = "Arial Narrow"), show.legend = FALSE, vjust = 1) +
  # ID = 2, eGFR = 90
  geom_segment(aes(color = "ID = 2", x = as.Date("2006-01-01"), xend = as.Date("2006-01-01") + (90 - beta[1] - bs[2, 1]) / (beta[2] + bs[2, 2]) * 365.242, y = 90, yend = 90), linetype = "dotted") +
  geom_segment(aes(color = "ID = 2", x = as.Date("2006-01-01") + (90 - beta[1] - bs[2, 1]) / (beta[2] + bs[2, 2]) * 365.242, xend = as.Date("2006-01-01") + (90 - beta[1] - bs[2, 1]) / (beta[2] + bs[2, 2]) * 365.242, y = 90, yend = min(eGFR)), linetype = "dotted") +
  geom_text(data = summarise(fakedf, min_eGFR = min(eGFR)), aes(color = "ID = 2", x = as.Date("2006-01-01") + (90 - beta[1] - bs[2, 1]) / (beta[2] + bs[2, 2]) * 365.242, y = min_eGFR, label = "Index Date (ID=2, eGFR=90)", family = "Arial Narrow"), show.legend = FALSE, vjust = 1)
#
ggsave(plot = p_interpol, filename = "figures/20-interpolation-example.pdf", width = 8, height = 4, device = cairo_pdf)
ggsave(plot = p_interpol, filename = "figures/20-interpolation-example.png", width = 8, height = 4, device = ragg::agg_png, res = 600, units = "in")
