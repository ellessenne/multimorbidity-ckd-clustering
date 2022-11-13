### Script calculate multimorbidity conditions at each index date
library(tidyverse)
library(SCREAM)
library(matrixStats)
library(knitr)
library(hrbrthemes)

# 'data' here is the input dataset
data <- readRDS(file = "path/to/file")

# Import splits
splits <- readRDS(file = "data/02-interpolated-index-dates.RDS")

# Now each element of 'splits' is a dataframe with index dates
# Time to augment this with ICD codes
# SLV ->
# (SLV is hospital care)
slv <- readRDS(file = "path/to/file")
slv$datum <- as.Date(slv$bdat, format = "%Y-%m-%d")
slv$bdat <- NULL
slv_data <- map(.x = splits, .f = function(x) {
  out <- left_join(x, slv, by = "lopnr") %>%
    filter(datum <= index_date) %>%
    arrange(lopnr, datum)
  return(out)
})
rm(slv)
gc()
# 'datum' is date at which each code is recorded
# 'diagnosis' will be the column containing each code
# 'index_date' comes from the interpolation procedure

# OVR + KON ->
# (OVR is outpatient care, KON is GP care)
ovr <- readRDS(file = "path/to/file")
kon <- readRDS(file = "path/to/file")
ovr_kon <- bind_rows(ovr, kon) %>%
  distinct()
rm(ovr, kon)
gc()
ovr_kon$datum <- as.Date(ovr_kon$bdat, format = "%Y-%m-%d")
ovr_kon$bdat <- NULL
ovr_kon_data <- map(.x = splits, .f = function(x) {
  out <- left_join(x, ovr_kon, by = "lopnr") %>%
    filter(datum <= index_date) %>%
    arrange(lopnr, datum)
  return(out)
})
rm(ovr_kon)
gc()

# Drugs/medications data ->
drugs <- readRDS(file = "path/to/file")
ids <- map(.x = splits, .f = function(x) x$lopnr)
ids <- unlist(ids)
names(ids) <- NULL
ids <- unique(ids)
drugs <- filter(drugs, lopnr %in% ids)
drugs_data <- map(.x = splits, .f = function(x) {
  out <- left_join(x, drugs, by = "lopnr") %>%
    filter(edatum <= index_date) %>%
    arrange(lopnr, edatum) %>%
    group_by(lopnr, index_date, edatum, atc) %>%
    summarise(antal = sum(antal)) %>%
    ungroup() %>%
    rename(datum = edatum)
  return(out)
})
rm(drugs)
gc()
# - 'datum' here will be each dispensation date
# - 'antal' will be the number of packages dispensed
# - 'atc' will be the ATC code of each medication

# Now, we calculate the multimorbidity conditions
# We combine all cancers here
mm <- map(.x = seq_along(splits), .f = function(i) {
  d <- SCREAM::multimorbidity(data_hospitalisations = slv_data[[i]], data_claims = ovr_kon_data[[i]], data_drugs = drugs_data[[i]], id = "lopnr", code = "diagnosis", atc = "atc", npacks = "antal", date = "datum", index_date = "index_date")
  d[["cancer"]] <- pmax(d[["cancer_lymphoma"]], d[["cancer_metastatic"]], d[["cancer_nonmetastatic"]])
  d[["cancer_lymphoma"]] <- NULL
  d[["cancer_metastatic"]] <- NULL
  d[["cancer_nonmetastatic"]] <- NULL
  usethis::ui_done(x = "Done with eGFR = {names(splits)[i]}")
  d
})
# Then, restore subjects with no codes whatsoever
mm <- map(.x = seq_along(mm), .f = function(i) {
  out <- left_join(select(splits[[i]], lopnr), mm[[i]], by = "lopnr")
  out[is.na(out)] <- 0
  return(out)
})
names(mm) <- names(splits)
saveRDS(object = mm, file = "data/03-mm.RDS")

# Summarising that as a test
prop_conditions <- map_dfr(.x = seq_along(mm), .f = function(i) {
  mm[[i]] %>%
    dplyr::select(-lopnr) %>%
    summarise(across(.fns = mean)) %>%
    mutate(eGFR = names(splits)[i]) %>%
    pivot_longer(cols = -eGFR)
}) %>%
  mutate(label = paste0(formattable::percent(value, 1))) %>%
  mutate(name = factor(name, levels = c("alcohol_misuse", "asthma", "afib", "cancer", "ckd", "cpain", "cvhepatitis", "cirrhosis", "dementia", "depression", "diabetes", "epilepsy", "chf", "hypertension", "hypothyroidism", "ibd", "ibs", "multiple_sclerosis", "mi", "parkinson", "cpd", "pud", "pvd", "psoriasis", "rheum_arthritis", "schizofrenia", "severe_constipation", "stroke"), labels = c("Alcohol Misuse", "Asthma", "Atrial Fibrillation", "Cancer", "CKD", "Chronic Pain", "Chronic Viral Hepatitis B", "Cirrhosis", "Dementia", "Depression", "Diabetes", "Epilepsy", "Heart Failure", "Hypertension", "Hypothyroidism", "IBD", "IBS", "Multiple Sclerosis", "MI", "Parkinson's Disease", "Pulmonary Disease", "PUD", "PVD", "Psoriasis", "Rheumatoid Arthritis", "Schizophrenia", "Severe Constipation", "Stroke")))

p <- ggplot(prop_conditions, aes(x = eGFR, y = value)) +
  geom_hline(yintercept = 1, color = "grey66", linetype = "dashed") +
  geom_col() +
  geom_text(aes(label = label), angle = 90, hjust = -0.1, vjust = 0.5, size = 2.5) +
  facet_wrap(~name) +
  scale_y_sqrt(breaks = c(0, 0.1, 0.25, 0.5, 1)) +
  coord_cartesian(ylim = c(0, 2)) +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(y = "Proportion")
ggsave(plot = p, filename = "figures/03-conditions-plot.pdf", height = 7, width = 7 * 29.7 / 21, device = cairo_pdf)
ggsave(plot = p, filename = "figures/03-conditions-plot.png", height = 7, width = 7 * 29.7 / 21, device = ragg::agg_png, res = 600, units = "in")

p_no_ckd <- ggplot(filter(prop_conditions, name != "CKD"), aes(x = eGFR, y = value)) +
  geom_hline(yintercept = 1, color = "grey66", linetype = "dashed") +
  geom_col() +
  geom_text(aes(label = label), angle = 90, hjust = -0.1, vjust = 0.5, size = 2.5) +
  facet_wrap(~name) +
  scale_y_sqrt(breaks = c(0, 0.1, 0.25, 0.5, 1)) +
  coord_cartesian(ylim = c(0, 2)) +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(y = "Proportion")
ggsave(plot = p_no_ckd, filename = "figures/03-conditions-plot-no-ckd.pdf", height = 7, width = 7 * 29.7 / 21, device = cairo_pdf)
ggsave(plot = p_no_ckd, filename = "figures/03-conditions-plot-no-ckd.png", height = 7, width = 7 * 29.7 / 21, device = ragg::agg_png, res = 600, units = "in")

# Making a stacked barplot with number of conditions
sbp <- map_dfr(.x = seq_along(mm), .f = function(i) {
  df <- tibble(
    lopnr = mm[[i]]$lopnr,
    mm = mm[[i]] %>%
      select(-lopnr, -ckd) %>%
      as.matrix() %>%
      rowSums2()
  ) %>%
    mutate(mmc = case_when(
      mm >= 8 ~ "8+",
      TRUE ~ as.character(mm)
    )) %>%
    mutate(mmc = factor(mmc, levels = c("0", "1", "2", "3", "4", "5", "6", "7", "8+")))
  table(df$mmc) %>%
    as.data.frame() %>%
    mutate(eGFR = names(splits)[i])
})
# Stacked proportions:
sbpp <- ggplot(sbp, aes(x = eGFR, y = Freq, fill = Var1)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_viridis_d() +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  labs(y = "Proportion", fill = "", caption = "Proportion of Morbidities by eGFR")
ggsave(plot = sbpp, filename = "figures/03-pmorbidities-plot.pdf", height = 3, width = 3 * sqrt(2), device = cairo_pdf)
ggsave(plot = sbpp, filename = "figures/03-pmorbidities-plot.png", height = 3, width = 3 * sqrt(2), device = ragg::agg_png, res = 600, units = "in")
# Stacked numbers
sbpn <- ggplot(sbp, aes(x = eGFR, y = Freq, fill = Var1)) +
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_viridis_d() +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  labs(y = "Count", fill = "", caption = "Number of Morbidities by eGFR")
ggsave(plot = sbpn, filename = "figures/03-nmorbidities-plot.pdf", height = 3, width = 3 * sqrt(2), device = cairo_pdf)
ggsave(plot = sbpn, filename = "figures/03-nmorbidities-plot.png", height = 3, width = 3 * sqrt(2), device = ragg::agg_png, res = 600, units = "in")
# Tables
sink(file = "data/03-nmorbidities-by-egfr.md")
sbp %>%
  pivot_wider(names_from = eGFR, values_from = Freq) %>%
  rename(`N. Morbidities` = Var1) %>%
  kable(format = "markdown", caption = "N. of morbidities by eGFR") %>%
  print()
sink()
sink(file = "data/03-pmorbidities-by-egfr.md")
sbp %>%
  group_by(eGFR) %>%
  mutate(Freq = Freq / sum(Freq)) %>%
  pivot_wider(names_from = eGFR, values_from = Freq) %>%
  rename(`Prop. Morbidities` = Var1) %>%
  kable(format = "markdown", caption = "Prop. of morbidities by eGFR") %>%
  print()
sink()

# Line-plot data (and single plot for SCREAM)
lpdf <- map_dfr(.x = seq_along(mm), .f = function(i) {
  .in <- tibble(
    lopnr = mm[[i]]$lopnr,
    mm = mm[[i]] %>%
      select(-lopnr, -ckd) %>%
      as.matrix() %>%
      rowSums2()
  )
  tibble(
    min = min(.in$mm),
    q1 = fivenum(.in$mm)[2],
    median = fivenum(.in$mm)[3],
    q3 = fivenum(.in$mm)[4],
    max = max(.in$mm),
    mean = mean(.in$mm),
    var = var(.in$mm),
    eGFR = as.numeric(names(mm)[i])
  )
})
readr::write_csv(x = lpdf, file = "data/03-line-plot-conditions-data.csv")
plpdf <- ggplot(data = lpdf, aes(x = eGFR, y = median)) +
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 1 / 3) +
  geom_errorbar(aes(ymin = q1, ymax = q3), width = 2 / 3) +
  geom_line(linetype = "dashed") +
  geom_point() +
  scale_y_continuous(breaks = 0:16) +
  scale_x_reverse(breaks = as.numeric(names(mm))) +
  coord_cartesian(ylim = c(0, 16)) +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(panel.grid.minor = element_blank()) +
  labs(x = expression(eGFR ~ (mL / min / 1.73 ~ m^2)), y = "Median Conditions (with IQI)")
ggsave(plot = plpdf, filename = "figures/03-line-plot-conditions.pdf", height = 3, width = 3 * sqrt(2), device = cairo_pdf)
ggsave(plot = plpdf, filename = "figures/03-line-plot-conditions.png", height = 3, width = 3 * sqrt(2), device = ragg::agg_png, res = 600, units = "in")

# Same as above, but by age
ipdata <- readRDS(file = "data/02-interpolated-index-dates.RDS")
lpdf <- map_dfr(.x = seq_along(mm), .f = function(i) {
  tibble(
    lopnr = mm[[i]]$lopnr,
    mm = mm[[i]] %>%
      select(-lopnr, -ckd) %>%
      as.matrix() %>%
      rowSums2()
  ) %>%
    left_join(ipdata[[i]], by = "lopnr") %>%
    left_join(distinct(data, lopnr, dob), by = "lopnr") %>%
    mutate(age = as.numeric(index_date - dob) / 365.242) %>%
    mutate(agec = ifelse(age <= 65, 0, 1)) %>%
    mutate(agec = factor(agec, levels = c(0, 1), labels = c("Age <= 65", "Age > 65"))) %>%
    group_by(agec) %>%
    summarise(
      age = mean(age),
      min = min(mm),
      q1 = fivenum(mm)[2],
      median = fivenum(mm)[3],
      q3 = fivenum(mm)[4],
      max = max(mm),
      mean = mean(mm),
      var = var(mm)
    ) %>%
    mutate(eGFR = as.numeric(names(mm)[i]))
})
readr::write_csv(x = lpdf, file = "data/03-line-plot-conditions-data-agec.csv")
plpdf <- ggplot(data = lpdf, aes(x = eGFR, y = median)) +
  geom_ribbon(aes(ymin = min, ymax = max, fill = agec), alpha = 1 / 5, position = position_dodge(width = 3)) +
  geom_errorbar(aes(ymin = q1, ymax = q3, color = agec), position = position_dodge(width = 3), width = 2) +
  geom_line(aes(color = agec), linetype = "dashed", position = position_dodge(width = 3)) +
  geom_point(aes(color = agec), position = position_dodge(width = 3)) +
  scale_y_continuous(breaks = 0:16) +
  scale_x_reverse(breaks = as.numeric(names(mm))) +
  scale_color_ipsum() +
  scale_fill_ipsum() +
  coord_cartesian(ylim = c(0, 16)) +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(panel.grid.minor = element_blank(), legend.position = c(1, 1), legend.justification = c(1, 1)) +
  labs(x = expression(eGFR ~ (mL / min / 1.73 ~ m^2)), y = "Median Conditions (with IQI)", color = "", fill = "")
plpdf
ggsave(plot = plpdf, filename = "figures/03-line-plot-conditions-agec.pdf", height = 3, width = 3 * sqrt(2), device = cairo_pdf)
ggsave(plot = plpdf, filename = "figures/03-line-plot-conditions-agec.png", height = 3, width = 3 * sqrt(2), device = ragg::agg_png, res = 600, units = "in")
