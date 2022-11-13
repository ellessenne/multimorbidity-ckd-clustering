### Here we analyse association of clusters with MACE
library(tidyverse)
library(survival)
library(broom)
library(aod)
library(knitr)
library(stdReg)
library(cowplot)
library(readxl)

# Base data, used before
base_data <- readRDS("path/to/file")
# File with causes of death information
cod <- readRDS("path/to/file")
# Hospitalisation data
hosp <- readRDS("path/to/file")
mm <- readRDS("data/03-mm.RDS")
fitted_clusters <- readRDS("data/04-fitted-clusters.RDS")
iid <- readRDS("data/02-interpolated-index-dates.RDS")
n_to_pick <- c("7", "9", "7", "6", "6", "5")

# Fix dates
cod$dodsdat <- as.Date(cod$dodsdat, format = "%Y-%m-%d")
# 'dodsdat' = death date
hosp$bdat <- as.Date(hosp$bdat, format = "%Y-%m-%d")
# 'bdat' = date of each hospitalisation

# Define MACE3, MACE4, and all-cause mortality
ooo <- map(.x = seq_along(mm), .f = function(i) {
  tmp <- mm[[i]] %>%
    select(lopnr) %>%
    left_join(iid[[i]], by = "lopnr")
  tmp$cluster <- fitted_clusters[[i]][[n_to_pick[i]]]$cluster
  tmp <- tmp %>%
    left_join(select(base_data, lopnr, dodsdat, eof) %>% distinct(), by = "lopnr")
  cvd <- filter(cod, lopnr %in% tmp$lopnr) %>%
    filter(grepl(pattern = "^G45|^G46|^H341|^I", x = diagnos)) %>%
    select(lopnr, dodsdat) %>%
    mutate(cvd = 1) %>%
    rename(cvd_date = dodsdat)
  mi <- select(tmp, lopnr, index_date) %>%
    left_join(hosp, by = "lopnr") %>%
    filter(bdat >= index_date) %>%
    filter(grepl(pattern = "^I21|^I22|^I23", x = diagnosis)) %>%
    mutate(mi = 1) %>%
    rename(mi_date = bdat) %>%
    select(-diagnosis, -index_date) %>%
    group_by(lopnr) %>%
    filter(mi_date == min(mi_date)) %>%
    ungroup() %>%
    distinct()
  hf <- select(tmp, lopnr, index_date) %>%
    left_join(hosp, by = "lopnr") %>%
    filter(bdat >= index_date) %>%
    filter(grepl(pattern = "^I099|^I110|^I130|^I132|^I255|^I420|^I425|^I426|^I427|^I428|^I429|^I43|^I50", x = diagnosis)) %>%
    mutate(hf = 1) %>%
    rename(hf_date = bdat) %>%
    select(-diagnosis, -index_date) %>%
    group_by(lopnr) %>%
    filter(hf_date == min(hf_date)) %>%
    ungroup() %>%
    distinct()
  stroke <- select(tmp, lopnr, index_date) %>%
    left_join(hosp, by = "lopnr") %>%
    filter(bdat >= index_date) %>%
    filter(grepl(pattern = "^H341|^G45|^G46|^I60|^I61|^I63|^I64", x = diagnosis)) %>%
    mutate(stroke = 1) %>%
    rename(stroke_date = bdat) %>%
    select(-diagnosis, -index_date) %>%
    group_by(lopnr) %>%
    filter(stroke_date == min(stroke_date)) %>%
    ungroup() %>%
    distinct()
  mace <- cvd %>%
    full_join(mi, by = "lopnr") %>%
    full_join(hf, by = "lopnr") %>%
    full_join(stroke, by = "lopnr")
  tmp <- left_join(tmp, mace, "lopnr") %>%
    replace_na(list(cvd = 0, mi = 0, hf = 0, stroke = 0)) %>%
    mutate(
      mace4 = pmax(cvd, mi, hf, stroke),
      mace4_date = pmin(cvd_date, mi_date, hf_date, stroke_date, eof, na.rm = TRUE)
    ) %>%
    mutate(
      mace3 = pmax(cvd, mi, stroke),
      mace3_date = pmin(cvd_date, mi_date, stroke_date, eof, na.rm = TRUE)
    ) %>%
    mutate(
      acm = ifelse(dodsdat <= eof & !is.na(dodsdat), 1, 0),
      acm_date = pmin(dodsdat, eof, na.rm = TRUE)
    ) %>%
    select(lopnr, index_date, cluster, mace3, mace3_date, mace4, mace4_date, acm, acm_date) %>%
    left_join(distinct(base_data, lopnr, female, dob), by = "lopnr") %>%
    mutate(age = as.numeric(index_date - dob) / 365.252) %>%
    mutate(age5 = age / 5) %>%
    select(-dob) %>%
    mutate(
      t_mace4 = as.numeric(mace4_date - index_date) / 365.242,
      t_mace3 = as.numeric(mace3_date - index_date) / 365.242,
      t_acm = as.numeric(acm_date - index_date) / 365.242
    )
  return(tmp)
})
names(ooo) <- names(mm)

# Describe clusters
# This was subjectively defined by the analysts in collaboration with the clinical experts
cdesc <- read_excel("data/05-SCREAM-cluster-descriptions-v5.xlsx") %>%
  mutate(
    description = tools::toTitleCase(description),
    system = tools::toTitleCase(system),
    eGFR = as.character(eGFR)
  )
ooo <- map(.x = seq_along(ooo), .f = function(i) {
  ooo[[i]] %>%
    mutate(eGFR = names(ooo)[i]) %>%
    left_join(cdesc, by = c("eGFR", "cluster")) %>%
    mutate(
      description = factor(description),
      system = factor(system),
      cluster = factor(cluster)
    ) %>%
    mutate(description = relevel(description, ref = "No Prominent Condition")) %>%
    mutate(cluster = relevel(cluster, ref = cdesc$cluster[cdesc$ref == 1 & cdesc$eGFR == names(ooo)[i]])) %>%
    select(-ref)
})
names(ooo) <- names(mm)
saveRDS(object = ooo, file = "data/06-outcome-dfs.RDS")

# Fit some models...
fits_mace4 <- map(.x = ooo, .f = function(x) {
  f <- coxph(Surv(t_mace4, mace4 == 1) ~ cluster + age5 + factor(female), data = x, x = FALSE, y = FALSE, ties = "breslow")
  f$tvar <- "t_mace4"
  f$dvar <- "mace4"
  f$range_times <- range(x[[f$tvar]][x[[f$dvar]] == 1])
  f$avg_age <- mean(x$age)
  f$avg_female <- mean(x$female)
  f
})
names(fits_mace4) <- names(ooo)
saveRDS(object = fits_mace4, file = "data/06-fits-mace4.RDS")

fits_mace3 <- map(.x = ooo, .f = function(x) {
  f <- coxph(Surv(t_mace3, mace3 == 1) ~ cluster + age5 + factor(female), data = x, x = FALSE, y = FALSE, ties = "breslow")
  f$tvar <- "t_mace3"
  f$dvar <- "mace3"
  f$range_times <- range(x[[f$tvar]][x[[f$dvar]] == 1])
  f$avg_age <- mean(x$age)
  f$avg_female <- mean(x$female)
  f
})
names(fits_mace3) <- names(ooo)
saveRDS(object = fits_mace3, file = "data/06-fits-mace3.RDS")

fits_acm <- map(.x = ooo, .f = function(x) {
  f <- coxph(Surv(t_acm, acm == 1) ~ cluster + age5 + factor(female), data = x, x = FALSE, y = FALSE, ties = "breslow")
  f$tvar <- "t_acm"
  f$dvar <- "acm"
  f$range_times <- range(x[[f$tvar]][x[[f$dvar]] == 1])
  f$avg_age <- mean(x$age)
  f$avg_female <- mean(x$female)
  f
})
names(fits_acm) <- names(ooo)
saveRDS(object = fits_acm, file = "data/06-fits-acm.RDS")

# Summarise models with forest plot
fit_summ <- map_dfr(.x = seq_along(ooo), .f = function(i) {
  d1 <- tidy(fits_mace3[[i]], exponentiate = TRUE, conf.int = TRUE)
  opts <- paste0("cluster", seq(max(cdesc$cluster[cdesc$eGFR == names(ooo)[i]])))
  opts <- opts[!(opts %in% d1$term[grepl(pattern = "cluster", x = d1$term)])]
  .to_add <- tibble(term = c(opts, "factor(female)0"), estimate = 1, std.error = 0)
  d1 <- bind_rows(d1, .to_add)
  d1$outcome <- "MACE3"
  d2 <- tidy(fits_mace4[[i]], exponentiate = TRUE, conf.int = TRUE)
  d2 <- bind_rows(d2, .to_add)
  d2$outcome <- "MACE4"
  d3 <- tidy(fits_acm[[i]], exponentiate = TRUE, conf.int = TRUE)
  d3 <- bind_rows(d3, .to_add)
  d3$outcome <- "ACM"
  d <- bind_rows(d1, d2, d3)
  d$eGFR <- names(ooo)[i]
  lvls <- c("age5", "factor(female)0", "factor(female)1", "cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6", "cluster7", "cluster8", "cluster9")
  lbls <- c("Age (Per 5 Years)", "Sex = Male (ref)", "Sex = Female", "Cluster: 1", "Cluster: 2", "Cluster: 3", "Cluster: 4", "Cluster: 5", "Cluster: 6", "Cluster: 7", "Cluster: 8", "Cluster: 9")
  lbls[lvls == opts] <- paste(lbls[lvls == opts], "(ref)")
  d$term <- as.character(factor(d$term, levels = lvls, labels = lbls))
  d$cluster <- as.numeric(substr(d$term, 10, 11))
  d$cluster[!grepl(pattern = "Cluster", x = d$term)] <- NA
  d
})
fit_summ <- mutate(fit_summ,
  term = factor(term),
  outcome = factor(outcome, levels = c("ACM", "MACE3", "MACE4")),
  eGFR = factor(eGFR)
) %>%
  rename(Outcome = outcome) %>%
  mutate(sig = p.value <= 0.05) %>%
  replace_na(list(sig = 0)) %>%
  arrange(eGFR, Outcome, term) %>%
  left_join(cdesc, by = c("eGFR", "cluster")) %>%
  mutate(description = ifelse(Outcome == "MACE4", description, NA))
p_hr <- ggplot(fit_summ) +
  geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
  geom_errorbar(aes(x = fct_rev(term), ymin = conf.low, ymax = conf.high, color = factor(sig)), width = 1 / 3, show.legend = FALSE) +
  geom_point(aes(x = fct_rev(term), y = estimate, color = factor(sig)), show.legend = FALSE) +
  facet_grid(eGFR ~ Outcome, labeller = label_both, scales = "free_y") +
  scale_color_manual(values = c("grey75", "black")) +
  coord_flip() +
  labs(y = "Hazard Ratio (95% C.I.)", x = "") +
  theme_minimal(base_size = 12, base_family = "Arial Narrow")
ggsave(plot = p_hr, filename = "figures/06-forest-plot.pdf", width = 6, height = 6 * sqrt(2), device = cairo_pdf)
ggsave(plot = p_hr, filename = "figures/06-forest-plot.png", width = 6, height = 6 * sqrt(2), device = ragg::agg_png, res = 600, units = "in")

# Same as above, but with cluster descriptions on the side
p_hr_desc <- ggplot(fit_summ) +
  geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
  geom_errorbar(aes(x = fct_rev(term), ymin = conf.low, ymax = conf.high, color = factor(sig)), width = 1 / 3, show.legend = FALSE) +
  geom_point(aes(x = fct_rev(term), y = estimate, color = factor(sig)), show.legend = FALSE) +
  geom_text(aes(label = description, x = term, y = 6.5), hjust = 0, size = 3) +
  facet_grid(eGFR ~ Outcome, labeller = label_both, scales = "free_y") +
  scale_color_manual(values = c("grey75", "black")) +
  coord_flip(clip = FALSE, ylim = c(0, 5)) +
  labs(y = "Hazard Ratio (95% C.I.)", x = "") +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(plot.margin = unit(c(0.1, 15, 0.1, 0.1), "lines"))
ggsave(plot = p_hr_desc, filename = "figures/06-forest-plot-description.pdf", width = 9, height = 8, device = cairo_pdf)
ggsave(plot = p_hr_desc, filename = "figures/06-forest-plot-description.png", width = 9, height = 8, device = ragg::agg_png, res = 600, units = "in")

# Subset of forest plot
p_hr_subset <- fit_summ %>%
  filter(Outcome %in% c("ACM", "MACE4")) %>%
  filter(eGFR %in% c("30", "90")) %>%
  ggplot(aes(x = fct_rev(term), y = estimate, color = factor(sig))) +
  geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 1 / 3, show.legend = FALSE) +
  geom_point(show.legend = FALSE) +
  facet_grid(eGFR ~ Outcome, labeller = label_both, scales = "free_y") +
  scale_color_manual(values = c("grey75", "black")) +
  coord_flip(ylim = c(0, 4)) +
  labs(y = "Hazard Ratio (95% C.I.)", x = "") +
  theme_minimal(base_size = 12, base_family = "Arial Narrow")
p_hr_subset <- plot_grid(p_hr_subset, labels = "A") # Replace with B for SAIL
ggsave(plot = p_hr_subset, filename = "figures/06-forest-plot-subset.pdf", width = 6, height = 3, device = cairo_pdf)
ggsave(plot = p_hr_subset, filename = "figures/06-forest-plot-subset.png", width = 6, height = 3, device = ragg::agg_png, res = 600, units = "in")

# Subset of forest plot, but with cluster descriptions on the side
p_hr_subset_desc <- fit_summ %>%
  filter(Outcome %in% c("ACM", "MACE4")) %>%
  filter(eGFR %in% c("30", "90")) %>%
  ggplot() +
  geom_hline(yintercept = 1, color = "red", linetype = "dotted") +
  geom_errorbar(aes(x = fct_rev(term), ymin = conf.low, ymax = conf.high, color = factor(sig)), width = 1 / 3, show.legend = FALSE) +
  geom_point(aes(x = fct_rev(term), y = estimate, color = factor(sig)), show.legend = FALSE) +
  geom_text(aes(label = description, x = term, y = 6), hjust = 0, size = 3) +
  facet_grid(eGFR ~ Outcome, labeller = label_both, scales = "free_y") +
  scale_color_manual(values = c("grey75", "black")) +
  coord_flip(clip = FALSE, ylim = c(0, 5)) +
  labs(y = "Hazard Ratio (95% C.I.)", x = "") +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(plot.margin = unit(c(0.1, 15, 0.1, 0.1), "lines"))
ggsave(plot = p_hr_subset_desc, filename = "figures/06-forest-plot-subset-description.pdf", width = 9, height = 4, device = cairo_pdf)
ggsave(plot = p_hr_subset_desc, filename = "figures/06-forest-plot-subset-description.png", width = 9, height = 4, device = ragg::agg_png, res = 600, units = "in")

# Then, export the hazard ratios in table form too
sink(file = "data/06-tables-hr-by-outcome.md")
for (e in names(ooo)) {
  cat("eGFR:", e)
  kable(
    x = fit_summ %>%
      mutate(y = case_when(
        is.na(conf.low) ~ "---",
        TRUE ~ paste0(formattable::comma(estimate, 3), " (", formattable::comma(conf.low, 3), "-", formattable::comma(conf.high, 3), ")")
      )) %>%
      select(Outcome, eGFR, term, y) %>%
      pivot_wider(names_from = "Outcome", values_from = "y") %>%
      filter(eGFR == e) %>%
      select(-eGFR) %>%
      rename(Term = term),
    format = "markdown"
  ) %>%
    print()
  cat("\n\n")
}
sink()

# Wald tests for the cluster effect
sink(file = "data/06-tables-wald-tests.md")
for (e in names(ooo)) {
  cat("eGFR:", e)
  x <- fits_acm[[e]]
  wt <- wald.test(Sigma = vcov(x), b = coef(x), Terms = grep("cluster", names(coef(x))))
  cat("\nOutcome ACM: Chi2 =", wt$result$chi2["chi2"], "- df =", wt$result$chi2["df"], "- p.value =", wt$result$chi2["P"])
  x <- fits_mace3[[e]]
  wt <- wald.test(Sigma = vcov(x), b = coef(x), Terms = grep("cluster", names(coef(x))))
  cat("\nOutcome MACE3: Chi2 =", wt$result$chi2["chi2"], "- df =", wt$result$chi2["df"], "- p.value =", wt$result$chi2["P"])
  x <- fits_mace4[[e]]
  wt <- wald.test(Sigma = vcov(x), b = coef(x), Terms = grep("cluster", names(coef(x))))
  cat("\nOutcome MACE4: Chi2 =", wt$result$chi2["chi2"], "- df =", wt$result$chi2["df"], "- p.value =", wt$result$chi2["P"])
  cat("\n\n")
}
sink()

# Finally, standardised survival curves
get_std_curves <- function(x, .outcome, .progress = TRUE) {
  if (.progress) pb <- txtProgressBar(min = 0, max = length(x), style = 3)
  stdsurv <- map(.x = seq_along(x), .f = function(i) {
    stdf <- stdReg::stdCoxph(fit = x[[i]], data = ooo[[i]], X = "cluster", t = seq(x[[i]]$range_times[1], x[[i]]$range_times[2], length.out = 100))
    sumstdf <- summary(stdf)
    pdf <- map_dfr(.x = seq_along(sumstdf$est.table), .f = function(j) {
      sumstdf$est.table[[j]] %>%
        as_tibble() %>%
        mutate(t = sumstdf$input$t[j]) %>%
        mutate(cluster = rownames(sumstdf$est.table[[j]]))
    })
    names(pdf) <- make.names(names(pdf))
    pdf[["eGFR"]] <- names(ooo)[i]
    pdf[["Outcome"]] <- .outcome
    if (.progress) setTxtProgressBar(pb = pb, value = i)
    return(pdf)
  })
  if (.progress) close(pb)
  return(stdsurv)
}
stdsurv_mace4 <- get_std_curves(x = fits_mace4, .outcome = "MACE4")
stdsurv_mace3 <- get_std_curves(x = fits_mace3, .outcome = "MACE3")
stdsurv_acm <- get_std_curves(x = fits_acm, .outcome = "ACM")

stdsurv <- bind_rows(stdsurv_acm, stdsurv_mace3, stdsurv_mace4)
for (e in unique(stdsurv$eGFR)) {
  p_df <- filter(stdsurv, eGFR == e)
  p_surv <- ggplot(p_df, aes(x = t, y = Estimate)) +
    pammtools::geom_stepribbon(aes(ymin = lower.0.95, ymax = upper.0.95, fill = cluster), alpha = 1 / 3) +
    geom_step(aes(color = cluster)) +
    facet_wrap(~Outcome, labeller = label_both) +
    scale_x_continuous(breaks = seq(0, 10, by = 2)) +
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    coord_cartesian(ylim = c(0, 1), xlim = c(0, 10)) +
    labs(x = "Time", y = "Standardised Survival Probability", color = "Cluster", fill = "Cluster", title = paste0("Outcomes for eGFR = ", e)) +
    theme_minimal(base_size = 12, base_family = "Arial Narrow") +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))
  ggsave(plot = p_surv, filename = paste0("figures/06-survival-curves-plot-eGFR=", e, ".pdf"), width = 7, height = 7 / sqrt(3), device = cairo_pdf)
  ggsave(plot = p_surv, filename = paste0("figures/06-survival-curves-plot-eGFR=", e, ".png"), width = 7, height = 7 / sqrt(3), device = ragg::agg_png, res = 600, units = "in")
}

# Calculate crude event rates
crude_rates <- map_dfr(.x = seq_along(ooo), .f = function(i) {
  tmp <- ooo[[i]] %>%
    group_by(cluster) %>%
    summarise(
      C_mace3 = sum(t_mace3),
      n_mace3 = sum(mace3),
      C_mace4 = sum(t_mace4),
      n_mace4 = sum(mace4),
      C_acm = sum(t_acm),
      n_acm = sum(acm),
      N = n()
    )
  bind_rows(
    mutate(tmp,
      rate = n_mace3 / C_mace3,
      lower = (1 / C_mace3) * (n_mace3 - qnorm(1 - 0.05 / 2) * sqrt(n_mace3)),
      upper = (1 / C_mace3) * (n_mace3 + qnorm(1 - 0.05 / 2) * sqrt(n_mace3)),
      Outcome = "MACE3"
    ),
    mutate(tmp,
      rate = n_mace4 / C_mace4,
      lower = (1 / C_mace4) * (n_mace4 - qnorm(1 - 0.05 / 2) * sqrt(n_mace4)),
      upper = (1 / C_mace4) * (n_mace4 + qnorm(1 - 0.05 / 2) * sqrt(n_mace4)),
      Outcome = "MACE4"
    ),
    mutate(tmp,
      rate = n_acm / C_acm,
      lower = (1 / C_acm) * (n_acm - qnorm(1 - 0.05 / 2) * sqrt(n_acm)),
      upper = (1 / C_acm) * (n_acm + qnorm(1 - 0.05 / 2) * sqrt(n_acm)),
      Outcome = "ACM"
    )
  ) %>%
    mutate(eGFR = names(ooo)[i]) %>%
    select(cluster, eGFR, Outcome, rate, lower, upper)
}) %>%
  mutate(cluster = paste("Cluster:", cluster)) %>%
  mutate(cluster = factor(cluster)) %>%
  mutate(cluster = fct_rev(cluster))
p_rates <- ggplot(crude_rates, aes(x = cluster, y = rate * 1000)) +
  geom_errorbar(aes(ymin = lower * 1000, ymax = upper * 1000), width = 1 / 3) +
  geom_point() +
  facet_grid(eGFR ~ Outcome, scales = "free_y", labeller = label_both) +
  coord_flip() +
  labs(x = "", y = "Crude Event Rates (Per 1,000 Patient-Years)") +
  theme_minimal(base_size = 12, base_family = "Arial Narrow")
ggsave(plot = p_rates, filename = "figures/06-crude-event-rates.pdf", width = 7, height = 7 / sqrt(2), device = cairo_pdf)
ggsave(plot = p_rates, filename = "figures/06-crude-event-rates.png", width = 7, height = 7 / sqrt(2), device = ragg::agg_png, res = 600, units = "in")
