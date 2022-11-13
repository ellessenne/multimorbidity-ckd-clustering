library(tidyverse)
library(matrixStats)
library(survival)
library(glue)
library(riskRegression)

### Sensitivity analysis #1:
### Checking that predictive metrics for model with clusters vs
### models with just multimorbidity counts

# Assemble data
ooo <- readRDS(file = "data/06-outcome-dfs.RDS")
mm <- readRDS(file = "data/03-mm.RDS")
for (i in seq_along(mm)) {
  mm[[i]][["mmcount"]] <- rowSums2(x = as.matrix(select(mm[[i]], -lopnr)))
  mm[[i]] <- select(mm[[i]], lopnr, mmcount)
  ooo[[i]] <- left_join(ooo[[i]], mm[[i]], by = "lopnr")
}
rm(mm)
gc()

# Define function to fit either model to a given fold
fit_metrics <- function(data, outcome, egfr, .bootstrap = TRUE, .B = 100) {
  # Smaller subset of data
  # Should reduce memory usage
  data <- dplyr::select(data, -lopnr, -index_date, -mace3_date, -mace4_date, -acm_date, -age, -description, -system)
  if (!(outcome %in% c("ACM", "MACE3", "MACE4"))) stop("Selected 'outcome' not supported, possible choices are 'ACM', 'MACE3', 'MACE4'.", call. = FALSE)
  # The following is not very elegant, but the Score function below had problems picking up the correct formula
  # when defining that outside the coxph function.
  if (outcome == "ACM") {
    td_var <- data[["t_acm"]][data[["acm"]] == 1]
    fit0 <- coxph(Surv(t_acm, acm) ~ age5 + female, data = data, x = TRUE, ties = "breslow")
    fit1 <- coxph(Surv(t_acm, acm) ~ cluster + age5 + female, data = data, x = TRUE, ties = "breslow")
    fit2 <- coxph(Surv(t_acm, acm) ~ mmcount + age5 + female, data = data, x = TRUE, ties = "breslow")
  } else if (outcome == "MACE3") {
    td_var <- data[["t_mace3"]][data[["mace3"]] == 1]
    fit0 <- coxph(Surv(t_mace3, mace3) ~ age5 + female, data = data, x = TRUE, ties = "breslow")
    fit1 <- coxph(Surv(t_mace3, mace3) ~ cluster + age5 + female, data = data, x = TRUE, ties = "breslow")
    fit2 <- coxph(Surv(t_mace3, mace3) ~ mmcount + age5 + female, data = data, x = TRUE, ties = "breslow")
  } else if (outcome == "MACE4") {
    td_var <- data[["t_mace4"]][data[["mace4"]] == 1]
    fit0 <- coxph(Surv(t_mace4, mace4) ~ age5 + female, data = data, x = TRUE, ties = "breslow")
    fit1 <- coxph(Surv(t_mace4, mace4) ~ cluster + age5 + female, data = data, x = TRUE, ties = "breslow")
    fit2 <- coxph(Surv(t_mace4, mace4) ~ mmcount + age5 + female, data = data, x = TRUE, ties = "breslow")
  }
  ### ---
  td <- quantile(x = td_var, probs = seq(0.1, 0.9, by = 0.1))
  mu <- riskRegression::Score(
    object = list(
      "Age + Gender" = fit0,
      "Age + Gender + Cluster" = fit1,
      "Age + Gender + Count" = fit2
    ),
    formula = update(fit0$formula, . ~ 1),
    se.fit = 0,
    times = td,
    data = data,
    contrasts = FALSE,
    progress.bar = NULL
  )
  mu_out <- bind_rows(
    mu$AUC$score %>% mutate(metric = "AUC") %>% rename(value = AUC),
    mu$Brier$score %>% mutate(metric = "Brier") %>% rename(value = Brier)
  )
  ### ---
  if (.bootstrap) {
    cat("Bootstrapping...\n")
    pb <- utils::txtProgressBar(min = 0, max = .B, style = 3)
    std.error <- map(.x = seq(.B), .f = function(i) {
      boot_data <- data[sample(x = seq(nrow(data)), size = nrow(data), replace = TRUE), ]
      boot_mu <- riskRegression::Score(
        object = list("Age + Gender" = fit0, "Age + Gender + Cluster" = fit1, "Age + Gender + Count" = fit2),
        formula = update(fit0$formula, . ~ 1),
        se.fit = 0,
        times = td,
        data = boot_data,
        contrasts = FALSE,
        progress.bar = NULL
      )
      boot_mu <- bind_rows(
        boot_mu$AUC$score %>% mutate(metric = "AUC") %>% rename(value = AUC),
        boot_mu$Brier$score %>% mutate(metric = "Brier") %>% rename(value = Brier)
      )
      boot_mu[["i"]] <- i
      utils::setTxtProgressBar(pb = pb, value = i)
      return(boot_mu)
    })
    close(pb)
    std.error_mu <- bind_rows(std.error) %>%
      group_by(model, times, metric) %>%
      summarise(
        boot_mu = mean(value),
        boot_se = sqrt(var(value))
      ) %>%
      ungroup()
    mu_out <- merge(mu_out, std.error_mu)
  }
  mu_out <- arrange(mu_out, metric, times, model) %>%
    mutate(Outcome = outcome, eGFR = egfr)
  # Return
  return(mu_out)
}

# Fit all models to each eGFR cohort
for (i in seq_along(ooo)) {
  n_acm <- glue("data/30-fitted-metrics-acm-eGFR={names(ooo)[i]}.RDS")
  if (!file.exists(n_acm)) {
    f_acm <- fit_metrics(data = ooo[[i]], outcome = "ACM", egfr = names(ooo)[i])
    saveRDS(object = f_acm, file = n_acm)
    rm(f_acm)
    gc()
  }
  n_mace3 <- glue("data/30-fitted-metrics-mace3-eGFR={names(ooo)[i]}.RDS")
  if (!file.exists(n_mace3)) {
    f_mace3 <- fit_metrics(data = ooo[[i]], outcome = "MACE3", egfr = names(ooo)[i])
    saveRDS(object = f_mace3, file = n_mace3)
    rm(f_acm)
    gc()
  }
  n_mace4 <- glue("data/30-fitted-metrics-mace4-eGFR={names(ooo)[i]}.RDS")
  if (!file.exists(n_mace4)) {
    f_mace4 <- fit_metrics(data = ooo[[i]], outcome = "MACE4", egfr = names(ooo)[i])
    saveRDS(object = f_mace4, file = n_mace4)
    rm(f_acm)
    gc()
  }
}

# Summarise metrics:
res <- map_dfr(.x = list.files(path = "data/", pattern = "30-fitted-metrics-", full.names = TRUE), .f = readRDS)
#
p_AUC <- res %>%
  filter(metric == "AUC") %>%
  ggplot(aes(x = times, y = value)) +
  geom_ribbon(aes(ymin = value - qnorm(1 - 0.05 / 2) * boot_se, ymax = value + qnorm(1 - 0.05 / 2) * boot_se, fill = model), alpha = 1 / 3) +
  geom_point(aes(color = model)) +
  geom_line(aes(color = model)) +
  facet_grid(eGFR ~ Outcome, labeller = label_both) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  coord_cartesian(ylim = c(0.5, 1)) +
  labs(x = "Time", y = "Time-Dependent AUC", color = "Model", fill = "Model") +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(legend.position = "bottom")
ggsave(plot = p_AUC, filename = "figures/30-AUC.pdf", height = 6 * sqrt(2), width = 6, device = cairo_pdf)
ggsave(plot = p_AUC, filename = "figures/30-AUC.png", height = 6 * sqrt(2), width = 6, device = ragg::agg_png, res = 600, units = "in")
#
p_Brier <- res %>%
  filter(metric == "Brier") %>%
  ggplot(aes(x = times, y = value)) +
  geom_ribbon(aes(ymin = value - qnorm(1 - 0.05 / 2) * boot_se, ymax = value + qnorm(1 - 0.05 / 2) * boot_se, fill = model), alpha = 1 / 3) +
  geom_point(aes(color = model)) +
  geom_line(aes(color = model)) +
  facet_grid(eGFR ~ Outcome, labeller = label_both) +
  scale_color_viridis_d() +
  scale_fill_viridis_d() +
  labs(x = "Time", y = "Brier Score", color = "Model", fill = "Model") +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(legend.position = "bottom")
ggsave(plot = p_Brier, filename = "figures/30-Brier.pdf", height = 6 * sqrt(2), width = 6, device = cairo_pdf)
ggsave(plot = p_Brier, filename = "figures/30-Brier.png", height = 6 * sqrt(2), width = 6, device = ragg::agg_png, res = 600, units = "in")
