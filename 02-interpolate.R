# Script to interpolate kidney function at given dates
library(tidyverse)
library(lme4)
library(SCREAM) # install from https://github.com/ellessenne/SCREAM
library(hrbrthemes)
library(cowplot)

# 'data' here is the input dataset
# This is the same dataset used before, but potentially with
#  extra information and covariates.
data <- readRDS(file = "path/to/file")
str(data)
head(data)

# Import fitted LMM
lmm <- readRDS("data/01-fitted-lmm.RDS")

# Interpolate
.egfr <- seq(15, 90, by = 15) # <- will interpolate at these values
interpolated_data <- SCREAM::interpolate_date(
  fit = lmm,
  egfr = .egfr,
  id = "lopnr",
  origin = as.Date("2006-01-01"),
  time_scale = "years"
)
interpolated_data$lopnr <- as.numeric(rownames(coef(lmm)$lopnr))
interpolated_data <- left_join(
  interpolated_data,
  distinct(data, lopnr, first, last),
  by = "lopnr"
)

# Interpolated dates must be between 'first' and 'last'
# These ('first', 'last') are columns in data identifying the
# first and last observation date for each subject.
for (level in .egfr) {
  interpolated_data[[paste0("time_", level)]] <- ifelse(
    interpolated_data[[paste0("time_", level)]] >= interpolated_data$first &
      interpolated_data[[paste0("time_", level)]] <= interpolated_data$last,
    interpolated_data[[paste0("time_", level)]],
    NA
  )
  interpolated_data[[paste0("time_", level)]] <- structure(interpolated_data[[paste0("time_", level)]], class = "Date")
}

# Make splits
splits <- map(.x = .egfr, .f = function(value) {
  filter1 <- interpolated_data[!is.na(interpolated_data[[paste0("time_", value)]]), ]
  data.frame(
    lopnr = filter1$lopnr,
    index_date = filter1[[paste0("time_", value)]]
  )
})
names(splits) <- .egfr

# Export
saveRDS(object = splits, file = "data/02-interpolated-index-dates.RDS")

# Compare included vs excluded subjects at eGFR level
rm(lmm, interpolated_data)
gc()
data <- distinct(data, lopnr, dob, female, N)
for (i in seq_along(splits)) {
  data[[paste0("included_eGFR", names(splits)[i])]] <- ifelse(data$lopnr %in% splits[[i]]$lopnr, 1, 0)
}
data <- pivot_longer(data, cols = starts_with("included")) %>%
  separate(col = "name", into = c("tmp", "eGFR"), sep = "_") %>%
  mutate(eGFR = substr(eGFR, 5, 6)) %>%
  mutate(value = factor(value, levels = 0:1, labels = c("Not Included", "Included"))) %>%
  mutate(eGFR = factor(eGFR, levels = c("90", "75", "60", "45", "30", "15"))) %>%
  group_by(eGFR, value) %>%
  summarise(
    dob_median = median(dob),
    dob_lower = quantile(dob, 0.25, type = 1),
    dob_upper = quantile(dob, 0.75, type = 1),
    female_p = mean(female),
    female_n = sum(female),
    N_median = median(N),
    N_lower = fivenum(N)[2],
    N_upper = fivenum(N)[4]
  )
#
pA <- ggplot(data, aes(x = eGFR, y = dob_median, color = value)) +
  geom_errorbar(aes(ymin = dob_lower, ymax = dob_upper), position = position_dodge(width = 0.5), width = 0.5) +
  geom_point(position = position_dodge(width = 0.5)) +
  scale_y_date(breaks = c(as.Date("1930-01-01"), as.Date("1940-01-01"), as.Date("1950-01-01"), as.Date("1960-01-01"), as.Date("1970-01-01"), as.Date("1980-01-01"), as.Date("1990-01-01")), date_labels = "%Y") +
  scale_color_ipsum() +
  scale_fill_ipsum() +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(legend.position = "none") +
  labs(y = "Median Date of Birth (with I.Q.I.)")
pB <- ggplot(data, aes(x = eGFR, y = N_median, color = value)) +
  geom_errorbar(aes(ymin = N_lower, ymax = N_upper), position = position_dodge(width = 0.5), width = 0.5) +
  geom_point(position = position_dodge(width = 0.5)) +
  scale_color_ipsum() +
  scale_fill_ipsum() +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(legend.position = "bottom") +
  labs(y = "Median N. of eGFR Measurements (with I.Q.I.)", color = "", fill = "")
pC <- ggplot(data, aes(x = eGFR, y = female_p, color = value, fill = value)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6) +
  scale_y_percent(breaks = seq(0, 0.6, by = 0.1)) +
  coord_cartesian(ylim = c(0, 0.6)) +
  scale_color_ipsum() +
  scale_fill_ipsum() +
  theme_minimal(base_size = 12, base_family = "Arial Narrow") +
  theme(legend.position = "none") +
  labs(y = "Proportion of Sex = Female")
pABC <- plot_grid(pA, pB, pC, nrow = 1, labels = LETTERS, align = "hv", axis = "tblr")
ggsave(plot = pABC, filename = "figures/02-included-vs-excluded.pdf", width = 8.5, height = 4.25, device = cairo_pdf)
ggsave(plot = pABC, filename = "figures/02-included-vs-excluded.png", width = 8.5, height = 4.25, device = ragg::agg_png, res = 600, units = "in")
