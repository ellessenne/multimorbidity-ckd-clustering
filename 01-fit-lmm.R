### This script fit a linear mixed model for eGFR over time
library(tidyverse)
library(lme4)
library(formattable)
library(ggtext)

# 'data' here is the input dataset
# Relevant columns are:
# - 'egfr', observed eGFR values
# - 'time', at which eGFR values are recorded
# - 'lopnr', unique ID number identifying each study subject
data <- readRDS(file = "path/to/file") %>%
  select(egfr, time, lopnr)
str(data)
head(data)

# Fit a linear mixed model on the input dataset
# By default, it assumes an unstructured correlation structure,
# that is, the two random effects can be correlated (with covariance
# estimated from data).
lmm <- lmer(egfr ~ time + (time | lopnr), data = data, verbose = TRUE)
# Correlation matrix of the random effects:
vcov(lmm, correlation = TRUE)@factors$correlation

# Check optimiser fit
lmm_af <- allFit(lmm)
lmm_af # <- ok if difference in -logl is small!
rm(lmm_af)
gc()

# Final check, fit model after re-scaling time
# Should converge to same likelihood value if not a problem
data$time2 <- scale(data$time)
lmm2 <- lmer(egfr ~ time2 + (time2 | lopnr), data = data, verbose = TRUE)
# Converges to - essentially - the same REML value
# Additionally, fits the same variance of random intercept and residual variance
rm(lmm2)
gc

# Export LMM for later use
saveRDS(object = lmm, file = "data/01-fitted-lmm.RDS")
