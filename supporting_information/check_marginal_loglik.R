## =============================================================================
## Check: marginal_logLik() against lmer()
##
## scam/gam's logLik() is CONDITIONAL on the fitted random effects b.
## marginal_logLik() (defined in setup.R; ../analysis/main_analysis.R has an
## identical copy) returns the MARGINAL log-likelihood, with the b's integrated
## out: the multivariate normal density of y with mean mu (the fixed part) and
## covariance V = Z G Z' + sigma^2 I. That is the quantity lmer() reports.
##
## Each test fits the same model as a gam (method = "ML") and as an lmer
## (REML = FALSE) and checks that marginal_logLik(gam) equals logLik(lmer).
## mgcv's s(g, x, bs = "re") is an UNCORRELATED random slope, so the lmer twin
## uses (1 | g) + (0 + x | g). The tests cover the random-effect terms the
## analyses use (item intercepts and slopes, language intercepts and slopes),
## alone and combined:
##   T1-T2  German by-item means (as in the per-experiment fits)
##   T2b    simulated data with clearly non-zero item intercept and slope
##          variances (the real item-slope variance can be estimated at zero)
##   T3-T5  by-item means pooled over all nine experiments (as in the pooled fits)
##
## Small fits only. Not part of run_all.R; run it on its own:
##   LC_ALL=en_US.UTF-8 Rscript supporting_information/check_marginal_loglik.R
## It stops with an error if any test fails.
## =============================================================================

library(tidyverse)
library(lme4)
library(mgcv)

## Set plotting preferences
ggplot2::theme_set(ggplot2::theme_bw(base_size=18))
ggplot2::theme_update(panel.grid = ggplot2::element_blank(),
                      strip.background = ggplot2::element_blank(),
                      legend.key = ggplot2::element_blank(),
                      panel.border = ggplot2::element_blank(),
                      axis.line = ggplot2::element_line(),
                      strip.text = ggplot2::element_text(face = "bold"),
                      plot.title = element_text(hjust = 0.5))
options(ggplot2.discrete.colour = c("#A31F34", "#8A8B8C"))
options(ggplot2.discrete.fill = c("#A31F34", "#8A8B8C"))

## Set working directory
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg)) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  }
  getwd()
}
setwd(get_script_dir())

## marginal_logLik() and the by-item data frame `crit_by_items` come from
## setup.R (reading the data only; it fits no models)
source("setup.R")

## ---- Test helper ----
## Absolute tolerance in log-likelihood units
tolerance <- 0.05

check_marginal <- function(name, gam_formula, lmer_formula, data) {
  mod_gam  <- mgcv::gam(gam_formula, data = data, method = "ML")
  mod_lmer <- lme4::lmer(lmer_formula, data = data, REML = FALSE,
                         control = lme4::lmerControl(check.conv.singular = "ignore"))
  tibble(
    test = name,
    marginal_gam = marginal_logLik(mod_gam),
    lmer = as.numeric(logLik(mod_lmer))
  ) %>%
    mutate(diff = marginal_gam - lmer,
           pass = abs(diff) < tolerance)
}

## ---- Test data ----
d_item <- crit_by_items %>%
  filter(language == "de") %>%
  droplevels()

d_comb <- crit_by_items %>%
  droplevels()

set.seed(1)
n_items <- 60
n_reps <- 6
sim <- tibble(
  item = factor(rep(1:n_items, each = n_reps)),
  x = rnorm(n_items * n_reps)
)
sim$y <- 2 + 0.5 * sim$x + rnorm(n_items, 0, 1.2)[sim$item] +
  rnorm(n_items, 0, 0.8)[sim$item] * sim$x + rnorm(n_items * n_reps, 0, 0.7)

## ---- Tests ----
results <- bind_rows(
  check_marginal(
    "T1  item intercept",
    response_num ~ sum_surprisal + s(item, bs = "re"),
    response_num ~ sum_surprisal + (1 | item),
    d_item),
  check_marginal(
    "T2  item intercept + item slope",
    response_num ~ sum_surprisal + s(item, bs = "re") + s(item, sum_surprisal, bs = "re"),
    response_num ~ sum_surprisal + (1 | item) + (0 + sum_surprisal | item),
    d_item),
  check_marginal(
    "T2b item intercept + item slope (simulated)",
    y ~ x + s(item, bs = "re") + s(item, x, bs = "re"),
    y ~ x + (1 | item) + (0 + x | item),
    sim),
  check_marginal(
    "T3  crossed item + language intercepts",
    response_num ~ n_words + s(item, bs = "re") + s(language, bs = "re"),
    response_num ~ n_words + (1 | item) + (1 | language),
    d_comb),
  check_marginal(
    "T4  item + language intercepts, language slope",
    response_num ~ n_words + s(item, bs = "re") + s(language, bs = "re") +
      s(language, n_words, bs = "re"),
    response_num ~ n_words + (1 | item) + (1 | language) + (0 + n_words | language),
    d_comb),
  check_marginal(
    "T5  item + language intercepts, two language slopes",
    response_num ~ weighted_dl_distance + n_words + s(item, bs = "re") +
      s(language, bs = "re") + s(language, weighted_dl_distance, bs = "re") +
      s(language, n_words, bs = "re"),
    response_num ~ weighted_dl_distance + n_words + (1 | item) + (1 | language) +
      (0 + weighted_dl_distance | language) + (0 + n_words | language),
    d_comb)
)

cat("\n==== marginal_logLik() vs lmer() ====\n")
print(as.data.frame(results), digits = 8, row.names = FALSE)

if (!all(results$pass)) {
  stop("marginal_logLik() differs from lmer() by more than ", tolerance,
       " in: ", paste(results$test[!results$pass], collapse = "; "))
}
cat("\nAll", nrow(results), "tests pass (|difference| <", tolerance, ").\n")
