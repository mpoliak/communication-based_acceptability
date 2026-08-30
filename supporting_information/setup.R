## =============================================================================
## Shared setup for the Supporting Information analyses
##
## Sourced by every section script. Loads the libraries and plotting theme,
## defines the helper functions the sections share, reads the data and builds
## `crit`, and prepares the by-item data frames and plotting helpers that
## Sections 1-3 and the diagnostics all draw on.
##
## Run in a UTF-8 locale (default in RStudio; for Rscript use e.g.
## LC_ALL=en_US.UTF-8). Fitted models are cached under models/; delete that
## folder's contents to refit from scratch.
## =============================================================================

library(tidyverse)
library(lme4)
library(lmerTest)
library(mgcv)
library(scam)
library(patchwork)

## Set plotting preferences
ggplot2::theme_set(ggplot2::theme_bw(base_size=18))
ggplot2::theme_update(panel.grid = ggplot2::element_blank(), 
                      strip.background = ggplot2::element_blank(),
                      legend.key = ggplot2::element_blank(),
                      panel.border = ggplot2::element_blank(),
                      axis.line = ggplot2::element_line(),
                      strip.text = ggplot2::element_text(face = "bold"),
                      plot.title = element_text(hjust = 0.5))

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

## Create the output folders if they are absent. models/ holds git-ignored model
## caches, so it does not exist in a fresh clone and save() below would fail.
dir.create("models", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)

## Helper functions
label_language <- function(language_code) {
  factor(dplyr::case_when(
    language_code == "da" ~ "Danish",
    language_code == "de" ~ "German",
    language_code == "en" ~ "English",
    language_code == "hi" ~ "Hindi",
    language_code == "fr" ~ "French",
    language_code == "ko" ~ "Korean",
    language_code == "ru" ~ "Russian",
    language_code == "zh" ~ "Mandarin",
    language_code == "en_exp" ~ "English\n(controlled)"
  ),
  ## Use the reverse of this order: levels=c("English\n(controlled)", "Russian", "Mandarin", "Korean", "Hindi", "German", "French", "English", "Danish")
  levels=c("Danish", "English", "French", "German", "Hindi", "Korean", 
           "Mandarin", "Russian", "English\n(controlled)")
  )
}


## ---- Marginal log-likelihood for scam/gam models with random effects ----
## scam/gam's logLik() is CONDITIONAL on the fitted random effects b; this returns
## the MARGINAL log-likelihood (integrating the b's out), i.e. the quantity that
## lmer() reports, so that model comparison does not credit a model for its fitted
## item- or language-specific random effects. Validated against lmer() for every
## random-effects structure used here (item/language intercepts and slopes) in
## ../marginal_loglik_test.R (outside this repo).
marginal_logLik <- function(m) {
  scale  <- m$sig2
  y      <- if (!is.null(m$y)) m$y else model.response(m$model)
  n      <- length(y)
  re_idx <- which(vapply(m$smooth, function(s) inherits(s, "random.effect"), logical(1)))
  if (!length(re_idx)) return(as.numeric(logLik(m)))       # no REs: nothing to marginalize
  re_lab <- vapply(m$smooth[re_idx], function(s) s$label, character(1))
  mu <- as.numeric(predict(m, exclude = re_lab))           # fixed part with random effects off
  Xp <- predict(m, type = "lpmatrix")                      # only the random-effect columns used below
  V  <- diag(scale, n)                                     # residual variance
  for (k in re_idx) {
    sm   <- m$smooth[[k]]
    Zk   <- Xp[, sm$first.para:sm$last.para, drop = FALSE]
    vark <- scale / (m$sp[[sm$label]] * diag(sm$S[[1]]))   # RE coefficient variance (accounts for
    V    <- V + Zk %*% (vark * t(Zk))                      #   mgcv's penalty rescaling of slopes)
  }
  r  <- y - mu
  ch <- chol(V)
  -0.5 * (n * log(2 * pi) + 2 * sum(log(diag(ch))) + sum(backsolve(ch, r, transpose = TRUE)^2))
}

## ---- By-language stats tibble with MARGINAL log-likelihoods ----
## Rebuilds a per-language stats tibble (for the likelihood-difference figures)
## from a fitted model list, using marginal_logLik() rather than the conditional
## logLik(). Called after each cache block so the figures use marginal values
## whether the models were just fit or loaded from disk.
marginal_stats_bylang <- function(model_list, model_map) {
  do.call(rbind, lapply(names(model_list), function(lg) {
    M <- model_list[[lg]]
    tibble(language = lg, model = names(M),
           log_likelihood = vapply(M, marginal_logLik, numeric(1)))
  })) %>%
    mutate(language = label_language(language),
           model = factor(model_map[model], levels = unname(model_map)))
}

## ---- Read and Prepare Data ----
## Read full data
dat <- read_csv("../data/full_data.csv")
surprisals <- read_csv("goldfish_surprisals/surprisal_goldfish.csv")

## Examine comprehension questions
comprehension <- dat %>% filter(grepl("comprehension", experiment))

dat <- dat %>% 
  left_join(
    surprisals %>% select(-sentence) %>% mutate(language = ifelse(language == "en_experimental", "en_exp", language)),
    by=c("language", "item")
  )

## Exclude participants with a low comprehension score
exclusions <- comprehension %>%
  group_by(participant_id) %>%
  mutate(
    correct_response = literal_response == response
    ) %>%
  summarize(
    accuracy = mean(correct_response)
  ) %>% 
  filter(
    accuracy < .8
  ) %>%
  pull(participant_id)

dat <- dat %>% 
  filter(!participant_id %in% exclusions)

## Get critical items
crit <- dat %>% 
  filter(experiment == "exp_b") %>%
  mutate(
    participant_id = as.factor(participant_id),
    item = as.factor(item),
    language = as.factor(language)
  )

## Remove nonsense conditions
crit <- crit %>% 
  filter(!grepl("non_sense", condition))

## Count number of participants per language
n_participants <- crit %>% 
  group_by(language) %>%
  summarize(
    n = n_distinct(participant_id)
  )

n_participants

table(crit$participant_id, crit$language)

## Verify the number of sentences per each language
n_sentences <- crit %>%
  group_by(language) %>%
  summarize(
    n = n_distinct(original_sentence)
  )

n_sentences

## Convert responses to numbers
crit$response_num <- as.numeric(substr(crit$response,1,1))

## Prepare the data frame
crit <- crit %>% 
  mutate(
    weighted_dl_distance = crit$dl_distance / crit$n_words,
    weighted_dl_distance_surprisal = crit$dl_distance / crit$sum_surprisal,
    grammatical = factor(
      ifelse(dl_distance == 0, "grammatical", "ungrammatical"),
      levels=c("ungrammatical", "grammatical")),
    item = as.factor(item),
    language = as.factor(language),
    response = as.factor(response_num)
  ) %>%
  select(
    participant_id, language, item, response, n_words, dl_distance, 
    weighted_dl_distance, grammatical, response_num, first_language,
    english_first_language, french_first_language, sum_surprisal, condition,
    weighted_dl_distance_surprisal
  )



## ---- By-item data frames ----
## Item means used by the shape-constrained GAMs in Sections 1-3 and by the
## diagnostics. `crit_by_items` is the dl_distance <= 3 subset the main
## analyses use; `crit_by_items_full` keeps the shuffled condition too.
crit_by_items_full <- crit %>% 
  mutate(
    weighted_dl_distance = crit$dl_distance / crit$n_words,
    grammatical = factor(
      ifelse(dl_distance == 0, "grammatical", "ungrammatical"),
      levels=c("ungrammatical", "grammatical")),
    item = as.factor(item),
    language = as.factor(language),
    response = as.factor(response),
    response_num = as.numeric(as.character(response)),
    participant_id = as.factor(participant_id),
    item = as.factor(item),
    language = as.factor(language)
  ) %>%
  group_by(language, item, n_words, dl_distance, weighted_dl_distance, 
           grammatical, sum_surprisal, condition, weighted_dl_distance_surprisal) %>%
  summarize(
    response_num = mean(response_num)
  ) %>%
  mutate(
    grammatical_numeric = ifelse(grammatical == "grammatical", .5, -.5),
    item = as.factor(gsub("[sml]$", "", item))
  )

crit_by_items <- crit_by_items_full %>%
  filter(
    dl_distance <= 3
  )


## ---- Shared plotting helpers ----
## Colour scale, DL-distance buckets, the model-prediction grid builder and the
## by-item scatter data, all shared by the Section 1-3 figures.
#### ---- Helper Functions 
YLIM <- c(0, 8)

DL_COLORS <- c("0"="#EB7AA3", "1-2"="#1E88E5", "3"="#FFC107", "4+"="#BBBBBB")
DL_BUCKET_LEVELS <- c("0", "1-2", "3", "4+")

bucket_dl <- function(dl_distance) {
  dplyr::case_when(
    dl_distance == 0 ~ "0",
    dl_distance %in% 1:2 ~ "1-2",
    dl_distance == 3 ~ "3",
    dl_distance >= 4 ~ "4+"
  ) %>%
    factor(levels = DL_BUCKET_LEVELS)
}

scale_dl_color <- function() {
  ggplot2::scale_color_manual(
    name   = "Damerau-Levenshtein Distance",
    values = DL_COLORS,
    breaks = DL_BUCKET_LEVELS,
    limits = DL_BUCKET_LEVELS,
    drop   = FALSE,
    guide  = ggplot2::guide_legend(
      override.aes = list(alpha = 1, linewidth = 0, size = 3)
    )
  )
}


## `fits` is the by-language model list to take the item factor from. The item
## level is only a placeholder so predict() sees every term in the model -- the
## item random effect is excluded from the predictions themselves -- but it has
## to be a level the model knows, so it is read from the fits being predicted.
get_prediction_grid <- function(
    fits,
    lang,
    model_type,
    n_words=c(4, 40), 
    dl_distance=c(0, 1, 3)
) {
  expand.grid(
    n_words = seq(n_words[1] - .25, n_words[2] + .25, 1),
    dl_distance = dl_distance
  ) %>%
    mutate(
      weighted_dl_distance = dl_distance / n_words,
      grammatical_numeric = ifelse(dl_distance == 0, .5, -.5),
      item=factor(fits[[lang]][[model_type]]$model$item[1], levels=levels(fits[[lang]][[model_type]]$model$item))
    )
}


### ---- Baseline data for plotting (item means) ----
plot_data_by_item <- crit %>% 
  mutate(
    weighted_dl_distance = crit$dl_distance / crit$n_words,
    grammatical = factor(
      ifelse(dl_distance == 0, "grammatical", "ungrammatical"),
      levels=c("ungrammatical", "grammatical")),
    item = as.factor(item),
    language = as.factor(language),
    response = as.factor(response),
    response_num = as.numeric(as.character(response)),
    participant_id = as.factor(participant_id),
    item = as.factor(item),
    language = as.factor(language)
  ) %>%
  group_by(language, item, n_words, dl_distance, weighted_dl_distance, 
           grammatical, condition) %>%
  summarize(
    response_num = mean(response_num)
  ) %>%
  mutate(
    grammatical_numeric = ifelse(grammatical == "grammatical", .5, -.5)
  ) %>%
  dplyr::group_by(language, item, n_words, dl_distance) %>%
  dplyr::summarise(response = mean(response_num), .groups = "drop") %>%
  dplyr::mutate(
    dl_distance_bucket = bucket_dl(dl_distance),
    language           = label_language(language)
  )
