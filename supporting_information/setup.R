## =============================================================================
## Shared setup for the Supporting Information analyses
##
## Sourced by every section script. Loads the libraries and plotting theme,
## defines the helper functions the sections share, reads the data and builds
## `crit`, and prepares the by-item data frames (used by Sections 1-3 and
## check_marginal_loglik.R) and plotting helpers (Sections 1-5).
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

## Print warnings as they happen, so they show up next to the step that raised
## them in the run logs
options(warn = 1)

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
  ## Only while the path still resolves: once a script has moved into its
  ## folder, a relative --file path no longer does, and getwd() is the folder
  if (length(file_arg) && file.exists(sub("^--file=", "", file_arg))) {
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
## caches, so it does not exist in a fresh clone and the save() calls of
## Sections 1-3 would fail.
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
  ## Figures that list the languages top to bottom reverse this order (limits = rev)
  levels=c("Danish", "English", "French", "German", "Hindi", "Korean", 
           "Mandarin", "Russian", "English\n(controlled)")
  )
}


## ---- Marginal log-likelihood for scam/gam models with random effects ----
## scam/gam's logLik() is CONDITIONAL on the fitted random effects b; this returns
## the MARGINAL log-likelihood (integrating the b's out), i.e. the quantity that
## lmer() reports, so that model comparison does not credit a model for its fitted
## item- or language-specific random effects. check_marginal_loglik.R checks it
## against lmer() for each random-effect term used here (item/language
## intercepts and slopes).
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

## ---- Random-effect smooth labels ----
## Labels of every random-effect smooth in a scam/gam model (item and language
## intercepts and slopes). Passed to predict(exclude = ...) so that predictions
## are population-level: no item effect (including the en_exp by-item slopes)
## and, in the pooled fits, no language effect.
random_effect_labels <- function(m) {
  re_idx <- which(vapply(m$smooth, function(s) inherits(s, "random.effect"), logical(1)))
  vapply(m$smooth[re_idx], function(s) s$label, character(1))
}

## ---- MORCELA weights and score ----
## MORCELA (Tjuatja, Neubig, Linzen & Hao 2025) is SLOR with the unigram and
## length corrections allowed to be weaker or stronger than SLOR assumes:
##
##     MORCELA = ( log p(sentence) - beta * log p_unigram(sentence) + gamma ) / length
##
## The weights are learned from the ratings, as in that paper: regress the ratings
## on the three terms and read beta and gamma off the coefficients (as in
## ../analysis/main_analysis.R). Split into two functions so that the weights
## themselves can be reported.
morcela_weights <- function(model_data) {
  morcela_regression <- lm(
    response_num ~ I(log_p_sentence / n_words) + I(log_p_unigram / n_words) + I(1 / n_words),
    data = model_data
  )
  ## Coefficients in formula order: intercept, sentence, unigram, 1/length
  weights <- unname(coef(morcela_regression))
  list(
    beta  = -weights[3] / weights[2],
    gamma =  weights[4] / weights[2]
  )
}

morcela_score <- function(model_data, weights) {
  (model_data$log_p_sentence - weights$beta * model_data$log_p_unigram + weights$gamma) /
    model_data$n_words
}

## ---- Read and Prepare Data ----
## Read full data
dat <- read_csv("../data/full_data.csv", show_col_types = FALSE)

## Sentence log-probabilities from the Goldfish models (log2, i.e. bits), one row
## per sentence per unigram model, written by data/get_slor_morcela_scores.py. We
## use the model-initial unigram throughout. The item ids there are bare, while
## full_data.csv prefixes them with the folder name, so rebuild that prefix (the
## controlled English experiment is en_exp in the data but lives in the
## en_experimental folder).
scores <- read_csv(
  "../data/slor_morcela_scores.csv",
  col_types = cols(item = col_character())
) %>%
  filter(unigram_model == "model_initial") %>%
  mutate(
    item = paste(ifelse(language == "en_exp", "en_experimental", language), item, sep = "_")
  )

## Comprehension-question rows (scored below for the exclusions)
comprehension <- dat %>% filter(grepl("comprehension", experiment))

## Information(s) is the surprisal of the intended (original) sentence in bits.
## It is a property of the item, so it is the same in every condition of an item.
## SLOR and MORCELA instead score each condition's own sentence, so they take
## that condition's log-probabilities.
dat <- dat %>%
  left_join(
    scores %>%
      filter(condition == "original") %>%
      transmute(language, item, sum_surprisal = -log_p_sentence),
    by=c("language", "item")
  ) %>%
  left_join(
    scores %>% select(language, item, condition, log_p_sentence, log_p_unigram),
    by=c("language", "item", "condition")
  ) %>%
  mutate(
    weighted_dl_distance_surprisal = dl_distance / sum_surprisal,
    ## SLOR (Lau, Clark & Lappin 2017), with length in tokens
    slor = (log_p_sentence - log_p_unigram) / n_words
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
## Each response label starts with its scale point (1-7)
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
    weighted_dl_distance_surprisal, log_p_sentence, log_p_unigram, slor
  )

## Every critical sentence must have a score
stopifnot(
  !anyNA(crit$sum_surprisal),
  !anyNA(crit$log_p_sentence),
  !anyNA(crit$log_p_unigram)
)


## ---- By-item data frames ----
## Item x condition means (averaged over participants) used by the
## shape-constrained GAMs in Sections 1-3 and by check_marginal_loglik.R.
## `crit_by_items` is the dl_distance <= 3 subset the main analyses use;
## `crit_by_items_full` (Section 2) also keeps the cells at DL >= 4. The filter
## is on DL distance, not condition: the 20 shuffled cells that came out at DL
## 2-3 stay in, the one three-exchange cell at DL 4 drops out, and the 9
## one-exchange cells that the exchange left unchanged (DL 0) count as
## error-free.
crit_by_items_full <- crit %>%
  group_by(language, item, n_words, dl_distance, weighted_dl_distance,
           grammatical, sum_surprisal, condition, weighted_dl_distance_surprisal,
           log_p_sentence, log_p_unigram, slor) %>%
  summarize(
    response_num = mean(response_num),
    .groups = "drop"
  ) %>%
  mutate(
    grammatical_numeric = ifelse(grammatical == "grammatical", .5, -.5),
    ## Strip the s/m/l length suffix of the controlled English item ids after
    ## the means are taken, so the three versions of a base item share one item
    ## level (72 instead of 216)
    item = as.factor(gsub("[sml]$", "", item))
  )

crit_by_items <- crit_by_items_full %>%
  filter(
    dl_distance <= 3
  )


## ---- Shared plotting helpers ----
## Colour scale and DL-distance buckets (Sections 1-3 and 5); the prediction-grid
## builders, the print theme and the by-item scatter data (Sections 1-3).
YLIM <- c(0, 8)

DL_COLORS <- c("0"="#EB7AA3", "1-2"="#1E88E5", "3"="#FFC107", "4+"="#BBBBBB")
DL_BUCKET_LEVELS <- c("0", "1-2", "3", "4+")

## The binary model predicts the same rating at every DL distance above 0, so
## its panels draw a single "any error" line, in a neutral dark grey that is no
## model's colour
DL_COLORS_BINARY <- c(DL_COLORS["0"], "1-2" = "grey25", "3" = "grey25", "4+" = "grey25")

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
      override.aes = list(alpha = 1, linewidth = 0, size = 1.6)
    )
  )
}

## The curve figures of Sections 1-3 (Figures S1, S3 and S5) are drawn at their
## print size, like the main text's Figure 1 (../analysis/main_analysis.R): 7 in
## wide, with 8-pt text. This is the theme above at 8 pt, with axis lines
## thinned to match.
theme_print_8pt <- theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.key = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5))


## `fits` is the by-language model list to take the item factor from. The item
## level is only a placeholder so predict() sees every term in the model -- the
## item random effect is excluded from the predictions themselves -- but it has
## to be a level the model knows, so it is read from the fits being predicted.
## `n_words` is the range to cover (the experiment's shortest and longest
## sentence). The grid starts 0.25 below the shortest and runs in steps of 1
## token up to at most 0.25 above the longest, so the curves span the data.
get_prediction_grid <- function(
    fits,
    lang,
    model_type,
    n_words,
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

## The same grid over sentence surprisal (bits) instead of length, for the
## models whose Information(s) is surprisal (steps of 1 bit over the range given
## in `sum_surprisal`, padded in the same way).
get_prediction_grid_surprisal <- function(
    fits,
    lang,
    model_type,
    sum_surprisal,
    dl_distance=c(0, 1, 3)
) {
  expand.grid(
    sum_surprisal = seq(sum_surprisal[1] - .25, sum_surprisal[2] + .25, 1),
    dl_distance = dl_distance
  ) %>%
    mutate(
      weighted_dl_distance_surprisal = dl_distance / sum_surprisal,
      grammatical_numeric = ifelse(dl_distance == 0, .5, -.5),
      item=factor(fits[[lang]][[model_type]]$model$item[1], levels=levels(fits[[lang]][[model_type]]$model$item))
    )
}

## Colours, shapes and legend labels of the five models in the model-comparison
## figures (as in ../analysis/main_analysis.R). SLOR and MORCELA get
## filled shapes (23 and 25 take both a colour and a fill), so all five markers
## are solid.
MODEL_LEVELS <- c("Binary", "Gradient", "Communication", "SLOR", "MORCELA")
MODEL_LABELS <- c("Binary\nGrammaticality", "Gradient\nGrammaticality",
                  "Communication", "SLOR", "MORCELA")
MODEL_COLORS <- c(
  "Binary"        = "#939DB0",
  "Gradient"      = "#5673C5",
  "Communication" = "#003736",
  "SLOR"          = "#F7B957",
  "MORCELA"       = "#A90318"
)
MODEL_SHAPES <- c(
  "Binary"        = 16,
  "Gradient"      = 17,
  "Communication" = 15,
  "SLOR"          = 23,
  "MORCELA"       = 25
)


### ---- Baseline data for plotting (item means) ----
## Item x condition means, then averaged over the conditions of an item that
## share a DL distance (one point per item and DL distance), and deliberately
## WITHOUT the filter(dl_distance <= 3): the scatter of the curve figures also
## shows the sentences at DL >= 4 (nearly all of the shuffled condition), in the
## grey "4+" bucket. Sentence surprisal and length are properties of the item,
## so the same points serve the figures over surprisal and over length.
plot_data_by_item <- crit %>%
  group_by(language, item, n_words, sum_surprisal, dl_distance, condition) %>%
  summarize(response_num = mean(response_num), .groups = "drop") %>%
  group_by(language, item, n_words, sum_surprisal, dl_distance) %>%
  summarize(response = mean(response_num), .groups = "drop") %>%
  mutate(
    dl_distance_bucket = bucket_dl(dl_distance),
    language           = label_language(language)
  )
