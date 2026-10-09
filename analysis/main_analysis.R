library(tidyverse)
library(mgcv)
library(scam)
library(lme4)
library(lmerTest)
library(patchwork)

## Print warnings as they happen, so that they land next to the step that raised
## them in the log
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
## outputs/ holds the tables of statistics reported in the paper.
dir.create("models", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)

## ---- Marginal log-likelihood for scam/gam models with random effects ----
## scam/gam's logLik() is CONDITIONAL on the fitted random effects b; this returns
## the MARGINAL log-likelihood (integrating the b's out), i.e. the quantity that
## lmer() reports, so that model comparison does not credit a model for its fitted
## item- or language-specific random effects.
## ../supporting_information/check_marginal_loglik.R checks it against lmer() for
## each random-effect term used here (item/language intercepts and slopes).
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

## The labels of every random-effect smooth in a scam/gam model: the item
## intercepts, the by-item slopes of the controlled English experiment, and, in
## the pooled fits, the language intercepts and slopes. Passed to predict()'s
## exclude argument to get what a model predicts from its predictors alone.
random_effect_labels <- function(m) {
  vapply(
    Filter(function(s) inherits(s, "random.effect"), m$smooth),
    function(s) s$label,
    character(1)
  )
}

## The labels of all smooths in a model, used to check that a cached model has
## the specification this script expects
smooth_labels <- function(m) {
  vapply(m$smooth, function(s) s$label, character(1))
}

## ---- Read and Prepare Data ----
## Read full data
dat <- read_csv("../data/full_data.csv")

## Examine comprehension questions
comprehension <- dat %>% filter(grepl("comprehension", experiment))

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
  filter(experiment == "exp_b")

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

## Count the number of sentences per each language
n_sentences <- crit %>%
  group_by(language) %>%
  summarize(
    n = n_distinct(original_sentence)
  )

n_sentences

## Convert responses to numbers
crit$response <- as.numeric(substr(crit$response,1,1))

english_by_bucket <- crit %>% filter(language == "en") %>%
  select(
    item, n_words
  ) %>%
  distinct() %>%
  mutate(
    item = as.numeric(gsub("en_", "", item)),
    bucket = ifelse(
      item <= 24,
      "short",
      ifelse(
        item <= 48,
        "medium",
        "long"
      )
    )
  )

## Get quantiles .0, .25, .5, .75, 1 by bucket
lapply(unique(english_by_bucket$bucket), function(current_bucket) {
  print(paste0("Processing bucket: ", current_bucket))
  quantile(
    english_by_bucket %>% filter(bucket == current_bucket) %>% pull(n_words),
    probs = c(0, .25, .5, .75, 1)
  )
})

## ---- Add Surprisal and the Language-Model Scores ----
## One source of surprisal for every model: data/slor_morcela_scores.csv, written
## by data/get_slor_morcela_scores.py. It holds log2 probabilities (bits, signed
## negative) from the Goldfish 1000mb model of each language, scored on the
## lowercased, unpunctuated perceived_sentence.
##
## The item ids there are bare ("1", or "1s"/"1m"/"1l" in the controlled English
## experiment), while full_data.csv prefixes them with the folder name, so we
## rebuild the same prefix before joining. Note that the controlled English
## experiment is called en_exp in the data but lives in the en_experimental folder.
##
## We use the model-initial unigram throughout. wordfreq is kept in the file as a
## robustness check: it gives the same answer everywhere, but it does not know 9%
## of the Danish word tokens or 10% of the Korean ones, so it cannot be the
## default in a cross-linguistic analysis.
scores <- read_csv(
  "../data/slor_morcela_scores.csv",
  col_types = cols(item = col_character())
) %>%
  filter(unigram_model == "model_initial") %>%
  mutate(
    item = paste(
      ifelse(language == "en_exp", "en_experimental", language),
      item,
      sep = "_"
    )
  )

## Information(s) in the communication account is the surprisal of the intended
## sentence, i.e. of the original. It is a property of the item, so it is the
## same in all of an item's conditions, just as length was.
original_surprisal <- scores %>%
  filter(condition == "original") %>%
  transmute(
    language,
    item,
    sum_surprisal = -log_p_sentence
  )

stopifnot(!anyDuplicated(original_surprisal[c("language", "item")]))

## SLOR and MORCELA instead score the sentence that was actually read, so they
## take the log probabilities of each condition.
n_crit_rows <- nrow(crit)

crit <- crit %>%
  left_join(
    original_surprisal,
    by = c("language", "item")
  ) %>%
  left_join(
    scores %>%
      select(language, item, condition, log_p_sentence, log_p_unigram),
    by = c("language", "item", "condition")
  )

## Every critical rating must have found exactly one score
stopifnot(
  nrow(crit) == n_crit_rows,
  !anyNA(crit$sum_surprisal),
  !anyNA(crit$log_p_sentence),
  !anyNA(crit$log_p_unigram)
)

### ---- Correlation between surprisal and length ----
## Over the original sentences, one per item (so 216 in the controlled English
## experiment, where each base item comes in 3 lengths). Length and surprisal are
## both properties of the item's original sentence, so we take them from all of
## an item's conditions: after the exclusions, a few controlled English items
## have no rating left in the original condition itself.
original_sentences <- crit %>%
  distinct(language, item, n_words, sum_surprisal)

stopifnot(!anyDuplicated(original_sentences[c("language", "item")]))

surprisal_length_correlation <- function(sentences, experiment_code) {
  correlation <- cor.test(sentences$n_words, sentences$sum_surprisal)
  tibble(
    experiment = experiment_code,
    n_sentences = nrow(sentences),
    r = unname(correlation$estimate),
    ci_low = correlation$conf.int[1],
    ci_high = correlation$conf.int[2],
    p = correlation$p.value
  )
}

surprisal_length_correlations <- bind_rows(
  lapply(sort(unique(original_sentences$language)), function(current_language) {
    surprisal_length_correlation(
      original_sentences %>% filter(language == current_language),
      current_language
    )
  }),
  surprisal_length_correlation(original_sentences, "all_experiments")
)

surprisal_length_correlations

write_csv(surprisal_length_correlations, "outputs/surprisal_length_correlations.csv")

## ---- Main Analyses ----
### ---- Linear Model ----
## en_exp is left out of this loop because it is fit below with its full by-item
## random-effects structure.
lmer_fits <- lapply(setdiff(unique(crit$language), "en_exp"), function(current_language) {
  print(paste0("Processing language: ", current_language))
  lmer(
  response ~ dl_distance * sum_surprisal +
    (dl_distance * sum_surprisal|participant_id) +
    (dl_distance|item),
  data=crit %>% filter(
    language == current_language,
    dl_distance <= 3
    ) %>%
    mutate(
      sum_surprisal = scale(sum_surprisal)[,1],
      dl_distance = scale(dl_distance)[,1]
    ),
  control = lmerControl(
    optimizer ="bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
  )
})

names(lmer_fits) <- setdiff(unique(crit$language), "en_exp")

## Fit full RE for en_exp.
## In this experiment each of the 72 base items was written in 3 lengths (short,
## medium, long), and each of those crossed with the 4 grammaticality conditions.
## The item id in the data encodes the length version (e.g. en_experimental_62m),
## so we strip the trailing s/m/l to recover the 72 base items. Surprisal then
## varies WITHIN item, which is what makes the by-item surprisal slope identifiable.
lmer_fits[["en_exp"]] <- lmer(
  response ~ dl_distance * sum_surprisal +
    (dl_distance * sum_surprisal|participant_id) +
    (dl_distance * sum_surprisal|item),
  data=crit %>% filter(
    language == "en_exp",
    dl_distance <= 3
  ) %>%
    mutate(
      item = as.factor(gsub("[sml]$", "", item)),
      sum_surprisal = scale(sum_surprisal)[,1],
      dl_distance = scale(dl_distance)[,1]
    ),
  ## nlminb rather than bobyqa: with this random-effects structure bobyqa stops
  ## ~27 log-likelihood units short of the optimum and fails the convergence
  ## check. allFit() shows nlminb and NLopt reaching the same optimum; nlminb is
  ## the one that passes the check. (In the other 8 experiments bobyqa already
  ## reaches the optimum.) The fit is still singular: the maximal structure is
  ## more than the data support.
  control = lmerControl(
    optimizer ="nlminbwrap",
    optCtrl = list(eval.max = 1e5, iter.max = 1e5)
  )
)

lapply(lmer_fits, summary)

## Convergence messages (NULL when there are none) and singular fits
lapply(lmer_fits, function(model_fit) {
  list(
    messages = model_fit@optinfo$conv$lme4$messages,
    singular = isSingular(model_fit)
  )
})

## Save the fixed effects of every experiment
lmm_dl_by_surprisal <- bind_rows(lapply(names(lmer_fits), function(current_language) {
  coefficients <- summary(lmer_fits[[current_language]])$coefficients
  tibble(
    experiment = current_language,
    term = rownames(coefficients),
    estimate = coefficients[, "Estimate"],
    se = coefficients[, "Std. Error"],
    t = coefficients[, "t value"],
    p = coefficients[, "Pr(>|t|)"]
  )
}))

write_csv(lmm_dl_by_surprisal, "outputs/lmm_dl_by_surprisal.csv")

## Get largest p-value
## DL Distance
lapply(lmer_fits, function(model_fit) {
  summary(model_fit)$coefficients[2,5]
}) %>% unlist %>% max

## Surprisal
lapply(lmer_fits, function(model_fit) {
  summary(model_fit)$coefficients[3,5]
}) %>% unlist %>% max

## Interaction
lapply(lmer_fits, function(model_fit) {
  summary(model_fit)$coefficients[4,5]
}) %>% unlist %>% max

### ---- Effects of Surprisal in Well-formed Materials ----
well_formed <- crit %>% filter(dl_distance == 0)

## sum_surprisal is a between-item predictor, so the by-item random intercept is
## needed for its standard error to be correct. For these 8 languages surprisal is
## constant within an item, so only an intercept (not a slope) can be estimated by
## item. Surprisal is z-scored within each experiment, so that the slopes are on
## the same scale and can be compared with the t-test below.
surprisal_effects <- lapply(setdiff(unique(well_formed$language), "en_exp"), function(current_language) {
  print(paste0("Processing language: ", current_language))
  lmer(
    response ~ sum_surprisal + (sum_surprisal | participant_id) + (1 | item),
    data=well_formed %>% filter(
      language == current_language
    ) %>%
      mutate(
        sum_surprisal = scale(sum_surprisal)[,1]
      ),
    control = lmerControl(
      optimizer ="bobyqa",
      optCtrl = list(maxfun = 1e5)
    )
  )
})

names(surprisal_effects) <- setdiff(unique(well_formed$language), "en_exp")

## en_exp manipulated length within item, so here surprisal varies within the 72
## base items and a by-item surprisal slope is estimable (see the note above on
## item ids).
surprisal_effects[["en_exp"]] <- lmer(
  response ~ sum_surprisal +
    (sum_surprisal | participant_id) +
    (sum_surprisal | item),
  data=well_formed %>% filter(
    language == "en_exp"
  ) %>%
    mutate(
      item = as.factor(gsub("[sml]$", "", item)),
      sum_surprisal = scale(sum_surprisal)[,1]
    ),
  control = lmerControl(
    optimizer ="bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)

lapply(surprisal_effects, summary)

lapply(names(surprisal_effects), function(current_language) {
  paste0(current_language,
         ifelse(
           summary(surprisal_effects[[current_language]])$coefficients[2,1] > 0,
           " + ",
           " - "
         ),
         ifelse(
           summary(surprisal_effects[[current_language]])$coefficients[2,5] < .05,
           "<.05",
           ">.05"
         )

  ) %>%
  print()
})

error_free_surprisal_slopes <- bind_rows(lapply(names(surprisal_effects), function(current_language) {
  coefficients <- summary(surprisal_effects[[current_language]])$coefficients
  tibble(
    experiment = current_language,
    beta = coefficients["sum_surprisal", "Estimate"],
    se = coefficients["sum_surprisal", "Std. Error"],
    p = coefficients["sum_surprisal", "Pr(>|t|)"]
  )
}))

write_csv(error_free_surprisal_slopes, "outputs/error_free_surprisal_slopes.csv")

surprisal_betas <- sapply(surprisal_effects,function(m) fixef(m)[2])
surprisal_ttest <- t.test(surprisal_betas)
surprisal_ttest

write_csv(
  tibble(
    mean = unname(surprisal_ttest$estimate),
    ci_low = surprisal_ttest$conf.int[1],
    ci_high = surprisal_ttest$conf.int[2],
    t = unname(surprisal_ttest$statistic),
    df = unname(surprisal_ttest$parameter),
    p = surprisal_ttest$p.value
  ),
  "outputs/error_free_surprisal_ttest.csv"
)

### ---- Generalized Additive Models: by-items ----
crit_by_items <- crit %>%
  mutate(
    weighted_dl_distance_surprisal = dl_distance / sum_surprisal,
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
  group_by(language, item, n_words, sum_surprisal, dl_distance,
           weighted_dl_distance_surprisal, grammatical, condition,
           log_p_sentence, log_p_unigram) %>%
  summarize(
    response_num = mean(response_num)
  ) %>%
  mutate(
    grammatical_numeric = ifelse(grammatical == "grammatical", .5, -.5)
  ) %>%
  ungroup() %>%
  filter(
    dl_distance <= 3
  )

### ---- Compute SLOR ----
## SLOR (Lau, Clark & Lappin 2017) is the sentence's log probability, corrected
## for how frequent its words are, divided by how long it is:
##
##     SLOR = ( log p_model(sentence) - log p_unigram(sentence) ) / length
##
## Subtracting the unigram term removes the penalty a sentence pays simply for
## containing rare words; dividing by length removes the penalty it pays for
## being long. Following its published definition, length is measured in words.
crit_by_items <- crit_by_items %>%
  mutate(
    slor = (log_p_sentence - log_p_unigram) / n_words
  )

### ---- Compute MORCELA ----
## MORCELA (Tjuatja, Neubig, Linzen & Hao 2025) is SLOR with the two corrections
## allowed to be weaker or stronger than SLOR assumes:
##
##     MORCELA = ( log p_model(sentence) - beta * log p_unigram(sentence) + gamma ) / length
##
## SLOR is the special case beta = 1, gamma = 0. The two weights are learned from
## the acceptability ratings themselves, by the recipe in that paper: regress the
## ratings on the three terms below and read the weights off the coefficients.
##
## Note that this fits beta and gamma on the very ratings the score is then used
## to predict. That is what the MORCELA paper does, and it is deliberately
## generous to MORCELA -- we want to compare against the strongest version of it.
crit_by_items <- crit_by_items %>%
  mutate(
    log_p_sentence_over_length = log_p_sentence / n_words,
    log_p_unigram_over_length  = log_p_unigram / n_words,
    one_over_length            = 1 / n_words
  )

morcela_weights <- function(model_data) {
  morcela_regression <- lm(
    response_num ~ log_p_sentence_over_length + log_p_unigram_over_length + one_over_length,
    data = model_data
  )
  weights <- coef(morcela_regression)
  c(
    beta  = -weights[["log_p_unigram_over_length"]] / weights[["log_p_sentence_over_length"]],
    gamma =  weights[["one_over_length"]]           / weights[["log_p_sentence_over_length"]]
  )
}

morcela_score <- function(model_data, weights) {
  model_data$log_p_sentence_over_length -
    weights[["beta"]] * model_data$log_p_unigram_over_length +
    weights[["gamma"]] * model_data$one_over_length
}

### ---- Prepare One Experiment's Data ----
## Used both when fitting below and when plotting further down, so that the two
## can never drift apart. MORCELA's weights are learned separately for each
## experiment, since there is no reason a single beta and gamma should suit all
## of them.
language_data <- function(language_iter) {
  language_specific <- crit_by_items %>%
    filter(language == language_iter) %>% droplevels()

  language_specific$morcela <- morcela_score(
    language_specific,
    morcela_weights(language_specific)
  )

  if (language_iter == "en_exp") {
    ## In this experiment each of the 72 base items was written in 3 lengths, so
    ## the item id encodes the length version (e.g. en_experimental_62m).
    ## Stripping the trailing s/m/l recovers the 72 base items, and surprisal then
    ## varies WITHIN item, which is what makes the by-item slopes identifiable.
    language_specific <- language_specific %>%
      mutate(
        item = as.factor(gsub("[slm]$", "", item))
      )
  }

  language_specific
}

### ---- Prepare the Data Pooled Over All Experiments ----
## The pooled fits include all nine experiments: the eight languages and the
## controlled English experiment. The language random effects are then really
## experiment random effects, with English and controlled English as two separate
## levels. As in the by-experiment fits, the controlled English item ids lose
## their length suffix, so the three length versions of each base item share one
## item intercept.
##
## MORCELA's two weights are estimated once here, on the pooled data, rather than
## separately per experiment. This is one model of all nine experiments, so it
## gets one beta and one gamma.
crit_by_items_all_experiments <- crit_by_items %>%
  mutate(
    item = as.factor(ifelse(language == "en_exp",
                            gsub("[slm]$", "", item), as.character(item)))
  ) %>%
  droplevels()

all_experiments_morcela_weights <- morcela_weights(crit_by_items_all_experiments)

crit_by_items_all_experiments$morcela <- morcela_score(
  crit_by_items_all_experiments,
  all_experiments_morcela_weights
)

## MORCELA's weights, pooled and for each experiment
morcela_weights_table <- bind_rows(
  tibble(
    experiment = "all_experiments",
    beta = all_experiments_morcela_weights[["beta"]],
    gamma = all_experiments_morcela_weights[["gamma"]]
  ),
  bind_rows(lapply(levels(crit_by_items$language), function(language_iter) {
    weights <- morcela_weights(crit_by_items %>% filter(language == language_iter))
    tibble(
      experiment = language_iter,
      beta = weights[["beta"]],
      gamma = weights[["gamma"]]
    )
  }))
)

morcela_weights_table

write_csv(morcela_weights_table, "outputs/morcela_weights.csv")

## ---- Analysis By Experiment ----
## The five models, fit separately to each experiment:
##   - Binary: grammatical vs. ungrammatical
##   - Gradient: a decreasing smooth of the number of errors (DL distance)
##   - Communication: a decreasing smooth of the errors per bit of information
##     (DL distance / surprisal of the intended sentence) plus a decreasing
##     smooth of that surprisal
##   - SLOR and MORCELA: an increasing smooth of the score. A higher score means
##     the sentence is more probable than its words alone would predict, and so
##     should be judged more acceptable. Fitting the decreasing version instead
##     collapses the smooth to a flat line in every language, because the
##     relationship really does go up.
if (file.exists("models/main_surprisal_by_language_5models.RData")) {
  load("models/main_surprisal_by_language_5models.RData")
} else {

  ## Fit a model for each experiment. Only the fitted models are cached; the
  ## log-likelihood table is built further below with marginal_logLik().
  models_by_language <- list()
  for (language_iter in levels(crit_by_items$language)) {
    print(
      paste0("Fitting models for ", language_iter)
    )

    ## Get specific data
    crit_by_items_language_specific <- language_data(language_iter)

    if (language_iter == "en_exp") {
      ## Fit models
      grammatical <- scam(
        response_num ~ grammatical_numeric +
          s(item, bs="re") + s(item, grammatical_numeric, bs="re"),
        data=crit_by_items_language_specific
      )

      distance <- scam(
        response_num ~ s(dl_distance, bs="mpd") +
          s(item, bs="re") + s(item, dl_distance, bs="re"),
        data=crit_by_items_language_specific
      )

      weighted_distance <- scam(
        response_num ~ s(weighted_dl_distance_surprisal, bs="mpd") + s(sum_surprisal, bs="mpd") +
          s(item, bs="re") + s(item, weighted_dl_distance_surprisal, bs="re") +
          s(item, sum_surprisal, bs="re"),
        data=crit_by_items_language_specific
      )

      slor <- scam(
        response_num ~ s(slor, bs="mpi") +
          s(item, bs="re") + s(item, slor, bs="re"),
        data=crit_by_items_language_specific
      )

      morcela <- scam(
        response_num ~ s(morcela, bs="mpi") +
          s(item, bs="re") + s(item, morcela, bs="re"),
        data=crit_by_items_language_specific
      )
    } else {
      ## Fit models
      grammatical <- scam(
        response_num ~ grammatical_numeric +
          s(item, bs="re") ,
        data=crit_by_items_language_specific
      )

      distance <- scam(
        response_num ~ s(dl_distance, bs="mpd") +
          s(item, bs="re"),
        data=crit_by_items_language_specific
      )

      weighted_distance <- scam(
        response_num ~ s(weighted_dl_distance_surprisal, bs="mpd") + s(sum_surprisal, bs="mpd") +
          s(item, bs="re"),
        data=crit_by_items_language_specific
      )

      slor <- scam(
        response_num ~ s(slor, bs="mpi") +
          s(item, bs="re"),
        data=crit_by_items_language_specific
      )

      morcela <- scam(
        response_num ~ s(morcela, bs="mpi") +
          s(item, bs="re"),
        data=crit_by_items_language_specific
      )
    }

    models_by_language[[language_iter]] <- list(
      grammatical = grammatical,
      distance = distance,
      weighted_distance = weighted_distance,
      slor = slor,
      morcela = morcela
    )

  }

  save(
    models_by_language,
    file="models/main_surprisal_by_language_5models.RData"
  )
}

## Check the models, whether they were just fit or loaded from disk: all nine
## experiments are there, and the communication model is the surprisal one
stopifnot(
  setequal(names(models_by_language), levels(crit_by_items$language)),
  all(vapply(models_by_language, function(M) {
    "s(sum_surprisal)" %in% smooth_labels(M$weighted_distance)
  }, logical(1)))
)

## Evaluate Convergence (one column per experiment)
sapply(models_by_language, function(M) sapply(M, function(m) m$conv))

## ---- Analysis Across All Experiments ----
## Each model is also fit once to all nine experiments at the same time, with a
## random effect for language (i.e. experiment) on top of the one for item, and
## reported as the "All Experiments" row of the model-comparison figure.
if (file.exists("models/main_surprisal_all_experiments_5models.RData")) {
  load("models/main_surprisal_all_experiments_5models.RData")
} else {

  print("Fitting models for all experiments together")

  mod_grammatical <- scam(
    response_num ~ grammatical_numeric +  # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, grammatical_numeric, bs="re"), # RE
    data=crit_by_items_all_experiments,
    control = scam.control(trace = TRUE)
  )

  mod_distance <- scam(
    response_num ~ s(dl_distance, bs="mpd") + # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, dl_distance, bs="re"), # RE
    data=crit_by_items_all_experiments,
    control = scam.control(trace = TRUE)
  )

  mod_weighted_distance <- scam(
    response_num ~ s(weighted_dl_distance_surprisal, bs="mpd") + s(sum_surprisal, bs="mpd") + # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, weighted_dl_distance_surprisal, bs="re") +
      s(language, sum_surprisal, bs="re"), # RE
    data=crit_by_items_all_experiments,
    control = scam.control(trace = TRUE)
  )

  mod_slor <- scam(
    response_num ~ s(slor, bs="mpi") + # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, slor, bs="re"), # RE
    data=crit_by_items_all_experiments,
    control = scam.control(trace = TRUE)
  )

  mod_morcela <- scam(
    response_num ~ s(morcela, bs="mpi") + # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, morcela, bs="re"), # RE
    data=crit_by_items_all_experiments,
    control = scam.control(trace = TRUE)
  )

  save(
    mod_grammatical,
    mod_distance,
    mod_weighted_distance,
    mod_slor,
    mod_morcela,
    file="models/main_surprisal_all_experiments_5models.RData"
  )

}

## Check the models, whether they were just fit or loaded from disk: all nine
## experiments are pooled, and the communication model is the surprisal one
stopifnot(
  setequal(levels(mod_weighted_distance$model$language), levels(crit_by_items$language)),
  "s(sum_surprisal)" %in% smooth_labels(mod_weighted_distance)
)

## Evaluate Convergence
mod_grammatical$conv
mod_distance$conv
mod_weighted_distance$conv
mod_slor$conv
mod_morcela$conv

## ---- Model Comparison ----
### ---- Log-likelihoods by experiment ----
## Build the by-experiment log-likelihood table using the MARGINAL log-likelihood
## (integrating the item random effects out) rather than the conditional value
## logLik() returns. scam's own logLik()/AIC()/BIC() are conditional on the
## fitted random effects; use marginal_logLik() throughout.
stats_models_by_language <- do.call(rbind, lapply(names(models_by_language), function(lg) {
  M <- models_by_language[[lg]]
  tibble(experiment = lg, language = lg, model = names(M),
         log_likelihood = vapply(M, marginal_logLik, numeric(1)))
})) %>%
  mutate(
    language = case_when(
      language == "da" ~ "Danish",
      language == "de" ~ "German",
      language == "en" ~ "English",
      language == "hi" ~ "Hindi",
      language == "fr" ~ "French",
      language == "ko" ~ "Korean",
      language == "ru" ~ "Russian",
      language == "zh" ~ "Mandarin",
      language == "en_exp" ~ "English\n(controlled)"
    ),
    language = factor(
      language,
      levels=c("English\n(controlled)", "Russian", "Mandarin",
               "Korean", "Hindi", "German", "French", "English", "Danish")
    ),
    model = case_when(
      model == "grammatical" ~ "Binary",
      model == "distance" ~ "Gradient",
      model == "weighted_distance" ~ "Communication",
      model == "slor" ~ "SLOR",
      model == "morcela" ~ "MORCELA"
    ) %>%
      factor(levels=c("Binary", "Gradient", "Communication", "SLOR", "MORCELA"))
  )

## The three models of the paper keep their colours, shapes and legend labels,
## and SLOR and MORCELA are added. SLOR and MORCELA get filled shapes (23 and 25
## take both a colour and a fill), so all five markers are solid.
MODEL_LEVELS <- c("Binary", "Gradient", "Communication", "SLOR", "MORCELA")
MODEL_LABELS <- c("Binary\nGrammaticality", "Gradient\nGrammaticality",
                  "Communication", "SLOR", "MORCELA")
MODEL_COLORS <- c(
  "Binary"        = "#EB7AA3",
  "Gradient"      = "#1E88E5",
  "Communication" = "#004D40",
  "SLOR"          = "#E08A2E",
  "MORCELA"       = "#A31F34"
)
MODEL_SHAPES <- c(
  "Binary"        = 16,
  "Gradient"      = 17,
  "Communication" = 15,
  "SLOR"          = 23,
  "MORCELA"       = 25
)

### ---- Log-likelihoods across all experiments ----
main_likelihoods <- tibble(
  model = c("Binary", "Gradient", "Communication", "SLOR", "MORCELA"),
  log_likelihood = c(
    marginal_logLik(mod_grammatical),
    marginal_logLik(mod_distance),
    marginal_logLik(mod_weighted_distance),
    marginal_logLik(mod_slor),
    marginal_logLik(mod_morcela)
  )
) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood,
    model = factor(model, levels=MODEL_LEVELS)
  )

main_likelihoods

write_csv(main_likelihoods, "outputs/likelihoods_all_experiments.csv")

p_main_likelihoods <- ggplot(main_likelihoods) +
  geom_point(
    aes(
      y="All Experiments",
      x=likelihood_diff,
      color=model,
      fill=model,
      shape=model
    ),
    size=3
  ) +
  scale_color_manual(
    name="",
    values=MODEL_COLORS,
    labels=MODEL_LABELS
  ) +
  scale_fill_manual(values=MODEL_COLORS) +
  scale_shape_manual(values=MODEL_SHAPES) +
  scale_x_continuous(
    name = ""
  ) +
  ## Five entries do not fit on one row at this width, so the paper's three
  ## models go on the first row and SLOR and MORCELA on the second.
  guides(color = guide_legend(nrow = 2, byrow = TRUE), shape = "none", fill = "none") +
  theme(legend.position="top") +
  scale_y_discrete(name="")

### ---- Figure based on individual-language likelihoods ----
stats_models_by_language <- stats_models_by_language %>%
  group_by(language) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood
  )

stats_models_by_language

write_csv(
  stats_models_by_language %>%
    ungroup() %>%
    select(experiment, model, log_likelihood, max_likelihood, likelihood_diff),
  "outputs/likelihoods_by_language.csv"
)

p_language_likelihoods <- ggplot(stats_models_by_language) +
  geom_point(
    aes(
      y=language,
      x=likelihood_diff,
      color=model,
      fill=model,
      shape=model
    ),
    size=3
  ) +
  scale_color_manual(
    name="",
    values=MODEL_COLORS
  ) +
  scale_fill_manual(values=MODEL_COLORS) +
  scale_shape_manual(values=MODEL_SHAPES) +
  scale_x_continuous(
    name = "Log-Likelihood Difference"
  ) +
  theme(legend.position="none") +
  scale_y_discrete(name="")

## Stack the combined "All Experiments" panel (top) over the by-language panel
## (bottom), matching the figure caption.
likelihood_differences_figure <- p_main_likelihoods / p_language_likelihoods +
  plot_layout(heights = c(1, 3))

ggsave("figures/likelihood_differences.png",
       plot = likelihood_differences_figure, width = 8, height = 8)

### ---- Smooths of the pooled models ----
## Save the figures
## The rug ticks are jittered at random; a fixed seed makes the figure reproducible
set.seed(2026)
png("figures/mod_gradient_dl.png", width=6, height=6, units="in", res=300)
plot(mod_distance, select = 1,
     shade = TRUE,            # optional: shaded CI instead of dashed lines
     shade.col = "lightblue", # color for shading
     main = "Gradient Grammaticality", # title
     xlab = "DL distance",    # x-axis label
     ylab = "Acceptability (centered)", # y-axis label
     rug = TRUE,              # keep the rug
     jit = TRUE               # jitter the rug ticks
)
dev.off()

## The rug ticks are jittered at random; a fixed seed makes the figure reproducible
set.seed(2026)
png("figures/mod_communication_dlweighted_surprisal.png", width=6, height=6, units="in", res=300)
plot(mod_weighted_distance, select = 1,
     shade = TRUE,            # optional: shaded CI instead of dashed lines
     shade.col = "lightblue", # color for shading
     main = "Communication", # title
     xlab = "DL distance / surprisal (bits)",    # x-axis label
     ylab = "Acceptability (centered)", # y-axis label
     rug = TRUE,              # keep the rug
     jit = TRUE               # jitter the rug ticks
)
dev.off()

## The rug ticks are jittered at random; a fixed seed makes the figure reproducible
set.seed(2026)
png("figures/mod_communication_sum_surprisal.png", width=6, height=6, units="in", res=300)
plot(mod_weighted_distance, select = 2,
         shade = TRUE,            # optional: shaded CI instead of dashed lines
         shade.col = "lightblue", # color for shading
         main = "Communication", # title
         xlab = "Sentence surprisal (bits)",    # x-axis label
         ylab = "Acceptability (centered)", # y-axis label
         rug = TRUE,              # keep the rug
         jit = TRUE               # jitter the rug ticks
)
dev.off()

## ---- Plotting Helpers ----
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

## ---- Main Figure: Observed vs. Predicted Acceptability ----
## Observed against predicted acceptability, one point per item and condition.
## One row for the models pooled over all experiments, then one row per
## experiment; one column per model.
OBSERVED_VS_PREDICTED_MODELS <- c(
  slor = "SLOR",
  morcela = "MORCELA",
  weighted_distance = "Communication",
  distance = "Gradient",
  grammatical = "Binary"
)

## The row of each panel: "All Experiments" first, then the experiments in the
## order label_language() gives them
label_row <- function(experiment_code) {
  factor(
    ifelse(
      experiment_code == "all_experiments",
      "All Experiments",
      as.character(label_language(experiment_code))
    ),
    levels = c("All Experiments", levels(label_language("da")))
  )
}

## Predicted acceptability with every random-effect term left OUT (item
## intercepts, the by-item slopes of the controlled English experiment, and in
## the pooled fits the language intercepts and slopes), so each prediction comes
## from the model's predictors alone. Leaving them in would let every model
## borrow each item's own average rating, and all five would hug the diagonal
## regardless of how good their predictors are.
observed_vs_predicted_data <- function(model_list, model_data, experiment_code) {
  bind_rows(lapply(names(OBSERVED_VS_PREDICTED_MODELS), function(model_type) {
    m <- model_list[[model_type]]
    tibble(
      experiment = experiment_code,
      language = as.character(model_data$language),
      item = as.character(model_data$item),
      condition = model_data$condition,
      n_words = model_data$n_words,
      sum_surprisal = model_data$sum_surprisal,
      dl_distance = model_data$dl_distance,
      observed = model_data$response_num,
      model = OBSERVED_VS_PREDICTED_MODELS[[model_type]],
      predicted = as.numeric(
        predict(m, newdata = model_data, exclude = random_effect_labels(m))
      )
    )
  }))
}

all_predictions <- bind_rows(
  observed_vs_predicted_data(
    list(
      grammatical = mod_grammatical,
      distance = mod_distance,
      weighted_distance = mod_weighted_distance,
      slor = mod_slor,
      morcela = mod_morcela
    ),
    crit_by_items_all_experiments,
    "all_experiments"
  ),
  bind_rows(lapply(names(models_by_language), function(language_iter) {
    observed_vs_predicted_data(
      models_by_language[[language_iter]],
      language_data(language_iter),
      language_iter
    )
  }))
) %>%
  mutate(
    model = factor(model, levels = OBSERVED_VS_PREDICTED_MODELS)
  )

write_csv(all_predictions, "outputs/observed_vs_predicted.csv")

all_predictions <- all_predictions %>%
  mutate(
    row_label = label_row(experiment)
  )

## Points are drawn for distances 0, 1 and 3. A few three-swap sentences come out
## at distance 2; they are in the models, but not drawn.
print(paste0("Item-condition cells at distance 2, not drawn: ",
             sum(all_predictions$dl_distance == 2 &
                   all_predictions$experiment != "all_experiments") /
               length(OBSERVED_VS_PREDICTED_MODELS)))

## Colour shows the DL distance, in the colours of the curve figure.
## Transparency shows sentence surprisal: low-surprisal sentences are faint,
## high-surprisal ones are solid. One transparency scale is shared by every
## panel, so a given shade means the same number of bits in every row. Higher-
## surprisal sentences are drawn last, so they sit on top.
plot_data <- all_predictions %>%
  filter(dl_distance %in% c(0, 1, 3)) %>%
  mutate(dl_label = factor(dl_distance, levels = c(0, 1, 3))) %>%
  arrange(sum_surprisal)

## The same range on both axes and in every panel, wide enough for every point,
## so that the dashed diagonal is always the line of perfect prediction
axis_range <- c(
  min(1, floor(min(plot_data$predicted) * 2) / 2),
  max(7, ceiling(max(plot_data$predicted) * 2) / 2)
)

## Pearson correlation between observed and predicted acceptability in each
## panel, over the points shown
panel_correlations <- plot_data %>%
  group_by(experiment, row_label, model) %>%
  summarize(r = cor(observed, predicted), n = n(), .groups = "drop") %>%
  mutate(label = paste0("r = ", sprintf("%.2f", r)))

print(panel_correlations %>%
        select(model, row_label, r) %>%
        pivot_wider(names_from = row_label, values_from = r),
      width = Inf)

write_csv(panel_correlations %>% select(experiment, model, r, n),
          "outputs/observed_vs_predicted_correlations.csv")

## MORCELA's two fitted weights for each row, printed in its panels. They are
## written as plotmath (parse = TRUE below) rather than as Greek characters,
## because R in a non-UTF-8 locale mangles non-ASCII text.
morcela_labels <- morcela_weights_table %>%
  mutate(
    row_label = label_row(experiment),
    model = factor("MORCELA", levels = levels(all_predictions$model)),
    label = paste0("list(beta == '", sprintf("%.2f", beta),
                   "', gamma == '", sprintf("%.1f", gamma), "')")
  )

p_observed_vs_predicted <- ggplot(plot_data, aes(x = predicted, y = observed)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey60") +
  geom_point(aes(colour = dl_label, alpha = sum_surprisal), size = 1.8) +
  scale_colour_manual(
    name = "DL distance",
    values = c("0" = "#EB7AA3", "1" = "#1E88E5", "3" = "#FFC107"),
    guide = guide_legend(override.aes = list(size = 4))
  ) +
  scale_alpha_continuous(
    name = "Sentence surprisal\n(bits)",
    range = c(0.2, 1)
  ) +
  geom_text(
    data = panel_correlations,
    aes(x = axis_range[1] + 0.1, y = axis_range[2], label = label),
    hjust = 0,
    vjust = 1,
    size = 5,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = morcela_labels,
    aes(x = axis_range[2], y = axis_range[1], label = label),
    parse = TRUE,
    hjust = 1,
    vjust = 0,
    size = 5,
    inherit.aes = FALSE
  ) +
  facet_grid(row_label ~ model) +
  coord_fixed(ratio = 1, xlim = axis_range, ylim = axis_range) +
  scale_x_continuous(breaks = seq(1, 7, 2)) +
  scale_y_continuous(breaks = seq(1, 7, 2)) +
  labs(
    x = "Predicted acceptability",
    y = "Observed mean acceptability",
    caption = "Predictions from each model's predictors alone (item and language random effects excluded).\nDashed line: perfect prediction. r: Pearson correlation over the points shown."
  ) +
  theme(plot.caption = element_text(size = 14, colour = "grey40"))

ggsave(
  p_observed_vs_predicted,
  filename = "figures/observed_vs_predicted.png",
  width    = 18,
  height   = 32.5,
  units    = "in",
  dpi      = 200
)

## Save as pdf
ggsave(
  p_observed_vs_predicted,
  filename = "figures/observed_vs_predicted.pdf",
  width    = 18,
  height   = 32.5,
  units    = "in"
)

## ---- SI Figure: Model Curves over Surprisal ----
## Human judgments and the predictions of the binary, gradient and communication
## models, drawn over sentence surprisal. Written to the SI figures folder.
SI_FIGURES <- "../supporting_information/figures"
dir.create(SI_FIGURES, showWarnings = FALSE)

#### ---- Helper Functions
## The binary, gradient and communication models are functions of surprisal and
## of the number of errors, so we can ask them what they predict for any
## combination of the two, whether or not a sentence like that was in the
## materials. That is what this grid is for.
get_prediction_grid <- function(
    lang,
    model_type,
    sum_surprisal=c(50, 400),
    dl_distance=c(0, 1, 3)
) {
  expand.grid(
    sum_surprisal = seq(sum_surprisal[1] - .25, sum_surprisal[2] + .25, 1),
    dl_distance = dl_distance
  ) %>%
    mutate(
      weighted_dl_distance_surprisal = dl_distance / sum_surprisal,
      grammatical_numeric = ifelse(dl_distance == 0, .5, -.5),
      item=factor(models_by_language[[lang]][[model_type]]$model$item[1], levels=levels(models_by_language[[lang]][[model_type]]$model$item))
    )
}


### ---- Baseline data for plotting (item means) ----
## Same aggregation as crit_by_items above, but deliberately WITHOUT the
## filter(dl_distance <= 3): the scatter shows the shuffled items as well, in the
## grey "4+" bucket, as context. The models (and hence the fitted curves drawn
## over these points) are still fit on distances up to 3 only.
plot_data_by_item <- crit %>%
  mutate(
    weighted_dl_distance_surprisal = dl_distance / sum_surprisal,
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
  group_by(language, item, sum_surprisal, dl_distance, weighted_dl_distance_surprisal,
           grammatical, condition) %>%
  summarize(
    response_num = mean(response_num)
  ) %>%
  mutate(
    grammatical_numeric = ifelse(grammatical == "grammatical", .5, -.5)
  ) %>%
  dplyr::group_by(language, item, sum_surprisal, dl_distance) %>%
  dplyr::summarise(response = mean(response_num), .groups = "drop") %>%
  dplyr::mutate(
    dl_distance_bucket = bucket_dl(dl_distance),
    language           = label_language(language)
  )



### ---- Baseline scatterplot (item means) ----
p_baseline <- ggplot(
  plot_data_by_item,
  aes(x = sum_surprisal, y = response, color = dl_distance_bucket)
) +
  geom_point(alpha = 0.3) +
  facet_wrap(~language, ncol = 1, strip.position = "left") +
  labs(x="Sentence Surprisal (bits)", y="Rating") +
  theme(
          legend.position = "top",
          plot.title = element_text(hjust = 0.5),

          # put facet strips on the outside + rotate text
          strip.placement = "outside",
          strip.text.y.left = element_text(angle = 90)
      ) +
  scale_y_continuous(
    limits = YLIM,
    breaks = NULL,
    labels = NULL
  ) +
  scale_dl_color()

### ---- GAM smooth over item means ----
p_gam <- p_baseline +
  geom_smooth(
    aes(x = sum_surprisal, y = response, color = dl_distance_bucket),
    method  = "gam",
    formula = y ~ s(x, bs = "tp"),
    se = TRUE,
    linewidth = 1.5,
    show.legend = FALSE
  ) +
  theme(
    legend.position = "none"
  ) +
  ggtitle("Human Judgments") +
  scale_y_continuous(
    limits = YLIM,
    breaks = seq(1, 7, 2),
    labels = seq(1, 7, 2),
    name = "Acceptability Ratings (1-7 scale)"
  ) +
  ggplot2::guides(color = "none")

ggsave(p_gam, filename = file.path(SI_FIGURES, "curves_surprisal_gam.png"),
       width    = 6,
       height   = 16,
       units    = "in",
       dpi      = 300)

### ---- Model-based plots ----
## Predictions leave out every random effect (random_effect_labels()), so that
## the controlled English curves do not carry the by-item slopes of whichever
## item the grid happens to name.
#### ---- Grammatical model ----
#### ---- Predictions from grammatical model ----
grammatical_predictions <- setNames(
    lapply(names(models_by_language), function(lang) {
      predict(
        models_by_language[[lang]]$grammatical,
        newdata = get_prediction_grid(lang, "grammatical", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal)),
        type = "response",
        se.fit = TRUE,
        exclude = random_effect_labels(models_by_language[[lang]]$grammatical)
      )
    }
  ),
  names(models_by_language)
)

## Combine list to tibble
grammatical_plot_data <- bind_rows(
  lapply(names(grammatical_predictions), function(lang) {
    preds <- grammatical_predictions[[lang]]
    grid  <- get_prediction_grid(lang, "grammatical", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal))

    tibble(
      language = label_language(lang),
      sum_surprisal = grid$sum_surprisal,
      dl_distance_bucket = bucket_dl(grid$dl_distance),
      response = preds$fit,
      se = preds$se.fit
    )
  })
)

##### ---- Plot grammatical model ----
p_grammatical <- p_baseline +
  ggnewscale::new_scale_color() +   # <- reset color scale AFTER baseline
  geom_ribbon(
    data=grammatical_plot_data,
    aes(
      x = sum_surprisal,
      ymin = response - 2 * se,
      ymax = response + 2 * se,
      group = dl_distance_bucket
    ),
    alpha = 0.25
    ) +
  geom_line(
    data=grammatical_plot_data,
    aes(
      x = sum_surprisal,
      y = response,
      color = dl_distance_bucket
    ),
    linewidth = 1.5
  ) +
  facet_wrap(~language, ncol = 1) +
  labs(
    x = "Sentence Surprisal (bits)",
    y = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  ggtitle("Binary Grammaticality Model") +
  scale_y_continuous(
    limits = YLIM,
    breaks = NULL,
    labels = NULL
  ) +
  scale_color_manual(
    values = c(
      "0" = "#EB7AA3",
      "1-2" = "#06402B",
      "3" = "#06402B",
      "4+" = "#06402B"
    )
  ) +
  ggplot2::guides(color = "none")

p_grammatical

ggsave(
  p_grammatical,
  filename = file.path(SI_FIGURES, "curves_surprisal_grammatical.png"),
  width    = 6,
  height   = 16,
  units    = "in",
  dpi      = 300
)

#### ---- Distance model ----
#### ---- Predictions from distance model ----
distance_predictions <- setNames(
  lapply(names(models_by_language), function(lang) {
    predict(
      models_by_language[[lang]]$distance,
      newdata = get_prediction_grid(lang, "distance", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal)),
      type = "response",
      se.fit = TRUE,
      exclude = random_effect_labels(models_by_language[[lang]]$distance)
    )
  }
  ),
  names(models_by_language)
)

## Combine list to tibble
distance_plot_data <- bind_rows(
  lapply(names(distance_predictions), function(lang) {
    preds <- distance_predictions[[lang]]
    grid  <- get_prediction_grid(lang, "distance", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal))

    tibble(
      language = label_language(lang),
      sum_surprisal = grid$sum_surprisal,
      dl_distance_bucket = bucket_dl(grid$dl_distance),
      response = preds$fit,
      se = preds$se.fit
    )
  })
)

##### ---- Plot distance model ----
p_distance <- p_baseline +
  geom_ribbon(
    data=distance_plot_data,
    aes(
      x = sum_surprisal,
      ymin = response - 2 * se,
      ymax = response + 2 * se,
      group = dl_distance_bucket
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(
    data=distance_plot_data,
    aes(
      x = sum_surprisal,
      y = response,
      color = dl_distance_bucket
    ),
    linewidth = 1.5
  ) +
  facet_wrap(~language, ncol = 1) +
  labs(
    x = "Sentence Surprisal (bits)",
    y = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  ggtitle("Gradient Grammaticality Model") +
  scale_y_continuous(
    limits = YLIM,
    breaks = NULL,
    labels = NULL
  ) +
  ggplot2::guides(color = "none")

p_distance

ggsave(
  p_distance,
  filename = file.path(SI_FIGURES, "curves_surprisal_distance.png"),
  width    = 6,
  height   = 16,
  units    = "in",
  dpi      = 300
)

#### ---- Weighted distance model ----
#### ---- Predictions from weighted distance model ----
weighted_distance_predictions <- setNames(
  lapply(names(models_by_language), function(lang) {
    predict(
      models_by_language[[lang]]$weighted_distance,
      newdata = get_prediction_grid(lang, "weighted_distance", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal)),
      type = "response",
      se.fit = TRUE,
      exclude = random_effect_labels(models_by_language[[lang]]$weighted_distance)
    )
  }
  ),
  names(models_by_language)
)

## Combine list to tibble
weighted_distance_plot_data <- bind_rows(
  lapply(names(weighted_distance_predictions), function(lang) {
    preds <- weighted_distance_predictions[[lang]]
    grid  <- get_prediction_grid(lang, "weighted_distance", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal))

    tibble(
      language = label_language(lang),
      sum_surprisal = grid$sum_surprisal,
      dl_distance_bucket = bucket_dl(grid$dl_distance),
      response = preds$fit,
      se = preds$se.fit
    )
  })
)

##### ---- Plot weighted distance (communication-based) model ----
p_weighted_distance <- p_baseline +
  geom_ribbon(
    data=weighted_distance_plot_data %>%
      mutate(
        lower = response - 2*se,
        upper = response + 2*se,
        upper = ifelse(upper > 7.5, 7.5, upper),
      ),
    aes(
      x = sum_surprisal,
      ymin = lower,
      ymax = upper,
      group = dl_distance_bucket
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(
    data=weighted_distance_plot_data %>%
      mutate(
        response = ifelse(response > 7.5, 7.5, response)
      ),
    aes(
      x = sum_surprisal,
      y = response,
      color = dl_distance_bucket
    ),
    linewidth = 1.5
  ) +
  facet_wrap(~language, ncol = 1) +
  labs(
    x = "Sentence Surprisal (bits)",
    y = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  ggtitle("Communication Model") +
  scale_y_continuous(
    limits = YLIM,
    breaks = NULL,
    labels = NULL
  ) +
  ggplot2::guides(color = "none")

p_weighted_distance

ggsave(
  p_weighted_distance,
  filename = file.path(SI_FIGURES, "curves_surprisal_weighted_distance.png"),
  width    = 6,
  height   = 16,
  units    = "in",
  dpi      = 300
)

### ---- Combine plots into the curve figure ----
curves_surprisal_figure = (p_gam |  p_weighted_distance | p_distance |p_grammatical) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
        )

curves_surprisal_figure

ggsave(
  filename = file.path(SI_FIGURES, "curves_surprisal.png"),
  plot     = curves_surprisal_figure,
  width    = 25,
  height   = 20,
  units    = "in",
  dpi      = 300
)

## Save as pdf
ggsave(
  filename = file.path(SI_FIGURES, "curves_surprisal.pdf"),
  plot     = curves_surprisal_figure,
  width    = 25,
  height   = 20,
  units    = "in"
)

print("Done.")
