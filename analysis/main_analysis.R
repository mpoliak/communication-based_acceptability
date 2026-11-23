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
## outputs/ holds the tables this script writes (statistics, log-likelihoods and
## model predictions).
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
## Kept as a function: written out, this Filter()/vapply() pair would be repeated
## in all four predict() calls below.
random_effect_labels <- function(m) {
  vapply(
    Filter(function(s) inherits(s, "random.effect"), m$smooth),
    function(s) s$label,
    character(1)
  )
}

## Full language names, for the tables and figures. The levels give the rows
## from top to bottom; the likelihood figure lists them in reverse, as its y axis
## runs from the bottom up.
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
  levels=c("Danish", "English", "French", "German", "Hindi", "Korean",
           "Mandarin", "Russian", "English\n(controlled)")
  )
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

## Every participant took part in one experiment only
stopifnot(!anyDuplicated(distinct(crit, participant_id, language)$participant_id))

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
## We use the model-initial unigram throughout. wordfreq is also in the file, but
## no model here uses it (SI section 4 compares unigram models in English only);
## it does not know some of the Danish and Korean word tokens, so it cannot be
## the default in a cross-linguistic analysis.
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
## same in all of an item's conditions, like its length (n_words).
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
## an item's conditions: a few controlled English items were never rated in the
## original condition itself.
original_sentences <- crit %>%
  distinct(language, item, n_words, sum_surprisal)

stopifnot(!anyDuplicated(original_sentences[c("language", "item")]))

## Token count vs. summed surprisal: Pearson correlation in each experiment, and
## in the eight language experiments together (the last row; controlled English
## is left out of the pool)
surprisal_length_correlations <- original_sentences %>%
  bind_rows(original_sentences %>%
              filter(language != "en_exp") %>%
              mutate(language = "all_experiments")) %>%
  group_by(experiment = language) %>%
  summarize(
    n_sentences = n(),
    r = cor(n_words, sum_surprisal),
    ci_low = cor.test(n_words, sum_surprisal)$conf.int[1],
    ci_high = cor.test(n_words, sum_surprisal)$conf.int[2],
    p = cor.test(n_words, sum_surprisal)$p.value
  ) %>%
  arrange(experiment == "all_experiments", experiment)

surprisal_length_correlations

write_csv(surprisal_length_correlations, "outputs/surprisal_length_correlations.csv")

## ---- Main Analyses ----
### ---- Linear Model ----
## Every experiment uses the nlminb optimizer. allFit() shows nlminb at the best
## log-likelihood of any optimizer (to within 0.0001) in all nine, while bobyqa
## stops ~27 log-likelihood units short in en_exp. In Danish, nlminb reports
## "singular convergence (7)" at the same optimum that bobyqa finds (identical
## estimates): the fit is singular, like most of these maximal fits.
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
    optimizer ="nlminbwrap",
    optCtrl = list(eval.max = 1e5, iter.max = 1e5)
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
  ## nlminb, as above. With this random-effects structure bobyqa stops ~27
  ## log-likelihood units short of the optimum and fails the convergence check.
  ## The fit is still singular: the maximal structure is more than the data
  ## support.
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
## Item x condition means, averaged over participants: the data of all five
## models and of MORCELA's weights. The DL <= 3 filter is on DL distance, not
## condition: shuffled sentences that came out at DL 2-3 stay in, the
## three-exchange sentence at DL 4 drops out, and one-exchange sentences the
## exchange left unchanged (DL 0) count as error-free.
crit_by_items <- crit %>%
  mutate(
    weighted_dl_distance_surprisal = dl_distance / sum_surprisal,
    grammatical = factor(
      ifelse(dl_distance == 0, "grammatical", "ungrammatical"),
      levels=c("ungrammatical", "grammatical")),
    item = as.factor(item),
    language = as.factor(language),
    response_num = response
  ) %>%
  group_by(language, item, n_words, sum_surprisal, dl_distance,
           weighted_dl_distance_surprisal, grammatical, condition,
           log_p_sentence, log_p_unigram) %>%
  summarize(
    response_num = mean(response_num),
    .groups = "drop"
  ) %>%
  mutate(
    grammatical_numeric = ifelse(grammatical == "grammatical", .5, -.5)
  ) %>%
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
## being long. Length is the number of tokens (space-separated; jiebaR segments
## in Mandarin), as in the paper.
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

## MORCELA's weights, once for each experiment and once for the data pooled over
## the eight languages. Each experiment gets its own weights, since there is no
## reason a single beta and gamma should suit all of them; the pooled fits below
## are one model of the eight language experiments, so they get one beta and
## one gamma. The weights come from the coefficients of the regression:
## beta = -(unigram coefficient) / (sentence coefficient) and
## gamma = (1/length coefficient) / (sentence coefficient).
morcela_weights_table <- bind_rows(
  crit_by_items %>%
    filter(language != "en_exp") %>%
    mutate(experiment = "all_experiments"),
  crit_by_items %>%
    mutate(experiment = as.character(language))
) %>%
  group_by(experiment) %>%
  summarize(
    weights = list(coef(lm(
      response_num ~ log_p_sentence_over_length + log_p_unigram_over_length + one_over_length,
      data = pick(everything())
    ))),
    .groups = "drop"
  ) %>%
  mutate(
    beta  = map_dbl(weights, ~ -.x[["log_p_unigram_over_length"]] / .x[["log_p_sentence_over_length"]]),
    gamma = map_dbl(weights, ~  .x[["one_over_length"]]           / .x[["log_p_sentence_over_length"]])
  ) %>%
  select(experiment, beta, gamma)

morcela_weights_table

write_csv(morcela_weights_table, "outputs/morcela_weights.csv")

### ---- Prepare Each Experiment's Data ----
## Used both when fitting below and when plotting further down, so that the two
## can never drift apart. MORCELA is scored with the experiment's own weights.
data_by_language <- lapply(levels(crit_by_items$language), function(language_iter) {
  weights <- morcela_weights_table %>% filter(experiment == language_iter)
  crit_by_items %>%
    filter(language == language_iter) %>%
    droplevels() %>%
    mutate(
      morcela = log_p_sentence_over_length -
        weights$beta * log_p_unigram_over_length +
        weights$gamma * one_over_length,
      ## In the controlled English experiment, strip the trailing s/m/l to
      ## recover the 72 base items (see the note on item ids at its mixed model
      ## above)
      item = if (language_iter == "en_exp") as.factor(gsub("[sml]$", "", item)) else item
    )
})

names(data_by_language) <- levels(crit_by_items$language)

### ---- Prepare the Data Pooled Over All Languages ----
## The pooled fits include the eight language experiments only. The controlled
## English experiment is left out of the pool and is fit only on its own, in the
## by-experiment fits below. The language random effects then have one level per
## language. The objects and output files keep their all_experiments names.
## MORCELA is scored with the pooled weights.
all_experiments_weights <- morcela_weights_table %>%
  filter(experiment == "all_experiments")

crit_by_items_all_experiments <- crit_by_items %>%
  filter(language != "en_exp") %>%
  droplevels() %>%
  mutate(
    morcela = log_p_sentence_over_length -
      all_experiments_weights$beta * log_p_unigram_over_length +
      all_experiments_weights$gamma * one_over_length
  )

## ---- Analysis By Experiment ----
## The five models, fit separately to each experiment:
##   - Binary: grammatical vs. ungrammatical
##   - Gradient: a decreasing smooth of the number of errors (DL distance)
##   - Communication: a decreasing smooth of the errors per bit of information
##     (DL distance / surprisal of the intended sentence) plus a decreasing
##     smooth of that surprisal
##   - SLOR and MORCELA: an increasing smooth of the score. A higher score means
##     the sentence is more probable than its words alone would predict, and so
##     should be judged more acceptable.
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
    crit_by_items_language_specific <- data_by_language[[language_iter]]

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
    "s(sum_surprisal)" %in% vapply(M$weighted_distance$smooth, function(s) s$label, character(1))
  }, logical(1)))
)

## Evaluate Convergence of scam's inner fitting loop (one column per experiment);
## the smoothing-parameter search is checked below
sapply(models_by_language, function(M) sapply(M, function(m) m$conv))

## ---- Analysis Across All Languages ----
## Each model is also fit once to the eight language experiments at the same
## time (controlled English is left out of the pool), with a random effect for
## language on top of the one for item, and reported as the "All Languages" row
## of the model-comparison figure.
if (file.exists("models/main_surprisal_all_experiments_5models.RData")) {
  load("models/main_surprisal_all_experiments_5models.RData")
} else {

  print("Fitting models for all languages together")

  mod_grammatical <- scam(
    response_num ~ grammatical_numeric +  # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, grammatical_numeric, bs="re"), # RE
    data=crit_by_items_all_experiments
  )

  mod_distance <- scam(
    response_num ~ s(dl_distance, bs="mpd") + # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, dl_distance, bs="re"), # RE
    data=crit_by_items_all_experiments
  )

  mod_weighted_distance <- scam(
    response_num ~ s(weighted_dl_distance_surprisal, bs="mpd") + s(sum_surprisal, bs="mpd") + # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, weighted_dl_distance_surprisal, bs="re") +
      s(language, sum_surprisal, bs="re"), # RE
    data=crit_by_items_all_experiments
  )

  mod_slor <- scam(
    response_num ~ s(slor, bs="mpi") + # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, slor, bs="re"), # RE
    data=crit_by_items_all_experiments
  )

  mod_morcela <- scam(
    response_num ~ s(morcela, bs="mpi") + # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, morcela, bs="re"), # RE
    data=crit_by_items_all_experiments
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

## Check the models, whether they were just fit or loaded from disk: the eight
## language experiments are pooled, without controlled English (so an old cache
## of all nine experiments stops the script here), and the communication model
## is the surprisal one
stopifnot(
  setequal(levels(mod_weighted_distance$model$language),
           setdiff(levels(crit_by_items$language), "en_exp")),
  "s(sum_surprisal)" %in% vapply(mod_weighted_distance$smooth, function(s) s$label, character(1))
)

## Evaluate Convergence of scam's inner fitting loop; the smoothing-parameter
## search is checked below
mod_grammatical$conv
mod_distance$conv
mod_weighted_distance$conv
mod_slor$conv
mod_morcela$conv

## How each fit's search for its smoothing parameters ended, per experiment and
## pooled. scam stops the search when the gradient of its smoothness score is
## small enough (gradtol.bfgs, 1e-6 by default), when it reaches its 200-step
## limit, or when its last step finds no lower score. `tolerance_met` is the
## smallest gradtol.bfgs at which the fit would pass scam's own stopping rule
## (in scam's bfgs_gcv.ubre, which compares the relative gradient with
## 6.0554 * gradtol.bfgs). Not reported in the paper.
all_fits <- c(
  unlist(models_by_language, recursive = FALSE),
  list(
    all_experiments.grammatical       = mod_grammatical,
    all_experiments.distance          = mod_distance,
    all_experiments.weighted_distance = mod_weighted_distance,
    all_experiments.slor              = mod_slor,
    all_experiments.morcela           = mod_morcela
  )
)

search_convergence <- tibble(
  fit = names(all_fits),
  bfgs_convergence = vapply(all_fits, function(m) m$bfgs.info$conv, character(1)),
  tolerance_met = vapply(all_fits, function(m) {
    max(abs(m$bfgs.info$grad)) * max(max(abs(log(m$sp))), 1) /
      max(abs(m$gcv.ubre), 1) / 6.0554
  }, numeric(1))
)

## Three fits stop at the 200-step limit (French MORCELA, Russian SLOR and
## MORCELA), each within about twice the default tolerance (tolerance_met about
## 2e-6, against gradtol.bfgs = 1e-6), so their search had all but ended. The
## closest margin they touch is Russian MORCELA over SLOR (1.78).
print(search_convergence %>% filter(bfgs_convergence != "Full convergence"), n = Inf)

## The fits whose last step found no lower score: refit each with scam's two
## other optimizers, and print how much its marginal log-likelihood changes.
## These are the Danish and Mandarin communication fits, which lead the
## next-best model by 5.2 and 18.8.
for (fit_name in search_convergence$fit[grepl("^Last step failed", search_convergence$bfgs_convergence)]) {
  m <- all_fits[[fit_name]]
  print(paste0(
    fit_name, ": marginal log-likelihood change with nlm ",
    round(marginal_logLik(scam(m$formula, data = m$model, optimizer = "nlm")) - marginal_logLik(m), 3),
    ", with optim ",
    round(marginal_logLik(scam(m$formula, data = m$model, optimizer = "optim")) - marginal_logLik(m), 3)
  ))
}

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
    ## Reversed, so that the likelihood figure lists the languages from the top
    language = fct_rev(label_language(language)),
    model = case_when(
      model == "grammatical" ~ "Binary",
      model == "distance" ~ "Gradient",
      model == "weighted_distance" ~ "Communication",
      model == "slor" ~ "SLOR",
      model == "morcela" ~ "MORCELA"
    ) %>%
      factor(levels=c("Binary", "Gradient", "Communication", "SLOR", "MORCELA"))
  )

## Colours, shapes and legend labels of the five models, as in the SI figures
## (setup.R). SLOR and MORCELA get filled shapes (23 and 25 take both a colour
## and a fill), so all five markers are solid.
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

### ---- Log-likelihoods across all languages ----
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
      y="All Languages",
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
  scale_fill_manual(name="", values=MODEL_COLORS, labels=MODEL_LABELS) +
  scale_shape_manual(name="", values=MODEL_SHAPES, labels=MODEL_LABELS) +
  scale_x_continuous(
    name = ""
  ) +
  ## Five entries do not fit on one row at this width, so the binary, gradient
  ## and communication models go on the first row and SLOR and MORCELA on the
  ## second. The colour, fill and shape scales share a name, labels and guide, so
  ## ggplot merges them into one legend whose keys show each model's marker.
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    fill = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
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

## ---- Figure 4: Model Comparison (Log-likelihoods) ----
## Stack the combined "All Languages" panel (top) over the by-language panel
## (bottom), matching the figure caption. It is drawn at 8 x 8 in, like the SI's
## likelihood figures; printed one column wide (3.42 in), its smallest text (tick
## labels and legend) comes out at about 6 pt, the PNAS minimum.
likelihood_differences_figure <- p_main_likelihoods / p_language_likelihoods +
  plot_layout(heights = c(1, 3))

ggsave("figures/likelihood_differences.png",
       plot = likelihood_differences_figure, width = 8, height = 8)

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

## Figures 1-3 are drawn at their print size, so that their text comes out at
## the size it is printed. This is the script's theme at the 8-pt print size of
## the full-width figures, with axis lines thinned to match.
theme_print_8pt <- theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.key = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5))

## ---- Predictions of the Five Models (Figures 2 and 3) ----
## Observed against predicted acceptability, one point per item and condition,
## for all five models: one row for the models pooled over the eight languages,
## then one row per experiment. The order of the models here is the row order
## of the saved tables.
OBSERVED_VS_PREDICTED_MODELS <- c(
  slor = "SLOR",
  morcela = "MORCELA",
  weighted_distance = "Communication",
  distance = "Gradient",
  grammatical = "Binary"
)

## Each row's fits and data: the pooled fits first, then each experiment
fits_by_row <- c(
  list(all_experiments = list(
    grammatical = mod_grammatical,
    distance = mod_distance,
    weighted_distance = mod_weighted_distance,
    slor = mod_slor,
    morcela = mod_morcela
  )),
  models_by_language
)

data_by_row <- c(
  list(all_experiments = crit_by_items_all_experiments),
  data_by_language
)

## Predicted acceptability with every random-effect term left OUT (item
## intercepts, the by-item slopes of the controlled English experiment, and in
## the pooled fits the language intercepts and slopes), so each prediction comes
## from the model's predictors alone. Leaving them in would let every model
## borrow each item's own average rating, and all five would hug the diagonal
## regardless of how good their predictors are.
all_predictions <- bind_rows(lapply(names(fits_by_row), function(experiment_code) {
  model_data <- data_by_row[[experiment_code]]
  bind_rows(lapply(names(OBSERVED_VS_PREDICTED_MODELS), function(model_type) {
    m <- fits_by_row[[experiment_code]][[model_type]]
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
})) %>%
  mutate(
    model = factor(model, levels = OBSERVED_VS_PREDICTED_MODELS)
  )

write_csv(all_predictions, "outputs/observed_vs_predicted.csv")

## The row of each panel: "All Languages" first, then the experiments in the
## order label_language() gives them
all_predictions <- all_predictions %>%
  mutate(
    row_label = factor(
      ifelse(
        experiment == "all_experiments",
        "All Languages",
        as.character(label_language(experiment))
      ),
      levels = c("All Languages", levels(label_language(experiment)))
    )
  )

## Points are drawn for distances 0, 1 and 3. A few sentences (43 with three
## exchanges, 2 shuffled) come out at distance 2; they are in the models, but
## not drawn.
print(paste0("Item-condition cells at distance 2, not drawn: ",
             sum(all_predictions$dl_distance == 2 &
                   all_predictions$experiment != "all_experiments") /
               length(OBSERVED_VS_PREDICTED_MODELS)))

## The points of Figures 2 and 3 (distances 0, 1 and 3). Higher-surprisal
## sentences are drawn last, so they sit on top.
plot_data <- all_predictions %>%
  filter(dl_distance %in% c(0, 1, 3)) %>%
  mutate(dl_label = factor(dl_distance, levels = c(0, 1, 3))) %>%
  arrange(sum_surprisal)

## Pearson correlation between observed and predicted acceptability in each
## panel, over the points shown, for all five models
panel_correlations <- plot_data %>%
  group_by(experiment, row_label, model) %>%
  summarize(r = cor(observed, predicted),
            ## Root mean squared error: how far the points sit from the dashed
            ## line, in rating points (r ignores offsets from it)
            rmse = sqrt(mean((observed - predicted)^2)),
            n = n(), .groups = "drop") %>%
  mutate(label = paste0("r = ", sprintf("%.2f", r)))

print(panel_correlations %>%
        select(model, row_label, r) %>%
        pivot_wider(names_from = row_label, values_from = r),
      width = Inf)

write_csv(panel_correlations %>% select(experiment, model, r, rmse, n),
          "outputs/observed_vs_predicted_correlations.csv")

## Mean prediction error (observed minus predicted) in each panel, over the
## points shown, for all five models: separately at each DL distance (0, 1, 3),
## and over the sentences with errors together ("1 or 3"). A positive error
## means the model underestimates acceptability, a negative one that it
## overestimates it. The figure caption gives the range of these errors over
## the rows: SLOR and MORCELA underestimate the error-free sentences and
## overestimate the sentences with errors in every row.
panel_errors <- bind_rows(
  plot_data %>% mutate(dl_group = as.character(dl_distance)),
  plot_data %>% filter(dl_distance > 0) %>% mutate(dl_group = "1 or 3")
) %>%
  group_by(experiment, model, dl_group) %>%
  summarize(mean_error = mean(observed - predicted),
            n = n(), .groups = "drop")

## The range of the mean errors over the 10 rows, for each model
print(panel_errors %>%
        group_by(model, dl_group) %>%
        summarize(min_error = min(mean_error),
                  max_error = max(mean_error), .groups = "drop"),
      n = Inf)

write_csv(panel_errors, "outputs/observed_vs_predicted_errors.csv")

## ---- Figures 2 and 3: All Five Models ----
## Two views of how the five models' predictions compare with the ratings, both
## in the main text. Both use the predictions computed above (random effects
## excluded), one row for the pooled fit and one per experiment, and the models
## in the order and with the names of Figures 1 and 4.

## The DL colours of Figure 1, for the distances drawn here
DL_COLORS_SHOWN <- setNames(DL_COLORS[c("0", "1-2", "3")], c("0", "1", "3"))

all_models_data <- plot_data %>%
  mutate(
    model = factor(as.character(model), levels = MODEL_LEVELS, labels = MODEL_LABELS),
    ## "All Languages" on two lines, so the row labels take less width
    row_label = fct_relabel(row_label, ~ sub("All Languages", "All\nLanguages", .x)),
    ## Observed minus predicted: positive when a sentence is rated higher than
    ## the model predicts
    residual = observed - predicted,
    ## Larger points in the per-experiment rows, which have far fewer of them
    point_size = ifelse(experiment == "all_experiments", 0.35, 0.7)
  )

all_models_correlations <- panel_correlations %>%
  mutate(
    model = factor(as.character(model), levels = MODEL_LEVELS, labels = MODEL_LABELS),
    row_label = fct_relabel(row_label, ~ sub("All Languages", "All\nLanguages", .x))
  )

## Shared layout for the two figures: language labels on the left, axis lines on
## every panel, tick labels on the outer panels only, legend at the bottom
all_models_layout <- theme_print_8pt +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1),
        panel.spacing.x = unit(6, "pt"),
        panel.spacing.y = unit(4, "pt"),
        legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-4, 0, 0, 0),
        plot.margin = margin(2, 4, 2, 2))

## Residuals against sentence surprisal, with a smooth line per DL distance. A
## line above zero means the model underpredicts those sentences; below zero,
## it overpredicts them. Cells whose mean rating is exactly 7 (or exactly 1)
## trace the model's prediction curve upside down, so some points line up.
p_residuals <- ggplot(
  all_models_data,
  aes(x = sum_surprisal, y = residual, colour = dl_label)
) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.3) +
  geom_point(aes(size = point_size), alpha = 0.35, stroke = 0) +
  scale_size_identity() +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5), se = FALSE, linewidth = 0.6) +
  scale_colour_manual(
    name = "DL distance",
    values = DL_COLORS_SHOWN,
    guide = guide_legend(override.aes = list(size = 2, alpha = 1, linewidth = 0.6))
  ) +
  ## One y scale for all panels, wide enough for every residual (-5.0 to 6.0)
  scale_y_continuous(breaks = seq(-4, 6, 2)) +
  facet_grid(row_label ~ model, switch = "y", axes = "all", axis.labels = "margins",
             scales = "free_x") +
  labs(x = "Sentence surprisal (bits)", y = "Observed minus predicted acceptability") +
  all_models_layout +
  ## Pull the y title toward the panels, into the space left of the shorter
  ## (right-aligned) language labels
  theme(axis.title.y = element_text(margin = margin(r = -8, unit = "pt")))

ggsave(p_residuals, filename = "figures/residuals_by_surprisal.png",
       width = 7, height = 7.9, units = "in", dpi = 300)

## Observed against predicted for all five models. Colour shows the DL distance
## and, within it, sentence surprisal: light = low, dark = high, on one shared
## range of bits. Each DL distance gets its own colour scale (ggnewscale).
x_range_all <- range(all_models_data$predicted)
surprisal_limits <- range(all_models_data$sum_surprisal)

p_observed_vs_predicted_all <- ggplot(all_models_data, aes(x = predicted, y = observed)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "grey60",
              linewidth = 0.3) +
  ## DL 0, on a colour scale of its own
  ggnewscale::new_scale_colour() +
  geom_point(
    data = filter(all_models_data, dl_distance == 0),
    aes(colour = sum_surprisal, size = point_size),
    stroke = 0
  ) +
  scale_colour_gradient(
    name = "DL 0: surprisal (bits)",
    low = colorspace::lighten(DL_COLORS_SHOWN[["0"]], 0.75),
    high = colorspace::darken(DL_COLORS_SHOWN[["0"]], 0.45),
    limits = surprisal_limits,
    breaks = c(100, 250, 400),
    guide = guide_colourbar(barwidth = unit(55, "pt"), barheight = unit(5, "pt"),
                            title.position = "top", order = 1)
  ) +
  ## DL 1, on a colour scale of its own
  ggnewscale::new_scale_colour() +
  geom_point(
    data = filter(all_models_data, dl_distance == 1),
    aes(colour = sum_surprisal, size = point_size),
    stroke = 0
  ) +
  scale_colour_gradient(
    name = "DL 1: surprisal (bits)",
    low = colorspace::lighten(DL_COLORS_SHOWN[["1"]], 0.75),
    high = colorspace::darken(DL_COLORS_SHOWN[["1"]], 0.45),
    limits = surprisal_limits,
    breaks = c(100, 250, 400),
    guide = guide_colourbar(barwidth = unit(55, "pt"), barheight = unit(5, "pt"),
                            title.position = "top", order = 2)
  ) +
  ## DL 3, on a colour scale of its own
  ggnewscale::new_scale_colour() +
  geom_point(
    data = filter(all_models_data, dl_distance == 3),
    aes(colour = sum_surprisal, size = point_size),
    stroke = 0
  ) +
  scale_colour_gradient(
    name = "DL 3: surprisal (bits)",
    low = colorspace::lighten(DL_COLORS_SHOWN[["3"]], 0.75),
    high = colorspace::darken(DL_COLORS_SHOWN[["3"]], 0.45),
    limits = surprisal_limits,
    breaks = c(100, 250, 400),
    guide = guide_colourbar(barwidth = unit(55, "pt"), barheight = unit(5, "pt"),
                            title.position = "top", order = 3)
  ) +
  scale_size_identity() +
  geom_label(
    data = all_models_correlations,
    aes(x = x_range_all[1], y = 7, label = label),
    hjust = 0,
    vjust = 1,
    size = 2.2,
    fill = "white",
    alpha = 0.8,
    label.size = 0,
    label.padding = unit(0.6, "pt"),
    inherit.aes = FALSE
  ) +
  facet_grid(row_label ~ model, switch = "y", axes = "all", axis.labels = "margins") +
  coord_fixed(ratio = 1, xlim = x_range_all, ylim = c(1, 7)) +
  scale_x_continuous(breaks = seq(1, 7, 2)) +
  scale_y_continuous(breaks = seq(1, 7, 2)) +
  labs(x = "Predicted acceptability", y = "Observed mean acceptability") +
  all_models_layout +
  theme(legend.title = element_text(size = 7))

## 7.6 in tall: the panels have a fixed aspect ratio (coord_fixed), so a taller
## figure only adds blank space above and below them
ggsave(p_observed_vs_predicted_all, filename = "figures/observed_vs_predicted.png",
       width = 7, height = 7.6, units = "in", dpi = 300)

## ---- Figure 1: Qualitative Predictions (Model Curves over Surprisal) ----
## Human judgments and the predictions of the binary, gradient and communication
## models, drawn over sentence surprisal. Like Figures 2 and 3, it is drawn at
## its print size: PNAS full width (17.8 cm, 7 in), with 8-pt text.

#### ---- Helper Functions
## The binary, gradient and communication models are functions of surprisal and
## of the number of errors, so we can ask them what they predict for any
## combination of the two, whether or not a sentence like that was in the
## materials. That is what this grid is for. It runs in 1-bit steps over the
## given surprisal range, padded by a quarter bit at each end. predict() needs
## an item column even though the item effects are left out, so the grid names
## the model's first item.
get_prediction_grid <- function(
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
      item=factor(models_by_language[[lang]][[model_type]]$model$item[1], levels=levels(models_by_language[[lang]][[model_type]]$model$item))
    )
}

### ---- Baseline data for plotting (item means) ----
## Item x condition means as in crit_by_items above, then averaged over the
## conditions of an item that share a DL distance (one point per item and DL
## distance), and deliberately WITHOUT the filter(dl_distance <= 3): the scatter
## also shows the sentences at DL >= 4 (nearly all of the shuffled condition), in
## the grey "4+" bucket, as context. The models (and hence the fitted curves
## drawn over these points) are still fit on distances up to 3 only.
plot_data_by_item <- crit %>%
  mutate(
    item = as.factor(item),
    language = as.factor(language)
  ) %>%
  group_by(language, item, sum_surprisal, dl_distance, condition) %>%
  summarize(
    response_num = mean(response),
    .groups = "drop"
  ) %>%
  group_by(language, item, sum_surprisal, dl_distance) %>%
  summarize(
    response = mean(response_num),
    .groups = "drop"
  ) %>%
  mutate(
    dl_distance_bucket = bucket_dl(dl_distance),
    language           = label_language(language)
  )

### ---- Baseline scatterplot (item means) ----
p_baseline <- ggplot(
  plot_data_by_item,
  aes(x = sum_surprisal, y = response, color = dl_distance_bucket)
) +
  geom_point(alpha = 0.3, size = 0.6, stroke = 0) +
  facet_wrap(~language, ncol = 1, strip.position = "left") +
  labs(x="Sentence surprisal (bits)", y="Rating") +
  theme_print_8pt +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    ## Facet strips on the outside, with the language names horizontal. Only
    ## this first column shows them; the other columns blank their strips.
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    panel.spacing = unit(4, "pt")
  ) +
  ## The model columns keep this y scale; the human-judgment column replaces it
  ## with one that has tick labels
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
    linewidth = 0.5,
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

### ---- Model-based plots ----
## Each model's predictions for each experiment, on a grid over the
## experiment's surprisal range. Predictions leave out every random effect
## (random_effect_labels()), so that the controlled English curves do not carry
## the by-item slopes of whichever item the grid happens to name.
#### ---- Binary model ----
grammatical_plot_data <- bind_rows(lapply(names(models_by_language), function(lang) {
  grid <- get_prediction_grid(lang, "grammatical", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal))
  preds <- predict(
    models_by_language[[lang]]$grammatical,
    newdata = grid,
    type = "response",
    se.fit = TRUE,
    exclude = random_effect_labels(models_by_language[[lang]]$grammatical)
  )
  tibble(
    language = label_language(lang),
    sum_surprisal = grid$sum_surprisal,
    dl_distance_bucket = bucket_dl(grid$dl_distance),
    response = preds$fit,
    se = preds$se.fit
  )
}))

## The binary model predicts one rating for any number of errors, so its line
## for sentences with errors is drawn in one neutral dark grey, a colour no
## model or DL distance uses
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
    linewidth = 0.5
  ) +
  facet_wrap(~language, ncol = 1) +
  labs(
    x = "Sentence surprisal (bits)",
    y = NULL
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  ggtitle("Binary\nGrammaticality") +
  scale_color_manual(
    values = c(
      "0" = DL_COLORS[["0"]],
      "1-2" = "grey25",
      "3" = "grey25",
      "4+" = "grey25"
    )
  ) +
  ggplot2::guides(color = "none")

#### ---- Gradient model ----
distance_plot_data <- bind_rows(lapply(names(models_by_language), function(lang) {
  grid <- get_prediction_grid(lang, "distance", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal))
  preds <- predict(
    models_by_language[[lang]]$distance,
    newdata = grid,
    type = "response",
    se.fit = TRUE,
    exclude = random_effect_labels(models_by_language[[lang]]$distance)
  )
  tibble(
    language = label_language(lang),
    sum_surprisal = grid$sum_surprisal,
    dl_distance_bucket = bucket_dl(grid$dl_distance),
    response = preds$fit,
    se = preds$se.fit
  )
}))

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
    linewidth = 0.5
  ) +
  facet_wrap(~language, ncol = 1) +
  labs(
    x = "Sentence surprisal (bits)",
    y = NULL
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  ggtitle("Gradient\nGrammaticality") +
  ggplot2::guides(color = "none")

#### ---- Communication model ----
weighted_distance_plot_data <- bind_rows(lapply(names(models_by_language), function(lang) {
  grid <- get_prediction_grid(lang, "weighted_distance", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal))
  preds <- predict(
    models_by_language[[lang]]$weighted_distance,
    newdata = grid,
    type = "response",
    se.fit = TRUE,
    exclude = random_effect_labels(models_by_language[[lang]]$weighted_distance)
  )
  tibble(
    language = label_language(lang),
    sum_surprisal = grid$sum_surprisal,
    dl_distance_bucket = bucket_dl(grid$dl_distance),
    response = preds$fit,
    se = preds$se.fit
  )
}))

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
    linewidth = 0.5
  ) +
  facet_wrap(~language, ncol = 1) +
  labs(
    x = "Sentence surprisal (bits)",
    y = NULL
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  ggtitle("Communication") +
  ggplot2::guides(color = "none")

### ---- Combine plots into the qualitative-predictions figure ----
## The models in the order and with the names of Figures 2-4. The four columns
## share one legend and one x-axis title. The column titles are at the base
## text size, with "Grammaticality" on a second line so they fit over their
## columns. The figure as a whole (its margins and the shared legend) also takes
## the print theme, rather than the script's 18-pt default.
qualitative_predictions_figure <- (p_gam | p_grammatical | p_distance | p_weighted_distance) +
  plot_layout(guides = "collect", axis_titles = "collect") +
  plot_annotation(theme = theme_print_8pt) &
  theme(legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"),
        plot.title = element_text(hjust = 0.5, size = rel(1))
        )

ggsave(
  filename = "figures/qualitative_predictions.png",
  plot     = qualitative_predictions_figure,
  width    = 7,
  height   = 7.5,
  units    = "in",
  dpi      = 300
)

print("Done.")
