## =============================================================================
## SI Section 5: SLOR and MORCELA with other language models and unigram models
##
## In the main text, SLOR and MORCELA take a sentence's log probability from the
## Goldfish 1000mb model of its language, and its words' unigram probabilities
## from the same model (each word scored alone at the start of a sentence: the
## "model-initial" unigram). Would either account do better with another
## language model or another unigram model? This section asks that for the
## English experiment (en), over 18 settings, 6 language models x 3 unigram
## models:
##   - language models: Goldfish 1000mb (the main text's), GPT-2, GPT-2 XL, and
##     interpolated Kneser-Ney 2-, 3- and 4-gram models trained on about 20M
##     words of WikiText-103 (the kind of n-gram model Lau, Clark & Lappin 2017
##     used);
##   - unigram models: model-initial (the main text's), wordfreq (word
##     frequencies from the wordfreq package), and WikiText-103 (word
##     frequencies in the n-gram models' training text).
## The scores are in ../data/slor_morcela_scores_english_settings.csv; see
## ../data/DATA_DICTIONARY.md for how they were computed.
##
## Under every setting, SLOR and MORCELA are fit exactly as
## ../analysis/main_analysis.R fits them for en -- the same item x condition
## means (DL distance <= 3), the same increasing smooth of the score with an
## item random intercept, MORCELA's two weights learned from those means -- and
## compared, by the same marginal log-likelihood, with the main text's Binary,
## Gradient and Communication models for en, which are refit here identically.
## Goldfish + model-initial is the main text's setting, so its two fits must
## reproduce the main text's SLOR and MORCELA results for en; the script stops
## if they do not.
##
## Takes the libraries, the plotting theme, the helper functions
## (marginal_logLik(), morcela_weights(), morcela_score()) and the five models'
## colours, shapes and legend labels from setup.R, but builds the English item
## means itself, exactly as ../analysis/main_analysis.R does, in place of
## setup.R's crit and crit_by_items (which hold all nine experiments). Each fit
## takes under a second, so nothing is cached.
##
## Produces figures/s5_slor_morcela_settings.png/.pdf and writes every
## log-likelihood and MORCELA weight to
## outputs/s5_slor_morcela_settings_likelihoods.csv.
## =============================================================================

## Print warnings as they happen, so they show up next to the step that raised
## them in the run logs
options(warn = 1)

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

## Libraries, plotting theme, helper functions and the five models' colours,
## shapes and legend labels all live in setup.R, which also creates the output
## folders. Its marginal_logLik() is the main analysis's marginal
## log-likelihood: scam/gam's logLik() is conditional on the fitted item random
## effects, and marginal_logLik() integrates them out, as lmer() does, so that
## model comparison does not credit a model for them.
source("setup.R")

## ---- Read and Prepare Data ----
## As in ../analysis/main_analysis.R: drop the participants who answered fewer
## than 80% of the comprehension questions correctly, keep the rating trials, and
## drop the nonsense conditions. Only the English experiment is used here.
dat <- read_csv("../data/full_data.csv", show_col_types = FALSE)

## Exclude participants with a low comprehension score
exclusions <- dat %>%
  filter(grepl("comprehension", experiment)) %>%
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

## The critical ratings of the English experiment, converted to numbers
crit <- dat %>%
  filter(
    !participant_id %in% exclusions,
    experiment == "exp_b",
    !grepl("non_sense", condition),
    language == "en"
  ) %>%
  mutate(
    response_num = as.numeric(substr(response, 1, 1))
  )

## ---- Surprisal of the Intended Sentence ----
## The Communication model's Information(s) is the main text's: the surprisal of
## the item's original sentence under the Goldfish model, from
## ../data/slor_morcela_scores.csv, the same in all of an item's conditions. The
## item ids there lack the en_ prefix of full_data.csv, so it is added.
main_scores <- read_csv(
  "../data/slor_morcela_scores.csv",
  col_types = cols(item = col_character())
) %>%
  filter(language == "en") %>%
  mutate(item = paste("en", item, sep = "_"))

original_surprisal <- main_scores %>%
  filter(unigram_model == "model_initial", condition == "original") %>%
  transmute(item, sum_surprisal = -log_p_sentence)

n_crit_rows <- nrow(crit)

crit <- crit %>%
  left_join(original_surprisal, by = "item")

## Every rating found exactly one surprisal
stopifnot(
  nrow(crit) == n_crit_rows,
  !anyNA(crit$sum_surprisal)
)

## ---- Item Means ----
## One row per sentence (item x condition) with its mean rating, keeping DL
## distances up to 3, as in ../analysis/main_analysis.R: 219 of the 288
## sentences.
crit_by_items <- crit %>%
  mutate(
    weighted_dl_distance_surprisal = dl_distance / sum_surprisal,
    grammatical = factor(
      ifelse(dl_distance == 0, "grammatical", "ungrammatical"),
      levels=c("ungrammatical", "grammatical")),
    item = as.factor(item)
  ) %>%
  group_by(item, n_words, sum_surprisal, dl_distance,
           weighted_dl_distance_surprisal, grammatical, condition) %>%
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

## Each sentence is one row
stopifnot(!anyDuplicated(crit_by_items[c("item", "condition")]))

nrow(crit_by_items)

## ---- Scores Under Every Setting ----
## log_p_sentence and log_p_unigram (bits) of each English sentence under each of
## the 18 settings, computed like ../data/slor_morcela_scores.csv: on the
## sentence the participant saw (perceived_sentence), with words split at
## spaces.
LANGUAGE_MODELS <- c("goldfish", "gpt2", "gpt2-xl", "kn2", "kn3", "kn4")
UNIGRAM_MODELS  <- c("model_initial", "wordfreq", "corpus")

settings_scores <- read_csv(
  "../data/slor_morcela_scores_english_settings.csv",
  col_types = cols(item = col_character())
) %>%
  mutate(item = paste("en", item, sep = "_"))

## Every setting scores each of the 288 sentences (72 items x 4 conditions) once
stopifnot(
  setequal(settings_scores$language_model, LANGUAGE_MODELS),
  setequal(settings_scores$unigram_model, UNIGRAM_MODELS),
  nrow(settings_scores) == length(LANGUAGE_MODELS) * length(UNIGRAM_MODELS) * 288,
  !anyDuplicated(settings_scores[c("language_model", "unigram_model", "item", "condition")]),
  all(is.finite(settings_scores$log_p_sentence)),
  all(is.finite(settings_scores$log_p_unigram))
)

## These scores were computed by an earlier script than
## ../data/get_slor_morcela_scores.py. Its Goldfish scores must equal the main
## text's, for both unigram models the two files share. (They agree to about
## 1e-13 bits.)
goldfish_check <- settings_scores %>%
  filter(
    language_model == "goldfish",
    unigram_model %in% c("model_initial", "wordfreq")
  ) %>%
  inner_join(
    main_scores %>% select(item, condition, unigram_model, log_p_sentence, log_p_unigram),
    by = c("item", "condition", "unigram_model"),
    suffix = c("", "_main_text")
  )

stopifnot(
  nrow(goldfish_check) == 2 * 288,
  max(abs(goldfish_check$log_p_sentence - goldfish_check$log_p_sentence_main_text)) < 1e-9,
  max(abs(goldfish_check$log_p_unigram - goldfish_check$log_p_unigram_main_text)) < 1e-9
)

## ---- The Main Text's Models ----
## Binary, Gradient and Communication, fit as ../analysis/main_analysis.R fits
## them for en, and under its names (grammatical, distance and
## weighted_distance). They do not depend on the SLOR and MORCELA settings.
grammatical <- scam(
  response_num ~ grammatical_numeric +
    s(item, bs="re"),
  data=crit_by_items
)

distance <- scam(
  response_num ~ s(dl_distance, bs="mpd") +
    s(item, bs="re"),
  data=crit_by_items
)

weighted_distance <- scam(
  response_num ~ s(weighted_dl_distance_surprisal, bs="mpd") + s(sum_surprisal, bs="mpd") +
    s(item, bs="re"),
  data=crit_by_items
)

main_model_likelihoods <- tibble(
  model = c("Binary", "Gradient", "Communication"),
  language_model = NA_character_,
  unigram_model = NA_character_,
  in_main_text = TRUE,
  log_likelihood = c(marginal_logLik(grammatical), marginal_logLik(distance),
                     marginal_logLik(weighted_distance)),
  converged = c(isTRUE(grammatical$conv), isTRUE(distance$conv),
                isTRUE(weighted_distance$conv)),
  morcela_beta = NA_real_,
  morcela_gamma = NA_real_
)

## ---- SLOR and MORCELA Under Every Setting ----
## As in ../analysis/main_analysis.R.
##
##     SLOR    = ( log p_model(sentence) -        log p_unigram(sentence)         ) / length
##     MORCELA = ( log p_model(sentence) - beta * log p_unigram(sentence) + gamma ) / length
##
## with length in words. MORCELA's two weights are learned from the ratings, by
## the recipe of its paper: regress the item means on log p_model(sentence),
## log p_unigram(sentence) and 1, each divided by length, and read the weights
## off the coefficients (setup.R's morcela_weights()). They are learned afresh
## under every setting, on the same item means the models are then fit to.
settings <- expand_grid(
  language_model = LANGUAGE_MODELS,
  unigram_model = UNIGRAM_MODELS
)

## Fit SLOR and MORCELA under each setting in turn, one row per fit
settings_likelihoods <- tibble()
for (setting_index in seq_len(nrow(settings))) {
  language_model_iter <- settings$language_model[setting_index]
  unigram_model_iter <- settings$unigram_model[setting_index]

  print(
    paste0("Fitting SLOR and MORCELA with ", language_model_iter, " + ", unigram_model_iter)
  )

  ## The item means with this setting's log probabilities joined on, and SLOR
  ## computed from them. The join turns `item` into a character vector, so it
  ## is made a factor again, with its original levels.
  crit_by_items_setting_specific <- crit_by_items %>%
    left_join(
      settings_scores %>%
        filter(
          language_model == language_model_iter,
          unigram_model == unigram_model_iter
        ) %>%
        select(item, condition, log_p_sentence, log_p_unigram),
      by = c("item", "condition")
    ) %>%
    mutate(
      item = factor(item, levels = levels(crit_by_items$item)),
      slor = (log_p_sentence - log_p_unigram) / n_words
    )

  ## Every sentence found exactly one score
  stopifnot(
    nrow(crit_by_items_setting_specific) == nrow(crit_by_items),
    !anyNA(crit_by_items_setting_specific$item),
    !anyNA(crit_by_items_setting_specific$log_p_sentence),
    !anyNA(crit_by_items_setting_specific$log_p_unigram)
  )

  ## MORCELA, with its two weights learned on this setting's item means
  setting_morcela_weights <- morcela_weights(crit_by_items_setting_specific)
  crit_by_items_setting_specific$morcela <- morcela_score(
    crit_by_items_setting_specific, setting_morcela_weights
  )

  ## Fit SLOR and MORCELA as ../analysis/main_analysis.R fits them for en: an
  ## increasing smooth of the score, since a higher score means the sentence is
  ## more probable than its words alone would predict, plus an item random
  ## intercept
  slor <- scam(
    response_num ~ s(slor, bs="mpi") +
      s(item, bs="re"),
    data=crit_by_items_setting_specific
  )

  morcela <- scam(
    response_num ~ s(morcela, bs="mpi") +
      s(item, bs="re"),
    data=crit_by_items_setting_specific
  )

  settings_likelihoods <- bind_rows(
    settings_likelihoods,
    tibble(
      model = c("SLOR", "MORCELA"),
      language_model = language_model_iter,
      unigram_model = unigram_model_iter,
      in_main_text = language_model_iter == "goldfish" & unigram_model_iter == "model_initial",
      log_likelihood = c(marginal_logLik(slor), marginal_logLik(morcela)),
      converged = c(isTRUE(slor$conv), isTRUE(morcela$conv)),
      morcela_beta = c(NA, setting_morcela_weights$beta),
      morcela_gamma = c(NA, setting_morcela_weights$gamma)
    )
  )
}

## ---- Log-Likelihoods ----
## Every log-likelihood also as its difference from the best model's, as in the
## main text's model-comparison figure
s5_likelihoods <- bind_rows(main_model_likelihoods, settings_likelihoods) %>%
  mutate(
    experiment = "en",
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood
  ) %>%
  select(experiment, model, language_model, unigram_model, in_main_text,
         log_likelihood, max_likelihood, likelihood_diff, converged,
         morcela_beta, morcela_gamma)

## Evaluate Convergence: the fits that did not converge (none, when this was
## written)
s5_likelihoods %>%
  filter(!converged)

## All 39 fits, best first
s5_likelihoods %>%
  arrange(desc(log_likelihood)) %>%
  print(n = Inf, width = Inf)

## Over the 18 settings: in how many SLOR and MORCELA fit better than each of
## the three main-text models, and their best n-gram (Kneser-Ney) and worst
## neural (Goldfish, GPT-2, GPT-2 XL) log-likelihoods
s5_likelihoods %>%
  filter(model %in% c("SLOR", "MORCELA")) %>%
  group_by(model) %>%
  summarize(
    n_settings = n(),
    n_better_than_binary = sum(log_likelihood > s5_likelihoods$log_likelihood[s5_likelihoods$model == "Binary"]),
    n_better_than_gradient = sum(log_likelihood > s5_likelihoods$log_likelihood[s5_likelihoods$model == "Gradient"]),
    n_better_than_communication = sum(log_likelihood > s5_likelihoods$log_likelihood[s5_likelihoods$model == "Communication"]),
    best_ngram_log_likelihood = max(log_likelihood[str_starts(language_model, "kn")]),
    worst_neural_log_likelihood = min(log_likelihood[!str_starts(language_model, "kn")])
  ) %>%
  print(width = Inf)

### ---- Check: the main setting reproduces the main text ----
## The three main-text models, and SLOR and MORCELA under Goldfish +
## model-initial, must give the log-likelihoods ../analysis/main_analysis.R
## reports for en, and MORCELA the weights it reports for en. On the machine
## that ran the main analysis, the three main-text models agree exactly and
## SLOR and MORCELA to within 1e-10, since their Goldfish scores here differ
## from the main text's by rounding (about 1e-13 bits). The tolerances allow
## for another machine's arithmetic.
main_text_likelihoods <- read_csv(
  "../analysis/outputs/likelihoods_by_language.csv",
  show_col_types = FALSE
) %>%
  filter(experiment == "en") %>%
  select(model, main_text_log_likelihood = log_likelihood)

main_text_morcela_weights <- read_csv(
  "../analysis/outputs/morcela_weights.csv",
  show_col_types = FALSE
) %>%
  filter(experiment == "en")

main_text_check <- s5_likelihoods %>%
  filter(in_main_text) %>%
  inner_join(main_text_likelihoods, by = "model") %>%
  mutate(difference = log_likelihood - main_text_log_likelihood)

main_text_check %>%
  select(model, log_likelihood, main_text_log_likelihood, difference)

main_text_morcela <- main_text_check %>%
  filter(model == "MORCELA")

stopifnot(
  nrow(main_text_check) == 5,
  all(abs(main_text_check$difference) < 1e-3),
  abs(main_text_morcela$morcela_beta - main_text_morcela_weights$beta) < 1e-6,
  abs(main_text_morcela$morcela_gamma - main_text_morcela_weights$gamma) < 1e-6
)

## One row per fit: first the three main-text models, which have no language
## model or unigram model of their own (NA; Communication's surprisal is the
## main text's, from Goldfish), then SLOR and MORCELA under each setting.
## in_main_text marks the five fits the main text reports for en.
## morcela_beta and morcela_gamma are MORCELA's two weights (NA for the other
## models; SLOR is beta = 1, gamma = 0).
write_csv(s5_likelihoods, "outputs/s5_slor_morcela_settings_likelihoods.csv")

## ---- Figure: Log-Likelihoods Under Every Setting ----
## Drawn like the main text's model-comparison figure, with its colours, markers
## and legend labels (setup.R's MODEL_LEVELS, MODEL_LABELS, MODEL_COLORS and
## MODEL_SHAPES): each log-likelihood as its difference from the best model's,
## one row per setting, grouped by language model, with SLOR and MORCELA as
## points. The three main-text models do not depend on the setting, so they are
## vertical lines.
LANGUAGE_MODEL_LABELS <- c(
  "goldfish" = "Goldfish",
  "gpt2"     = "GPT-2",
  "gpt2-xl"  = "GPT-2 XL",
  "kn2"      = "Kneser-Ney\n2-gram",
  "kn3"      = "Kneser-Ney\n3-gram",
  "kn4"      = "Kneser-Ney\n4-gram"
)
UNIGRAM_MODEL_LABELS <- c(
  "model_initial" = "Model-initial",
  "wordfreq"      = "wordfreq",
  "corpus"        = "WikiText-103"
)
MAIN_SETTING_LABEL <- "Model-initial\n(main text)"

## The rows of each language model run, top to bottom, in the order of
## UNIGRAM_MODEL_LABELS, and the main text's setting is labelled as such
settings_plot_data <- s5_likelihoods %>%
  filter(model %in% c("SLOR", "MORCELA")) %>%
  mutate(
    model = factor(model, levels = MODEL_LEVELS),
    language_model = factor(
      LANGUAGE_MODEL_LABELS[language_model],
      levels = LANGUAGE_MODEL_LABELS
    ),
    unigram_model = factor(
      ifelse(in_main_text, MAIN_SETTING_LABEL, UNIGRAM_MODEL_LABELS[unigram_model]),
      levels = rev(c(MAIN_SETTING_LABEL, UNIGRAM_MODEL_LABELS))
    )
  )

## Only the model and its difference: without a language_model column, the
## lines are drawn in every language model's panel
main_models_plot_data <- s5_likelihoods %>%
  filter(model %in% c("Binary", "Gradient", "Communication")) %>%
  transmute(
    model = factor(model, levels = MODEL_LEVELS),
    likelihood_diff
  )

p_settings_likelihoods <- ggplot(settings_plot_data) +
  geom_vline(
    data = main_models_plot_data,
    aes(
      xintercept=likelihood_diff,
      color=model
    ),
    linewidth=1
  ) +
  geom_point(
    aes(
      y=unigram_model,
      x=likelihood_diff,
      color=model,
      fill=model,
      shape=model
    ),
    size=3
  ) +
  facet_grid(
    language_model ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_color_manual(
    name="",
    values=MODEL_COLORS,
    breaks=MODEL_LEVELS,
    labels=MODEL_LABELS
  ) +
  scale_fill_manual(values=MODEL_COLORS) +
  scale_shape_manual(values=MODEL_SHAPES) +
  scale_x_continuous(
    name = "Log-Likelihood Difference"
  ) +
  ## As in the main figure, the three main-text models go on the first row of
  ## the legend and SLOR and MORCELA on the second. Their keys show what is
  ## drawn: a line for each main-text model, a marker for SLOR and MORCELA.
  guides(
    color = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(
        linetype = c("solid", "solid", "solid", "blank", "blank"),
        shape = c(NA, NA, NA, MODEL_SHAPES[["SLOR"]], MODEL_SHAPES[["MORCELA"]]),
        fill = c(NA, NA, NA, MODEL_COLORS[["SLOR"]], MODEL_COLORS[["MORCELA"]])
      )
    ),
    shape = "none",
    fill = "none"
  ) +
  ## The language models' names go to the left of the unigram models' names.
  ## The legend is centred on the whole figure rather than on the panels, which
  ## those names leave too narrow for it.
  theme(
    legend.position="top",
    legend.location = "plot",
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1)
  ) +
  scale_y_discrete(name="")

ggsave("figures/s5_slor_morcela_settings.png",
       plot = p_settings_likelihoods, width = 8, height = 10)

## Save as pdf, with the fonts embedded (cairo_pdf)
ggsave("figures/s5_slor_morcela_settings.pdf",
       plot = p_settings_likelihoods, width = 8, height = 10,
       device = cairo_pdf)
