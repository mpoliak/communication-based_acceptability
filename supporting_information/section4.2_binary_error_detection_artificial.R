## =============================================================================
## SI Section 4.2: Binary error detection and acceptability on the artificial materials
##
## en_artificial_binary_error_detection. The 27 constructed items of section 4.1, each
## in a short, medium and long version, each original, with 1 swap or with 3
## swaps: 9 conditions, every participant sees every item once. On each trial
## the participant first answers "Are there any errors in this sentence?"
## (Yes/No) and then, on the next screen, rates how acceptable the sentence is
## (1-7). The error question is correct when the answer matches
## literal_response: "No" for originals, "Yes" for every swap.
##
## The same questions as the binary error detection part of section4_binary_error_detection.R:
##   (1) How well do participants detect the errors?
##   (2) Do ratings depend on whether the error was judged correctly?
##   (3) Among originals judged error-free, are longer sentences rated lower?
##   (4) Which of the paper's three accounts fits best -- on all trials, on the
##       correctly judged trials only, and adding back originals wrongly flagged?
##       Two length baselines (gradient + length, length only) test whether a
##       plain length term does the communication model's work. Every model
##       with an information term is fitted twice: with sentence length and
##       with sentence surprisal (the main analysis's information measure).
## and, because length varies WITHIN item here, the within-item test of section 4:
##   (5) Does the same swap cost less in a longer version of the same sentence?
## The comparison with a comprehension-question version of the experiment is
## not repeated: no such version was run on these materials.
##
## Mixed models use the fullest random-effect structure that converges: each
## is tried from the maximal structure down, stopping at the first without a
## convergence failure (singular fits are ignored, as throughout).
##
## Figures: figures/artificial_binary_error_detection_*.png
##
## Section 4.3 re-runs this whole script on the sentences validated in the
## exchange identification task, by setting keep_sentences and figure_prefix and then sourcing
## this file. Run on its own, it analyses every sentence.
## =============================================================================

library(tidyverse)
library(lme4)
library(lmerTest)
library(scam)
library(stringdist)
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
dir.create("figures", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)

## ---- Colours and optimizer settings ----
condition_colors <- c("Original" = "#EB7AA3", "1 swap" = "#1E88E5", "3 swaps" = "#FFC107")
model_colors <- c("Binary grammaticality" = "#EB7AA3",
                  "Gradient grammaticality" = "#1E88E5",
                  "Communication-based" = "#004D40")
lmer_ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
glmer_ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

## ---- Which sentences to analyse ----
## Every sentence, unless section 4.3 has set keep_sentences (the item x length
## combinations to keep) and figure_prefix before sourcing this file.
if (!exists("keep_sentences")) keep_sentences <- NULL
if (!exists("figure_prefix")) figure_prefix <- "figures/artificial_binary_error_detection_"

## ---- Materials: one row per sentence version ----
## The error question is the anchor row and carries the sentence; the rating is
## its coordinated follow-up and has an empty sentence.
materials <- read_csv("../data/en_artificial_binary_error_detection/materials_en_artificial_binary_error_detection.csv",
                      show_col_types = FALSE) %>%
  filter(experiment == "exp_b_comprehension") %>%
  transmute(item, cond = condition, sentence,
            length_level = factor(str_extract(condition, "[^_]+$"),
                                  levels = c("short", "medium", "long")),
            condition_label = factor(recode(str_extract(condition, "^[^_]+"),
                                            "original" = "Original",
                                            "swap-1" = "1 swap",
                                            "swap-3" = "3 swaps"),
                                     levels = names(condition_colors)))

## Damerau-Levenshtein distance from each version to its original, computed as
## in data/preprocess_stimuli.R: lowercase, drop punctuation, write each
## distinct word of the item (at that length) as one letter, and take
## stringdist(method = "dl") between the letter strings.
sentence_letters <- materials %>%
  separate_rows(sentence, sep = " ") %>%
  mutate(word = str_remove_all(str_to_lower(sentence), "[^a-z0-9]")) %>%
  group_by(item, length_level) %>%
  mutate(letter = c(letters, LETTERS)[match(word, unique(word))]) %>%
  group_by(item, length_level, cond, condition_label) %>%
  summarize(letter_string = paste(letter, collapse = ""), n_words = n(), .groups = "drop")

## Sentence surprisal, the main analysis's information measure: the summed
## Goldfish surprisal (bits) of the original at each length, from
## data/slor_morcela_scores.csv (scored on the lowercased, unpunctuated
## sentence). Like length, it differs between an item's short, medium and long
## versions but is shared by the original and swapped versions of each.
sentence_surprisal <- read_csv("../data/slor_morcela_scores.csv",
                               col_types = cols(item = col_character()),
                               show_col_types = FALSE) %>%
  filter(language == "en_artificial", unigram_model == "model_initial",
         str_starts(condition, "original")) %>%
  transmute(item = as.integer(item),
            length_level = factor(str_extract(condition, "[^_]+$"),
                                  levels = c("short", "medium", "long")),
            sum_surprisal = -log_p_sentence)
stopifnot(nrow(sentence_surprisal) == 27 * 3,
          !anyDuplicated(sentence_surprisal[c("item", "length_level")]))

sentences <- sentence_letters %>%
  left_join(sentence_letters %>%
              filter(condition_label == "Original") %>%
              select(item, length_level, original_string = letter_string),
            by = c("item", "length_level")) %>%
  left_join(sentence_surprisal, by = c("item", "length_level")) %>%
  mutate(dl_distance = stringdist(letter_string, original_string, method = "dl"),
         weighted_dl_distance = dl_distance / n_words,
         weighted_dl_distance_surprisal = dl_distance / sum_surprisal,
         grammatical = factor(ifelse(dl_distance == 0, "grammatical", "ungrammatical"),
                              levels = c("ungrammatical", "grammatical"))) %>%
  select(item, cond, length_level, condition_label, n_words, dl_distance,
         weighted_dl_distance, sum_surprisal, weighted_dl_distance_surprisal, grammatical)
stopifnot(!anyNA(sentences))

cat("DL distance by condition (number of sentence versions):\n")
sentences %>% count(condition_label, dl_distance) %>% print()

## ---- Responses: one row per participant x sentence ----
responses <- read_csv("../data/en_artificial_binary_error_detection/en_artificial_binary_error_detection_data.csv",
                      show_col_types = FALSE)

## Participants who also took the exchange identification + acceptability
## ratings study saw the same sentences twice; they are excluded from both.
took_both <- intersect(
  str_trim(responses$participant_id),
  str_trim(read_csv("../data/en_artificial_exchange_identification_acceptability/en_artificial_exchange_identification_acceptability_data.csv",
                    show_col_types = FALSE)$participant_id))
cat("excluded for also taking exchange identification + acceptability ratings:",
    length(took_both), "participant(s)\n")
responses <- responses %>% filter(!str_trim(participant_id) %in% took_both)

judgments <- responses %>%
  filter(experiment == "exp_b_comprehension") %>%
  transmute(participant_id, item, cond, judgment_order = presentation_order,
            said_error = response == "Yes",
            detected = response == literal_response)

ratings <- responses %>%
  filter(experiment == "exp_b") %>%
  transmute(participant_id, item, cond, rating_order = presentation_order,
            acceptability = as.numeric(substr(response, 1, 1)))

cat("\n---- export integrity ----\n")
cat("participants:", n_distinct(responses$participant_id), "\n")
responses %>% count(participant_id, name = "rows") %>% count(rows, name = "participants") %>% print()
cat("(54 rows = 27 sentences x 2 questions; Surveyor writes no row for a skipped question)\n")
cat("judgments without a rating:",
    nrow(anti_join(judgments, ratings, by = c("participant_id", "item", "cond"))),
    "| ratings without a judgment:",
    nrow(anti_join(ratings, judgments, by = c("participant_id", "item", "cond"))), "\n")

trials <- judgments %>%
  inner_join(ratings, by = c("participant_id", "item", "cond")) %>%
  left_join(sentences, by = c("item", "cond")) %>%
  mutate(detection_label = factor(ifelse(detected, "Error judged correctly",
                                         "Error judged incorrectly"),
                                  levels = c("Error judged correctly",
                                             "Error judged incorrectly")))

if (!is.null(keep_sentences)) {
  trials <- trials %>% semi_join(keep_sentences, by = c("item", "length_level"))
  sentences <- sentences %>% semi_join(keep_sentences, by = c("item", "length_level"))
  cat("\n*** RESTRICTED to", nrow(keep_sentences),
      "validated item x length combinations (all three versions of each) ***\n")
}

## The rating should always be the screen right after its error question.
cat("rating shown right after its error question:",
    sum(trials$rating_order == trials$judgment_order + 1), "of", nrow(trials), "trials\n")

cells <- trials %>% count(item, cond)
cat("item x condition cells seen:", nrow(cells), "of", nrow(sentences),
    "| responses per cell: median", median(cells$n),
    "range", min(cells$n), "-", max(cells$n), "\n")

## Predictors z-scored for the mixed models; a linear rescaling leaves the
## log-likelihoods unchanged, so the model comparisons are unaffected.
model_trials <- trials %>%
  mutate(participant_id = factor(participant_id),
         item = factor(item),
         sentence_id = paste(item, length_level, sep = "_"),
         dl_distance_z = as.numeric(scale(dl_distance)),
         weighted_dl_distance_z = as.numeric(scale(weighted_dl_distance)),
         n_words_z = as.numeric(scale(n_words)),
         weighted_dl_distance_surprisal_z = as.numeric(scale(weighted_dl_distance_surprisal)),
         sum_surprisal_z = as.numeric(scale(sum_surprisal)),
         detected_c = ifelse(detected, 0.5, -0.5))

## ---- (1) Error detection ----
cat("\n===== (1) Error detection =====\n")
trials %>%
  group_by(condition_label) %>%
  summarize(trials = n(), accuracy = mean(detected), said_error = mean(said_error),
            mean_rating = mean(acceptability), .groups = "drop") %>%
  print()
trials %>%
  group_by(condition_label, length_level) %>%
  summarize(trials = n(), accuracy = mean(detected), mean_rating = mean(acceptability),
            .groups = "drop") %>%
  print(n = Inf)
trials %>%
  summarize(overall_accuracy = mean(detected),
            false_alarm_rate = mean(said_error[condition_label == "Original"]),
            hit_rate = mean(said_error[condition_label != "Original"])) %>%
  print()

## Does noticing an error depend on how much was displaced and how long the
## sentence is? Length varies within item, so it can take a by-item slope.
detection_formulas <- c(
  "detected ~ dl_distance_z * n_words_z + (dl_distance_z * n_words_z | participant_id) + (dl_distance_z * n_words_z | item)",
  "detected ~ dl_distance_z * n_words_z + (dl_distance_z + n_words_z | participant_id) + (dl_distance_z + n_words_z | item)",
  "detected ~ dl_distance_z * n_words_z + (dl_distance_z + n_words_z || participant_id) + (dl_distance_z + n_words_z || item)",
  "detected ~ dl_distance_z * n_words_z + (1 | participant_id) + (1 | item)")
for (f in detection_formulas) {
  m_detection <- suppressWarnings(glmer(as.formula(f), family = binomial,
                                        data = model_trials, control = glmer_ctrl))
  problems <- grep("singular", m_detection@optinfo$conv$lme4$messages,
                   ignore.case = TRUE, invert = TRUE, value = TRUE)
  if (length(problems) == 0) break
}
cat("\n---- correct judgment ~ DL distance x length ----\n")
cat("random effects:", f, "\n")
print(round(summary(m_detection)$coefficients, 4))

## ---- Figure 1: probability of a correct judgment, by length ----
## Points are sentences (one item at one length), jittered sideways so items of
## the same length stay visible; lines are linear fits to the trial-level 0/1
## data.
p_detection <- ggplot(mapping = aes(x = n_words, color = condition_label,
                                    fill = condition_label)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey70") +
  geom_point(data = trials %>%
               group_by(condition_label, item, length_level, n_words) %>%
               summarize(value = mean(detected), .groups = "drop"),
             aes(y = value), alpha = 0.4,
             position = position_jitter(width = 0.3, height = 0.02, seed = 2026)) +
  geom_smooth(data = trials, aes(y = as.integer(detected)), method = "lm", formula = y ~ x) +
  scale_color_manual(name = "Condition", values = condition_colors,
                     aesthetics = c("colour", "fill")) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Probability of correct judgment", labels = scales::percent,
                     breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1))
p_detection
ggsave(paste0(figure_prefix, "accuracy_by_length.png"), p_detection,
       width = 8, height = 5.5)

## ---- (2) Do ratings depend on whether the error was judged correctly? ----
## Acceptability from the communication predictors (DL distance per word, plus
## length) crossed with whether the error judgment was correct. Correctness is
## sum-coded so the predictor effects average over it.
cat("\n===== (2) Ratings and binary error detection =====\n")
awareness_formulas <- c(
  "acceptability ~ (weighted_dl_distance_z + n_words_z) * detected_c + (weighted_dl_distance_z + n_words_z + detected_c | participant_id) + (weighted_dl_distance_z + n_words_z + detected_c | item)",
  "acceptability ~ (weighted_dl_distance_z + n_words_z) * detected_c + (weighted_dl_distance_z + n_words_z + detected_c | participant_id) + (weighted_dl_distance_z | item)",
  "acceptability ~ (weighted_dl_distance_z + n_words_z) * detected_c + (weighted_dl_distance_z + n_words_z + detected_c || participant_id) + (1 | item)")
for (f in awareness_formulas) {
  m_awareness <- suppressWarnings(lmer(as.formula(f), data = model_trials, control = lmer_ctrl))
  problems <- grep("singular", m_awareness@optinfo$conv$lme4$messages,
                   ignore.case = TRUE, invert = TRUE, value = TRUE)
  if (length(problems) == 0) break
}
cat("---- acceptability ~ (DL/length + length) x judged correctly ----\n")
cat("random effects:", f, "\n")
print(round(summary(m_awareness)$coefficients, 4))

## ---- Figure 2: acceptability by length, split by the error judgment ----
## Points are sentence means (one item at one length); lines are linear fits.
p_by_detection <- trials %>%
  group_by(detection_label, condition_label, item, length_level, n_words) %>%
  summarize(value = mean(acceptability), .groups = "drop") %>%
  ggplot(aes(n_words, value, color = condition_label, fill = condition_label)) +
  geom_point(alpha = 0.4, position = position_jitter(width = 0.3, height = 0, seed = 2026)) +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(~ detection_label, nrow = 1) +
  scale_color_manual(name = "Condition", values = condition_colors,
                     aesthetics = c("colour", "fill")) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Acceptability", breaks = seq(1, 7, 2)) +
  coord_cartesian(ylim = c(1, 7))
p_by_detection
ggsave(paste0(figure_prefix, "acceptability_by_detection.png"), p_by_detection,
       width = 11, height = 5.5)

## ---- Figure 3: acceptability by length, all trials vs judged correctly ----
## The all-trials panel shows the ratings without any selection; the second
## keeps only trials whose error judgment was correct.
acceptability_panels <- bind_rows(
  trials %>% mutate(panel = "All trials"),
  trials %>% filter(detected) %>% mutate(panel = "Error judged correctly")) %>%
  mutate(panel = factor(panel, levels = c("All trials", "Error judged correctly")))

p_acceptability <- acceptability_panels %>%
  group_by(panel, condition_label, item, length_level, n_words) %>%
  summarize(value = mean(acceptability), .groups = "drop") %>%
  ggplot(aes(n_words, value, color = condition_label, fill = condition_label)) +
  geom_point(alpha = 0.4, position = position_jitter(width = 0.3, height = 0, seed = 2026)) +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(~ panel, nrow = 1) +
  scale_color_manual(name = "Condition", values = condition_colors,
                     aesthetics = c("colour", "fill")) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Acceptability", breaks = seq(1, 7, 2)) +
  coord_cartesian(ylim = c(1, 7))
p_acceptability
ggsave(paste0(figure_prefix, "acceptability.png"), p_acceptability,
       width = 11, height = 5.5)

## ---- (3) Error-free sentences: is longer still less acceptable? ----
## The communication account's length prediction on its own: among originals
## correctly judged error-free, longer versions should be rated lower.
cat("\n===== (3) Originals judged error-free: acceptability ~ length =====\n")
error_free <- model_trials %>%
  filter(condition_label == "Original", detected) %>%
  mutate(n_words_c = n_words - mean(n_words))
error_free_formulas <- c(
  "acceptability ~ n_words_c + (n_words_c | participant_id) + (n_words_c | item)",
  "acceptability ~ n_words_c + (n_words_c || participant_id) + (n_words_c || item)",
  "acceptability ~ n_words_c + (1 | participant_id) + (1 | item)")
for (f in error_free_formulas) {
  m_error_free <- suppressWarnings(lmer(as.formula(f), data = error_free, control = lmer_ctrl))
  problems <- grep("singular", m_error_free@optinfo$conv$lme4$messages,
                   ignore.case = TRUE, invert = TRUE, value = TRUE)
  if (length(problems) == 0) break
}
cat("random effects:", f, "\n")
print(round(summary(m_error_free)$coefficients, 4))

## ---- (4) The paper's three accounts, plus two length baselines ----
## Each account in both forms the paper uses: shape-constrained GAMs over
## sentence means, and linear mixed models over trials (REML = FALSE, because
## the models differ in their fixed effects). Two extra GAMs test whether a
## plain length term can stand in for the communication model: gradient +
## length, and length alone.
##
## Each model with an information term is fitted with length and again with
## sentence surprisal in its place (same smooths, same random effects). Binary
## and gradient grammaticality have no information term: one fit of each
## serves under both.
##
## Three views of the trials. Restricting to correctly judged trials selects on
## a variable correlated with the rating, so "All trials" is the baseline the
## others have to be read against.
datasets <- list(
  "All trials" = model_trials,
  "Judged correctly" = model_trials %>% filter(detected),
  "Judged correctly + flagged originals" = model_trials %>%
    filter(detected | condition_label == "Original"))

cat("\n===== (4) Model comparison =====\n")
bind_rows(datasets, .id = "dataset") %>%
  group_by(dataset = factor(dataset, levels = names(datasets))) %>%
  summarize(trials = n(), originals = sum(condition_label == "Original"),
            sentences = n_distinct(sentence_id), .groups = "drop") %>%
  print()

scam_models <- list()
lmer_models <- list()
for (dataset_name in names(datasets)) {
  d <- datasets[[dataset_name]] %>% droplevels()
  cat("fitting:", dataset_name, "\n")

  ## sentence means for the GAMs; item is the base item, so the by-item random
  ## slopes are estimated across its short, medium and long versions
  d_means <- d %>%
    group_by(item, length_level, cond, n_words, dl_distance, weighted_dl_distance,
             sum_surprisal, weighted_dl_distance_surprisal, grammatical) %>%
    summarize(response_num = mean(acceptability), .groups = "drop") %>%
    mutate(grammatical_numeric = ifelse(grammatical == "grammatical", 0.5, -0.5))

  binary_gam <- scam(
    response_num ~ grammatical_numeric + s(item, bs = "re") +
      s(item, grammatical_numeric, bs = "re"), data = d_means)
  gradient_gam <- scam(
    response_num ~ s(dl_distance, bs = "mpd") + s(item, bs = "re") +
      s(item, dl_distance, bs = "re"), data = d_means)
  scam_models[[dataset_name]] <- list(
    "Length" = list(
      "Binary grammaticality" = binary_gam,
      "Gradient grammaticality" = gradient_gam,
      "Communication-based" = scam(
        response_num ~ s(weighted_dl_distance, bs = "mpd") + s(n_words, bs = "mpd") +
          s(item, bs = "re") + s(item, weighted_dl_distance, bs = "re"), data = d_means),
      "Gradient + length" = scam(
        response_num ~ s(dl_distance, bs = "mpd") + s(n_words, bs = "mpd") +
          s(item, bs = "re") + s(item, dl_distance, bs = "re"), data = d_means),
      "Length only" = scam(
        response_num ~ s(n_words, bs = "mpd") + s(item, bs = "re"), data = d_means)),
    "Surprisal" = list(
      "Binary grammaticality" = binary_gam,
      "Gradient grammaticality" = gradient_gam,
      "Communication-based" = scam(
        response_num ~ s(weighted_dl_distance_surprisal, bs = "mpd") + s(sum_surprisal, bs = "mpd") +
          s(item, bs = "re") + s(item, weighted_dl_distance_surprisal, bs = "re"), data = d_means),
      "Gradient + surprisal" = scam(
        response_num ~ s(dl_distance, bs = "mpd") + s(sum_surprisal, bs = "mpd") +
          s(item, bs = "re") + s(item, dl_distance, bs = "re"), data = d_means),
      "Surprisal only" = scam(
        response_num ~ s(sum_surprisal, bs = "mpd") + s(item, bs = "re"), data = d_means)))

  ## lme4::lmer rather than lmerTest's: only the log-likelihoods are used here,
  ## and lmerTest's extra Hessian for its p-values can warn on singular fits
  ## without affecting the fit itself.
  binary_lmer <- lme4::lmer(
    acceptability ~ grammatical + (grammatical | participant_id) + (grammatical | item),
    data = d, REML = FALSE, control = lmer_ctrl)
  gradient_lmer <- lme4::lmer(
    acceptability ~ dl_distance_z + (dl_distance_z | participant_id) + (dl_distance_z | item),
    data = d, REML = FALSE, control = lmer_ctrl)
  lmer_models[[dataset_name]] <- list(
    "Length" = list(
      "Binary grammaticality" = binary_lmer,
      "Gradient grammaticality" = gradient_lmer,
      "Communication-based" = lme4::lmer(
        acceptability ~ weighted_dl_distance_z + n_words_z +
          (weighted_dl_distance_z + n_words_z | participant_id) +
          (weighted_dl_distance_z + n_words_z | item),
        data = d, REML = FALSE, control = lmer_ctrl)),
    "Surprisal" = list(
      "Binary grammaticality" = binary_lmer,
      "Gradient grammaticality" = gradient_lmer,
      "Communication-based" = lme4::lmer(
        acceptability ~ weighted_dl_distance_surprisal_z + sum_surprisal_z +
          (weighted_dl_distance_surprisal_z + sum_surprisal_z | participant_id) +
          (weighted_dl_distance_surprisal_z + sum_surprisal_z | item),
        data = d, REML = FALSE, control = lmer_ctrl)))
}

## Log-likelihoods. lmer's logLik() is already marginal over the random
## effects. scam's is conditional on the fitted random effects, so for the GAMs
## the marginal log-likelihood is computed directly: the fixed part of the fit
## with the random effects switched off, and a covariance that adds each random
## effect's variance to the residual variance (the same calculation as
## marginal_logLik() in setup.R).
likelihoods <- tibble()
for (dataset_name in names(scam_models)) {
  for (information in names(scam_models[[dataset_name]])) {
    for (model_name in names(scam_models[[dataset_name]][[information]])) {
      m <- scam_models[[dataset_name]][[information]][[model_name]]
      y <- m$y
      re_index <- which(sapply(m$smooth, inherits, "random.effect"))
      re_labels <- sapply(m$smooth[re_index], "[[", "label")
      mu <- as.numeric(predict(m, exclude = re_labels))
      X <- predict(m, type = "lpmatrix")
      V <- diag(m$sig2, length(y))
      for (k in re_index) {
        Z <- X[, m$smooth[[k]]$first.para:m$smooth[[k]]$last.para, drop = FALSE]
        re_variance <- m$sig2 / (m$sp[[m$smooth[[k]]$label]] * diag(m$smooth[[k]]$S[[1]]))
        V <- V + Z %*% (re_variance * t(Z))
      }
      R <- chol(V)
      marginal <- -0.5 * (length(y) * log(2 * pi) + 2 * sum(log(diag(R))) +
                            sum(backsolve(R, y - mu, transpose = TRUE)^2))
      likelihoods <- bind_rows(likelihoods, tibble(
        dataset = dataset_name, family = "Shape-constrained GAM (sentence means)",
        information = information, model = model_name, converged = m$conv,
        log_likelihood = marginal))
    }
    for (model_name in names(lmer_models[[dataset_name]][[information]])) {
      m <- lmer_models[[dataset_name]][[information]][[model_name]]
      problems <- grep("singular", m@optinfo$conv$lme4$messages,
                       ignore.case = TRUE, invert = TRUE, value = TRUE)
      likelihoods <- bind_rows(likelihoods, tibble(
        dataset = dataset_name, family = "Linear mixed model (trials)",
        information = information, model = model_name, converged = length(problems) == 0,
        log_likelihood = as.numeric(logLik(m))))
    }
  }
}

## Likelihoods are compared only within a family, a dataset and an information
## measure. The paper's three accounts are ranked among themselves; the two
## baselines are shown against the same best model.
likelihoods <- likelihoods %>%
  mutate(dataset = factor(dataset, levels = names(datasets)),
         information = factor(information, levels = c("Length", "Surprisal")),
         paper_model = model %in% names(model_colors)) %>%
  group_by(dataset, family, information) %>%
  mutate(likelihood_diff = log_likelihood - max(log_likelihood[paper_model])) %>%
  ungroup()

cat("\n---- log-likelihoods (difference from the best of the three accounts) ----\n")
likelihoods %>%
  arrange(family, dataset, information, desc(log_likelihood)) %>%
  mutate(log_likelihood = round(log_likelihood, 1), likelihood_diff = round(likelihood_diff, 1)) %>%
  select(family, dataset, information, model, converged, log_likelihood, likelihood_diff) %>%
  as.data.frame() %>%
  print(right = FALSE)
write_csv(likelihoods, paste0(sub("^figures/", "outputs/", figure_prefix), "likelihoods.csv"))

## ---- Figure 4: likelihood differences, shape-constrained GAMs ----
## The paper's three accounts only, and the two views that matter most, with
## length and with surprisal as the information measure; the linear mixed
## models, the baselines and the "+ flagged originals" view are in the printed
## table above.
p_likelihoods <- likelihoods %>%
  filter(family == "Shape-constrained GAM (sentence means)", paper_model,
         dataset %in% c("All trials", "Judged correctly")) %>%
  mutate(model = factor(model, levels = names(model_colors)),
         dataset = recode(dataset, "All trials" = "Binary error detection\n(all responses)",
                          "Judged correctly" = "Binary error detection\n(judged correctly)")) %>%
  ggplot(aes(x = likelihood_diff, y = dataset, color = model, shape = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
  geom_point(size = 3) +
  facet_wrap(~ information) +
  scale_color_manual(name = "", values = model_colors) +
  scale_shape_manual(name = "", values = c(16, 17, 15)) +
  scale_x_continuous(name = "Log-likelihood difference from the best model") +
  scale_y_discrete(name = "", limits = rev) +
  guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  theme(legend.position = "top", legend.text = element_text(size = 15),
        axis.text.y = element_text(size = 13, lineheight = 0.9))
p_likelihoods
ggsave(paste0(figure_prefix, "likelihood_differences.png"), p_likelihoods,
       width = 16, height = 5)

## ---- (5) Does the same swap cost less in a longer sentence? ----
## Within an item, every version of the same length has the same number of
## words -- a swap only reorders. So the penalty mean(original) - mean(swapped),
## taken item by item at each length, cancels any plain effect of length.
## Length aversion, binary and gradient grammaticality all predict a flat
## penalty; only the communication account predicts it shrinks with length.
## All trials, no selection on the error judgment.
cat("\n===== (5) Within-item swap penalty by length =====\n")
penalties <- trials %>%
  group_by(item, length_level, n_words, condition_label) %>%
  summarize(mean_rating = mean(acceptability), .groups = "drop") %>%
  pivot_wider(names_from = condition_label, values_from = mean_rating) %>%
  transmute(item, length_level, n_words,
            `1 swap` = Original - `1 swap`,
            `3 swaps` = Original - `3 swaps`) %>%
  pivot_longer(c(`1 swap`, `3 swaps`), names_to = "condition_label", values_to = "penalty") %>%
  filter(!is.na(penalty))

cat("penalty slope per extra word (sentence-level lm):\n")
penalties %>%
  group_by(condition_label) %>%
  summarize(sentences = n(),
            slope = coef(lm(penalty ~ n_words))[2],
            p = summary(lm(penalty ~ n_words))$coefficients[2, 4],
            .groups = "drop") %>%
  print()

## The same quantity at the trial level: the condition x length interaction.
## A positive interaction means the swap penalty shrinks as the sentence grows.
interaction_data <- model_trials %>%
  mutate(condition = factor(recode(as.character(condition_label), "Original" = "original",
                                   "1 swap" = "swap_1", "3 swaps" = "swap_3"),
                            levels = c("original", "swap_1", "swap_3")),
         n_words_c = n_words - mean(n_words))
interaction_formulas <- c(
  "acceptability ~ condition * n_words_c + (condition * n_words_c | participant_id) + (condition * n_words_c | item)",
  "acceptability ~ condition * n_words_c + (condition + n_words_c | participant_id) + (condition + n_words_c | item)",
  "acceptability ~ condition * n_words_c + (condition + n_words_c || participant_id) + (condition + n_words_c || item)",
  "acceptability ~ condition * n_words_c + (1 | participant_id) + (1 | item)")
for (f in interaction_formulas) {
  m_interaction <- suppressWarnings(lmer(as.formula(f), data = interaction_data, control = lmer_ctrl))
  problems <- grep("singular", m_interaction@optinfo$conv$lme4$messages,
                   ignore.case = TRUE, invert = TRUE, value = TRUE)
  if (length(problems) == 0) break
}
cat("\n---- acceptability ~ condition x length (trial level) ----\n")
cat("random effects:", f, "\n")
print(round(summary(m_interaction)$coefficients, 4))

## ---- Figure 5: the swap penalty against length ----
p_penalty <- penalties %>%
  ggplot(aes(n_words, penalty, color = condition_label, fill = condition_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(alpha = 0.4, position = position_jitter(width = 0.3, height = 0, seed = 2026)) +
  geom_smooth(method = "lm", formula = y ~ x) +
  scale_color_manual(name = "Condition", values = condition_colors[c("1 swap", "3 swaps")],
                     aesthetics = c("colour", "fill")) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Acceptability cost of the swap\n(original - swapped, same item)") +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
  theme(legend.position = "top")
p_penalty
ggsave(paste0(figure_prefix, "swap_penalty_by_length.png"), p_penalty,
       width = 9, height = 5.5)
