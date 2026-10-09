## =============================================================================
## SI Section 4.4: Acceptability alone on the artificial materials
##
## en_artificial_acceptability_ratings. The 27 constructed items of sections 4.1-4.3,
## each in a short, medium and long version, each original, with 1 swap or with
## 3 swaps: 9 conditions, every participant sees every item once. Participants
## only rate each sentence (1-7), with the original English experiment's prompt
## and scale -- no error question, no exchange identification task, no comprehension question.
##
##   (1) Ratings by condition and length, and the figure.
##   (2) Linear mixed model: condition x length interaction.
##   (3) The paper's three accounts as shape-constrained GAMs, with sentence
##       length and with sentence surprisal as the information measure.
##
## Figures: figures/artificial_acceptability_ratings_*.png
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

## ---- Materials: one row per sentence version ----
materials <- read_csv("../data/en_artificial_acceptability_ratings/materials_en_artificial_acceptability_ratings.csv",
                      show_col_types = FALSE) %>%
  transmute(item, cond = condition, sentence,
            length_level = factor(str_extract(condition, "[^_]+$"),
                                  levels = c("short", "medium", "long")),
            condition_label = factor(recode(str_extract(condition, "^[^_]+"),
                                            "original" = "Original",
                                            "swap-1" = "1 swap",
                                            "swap-3" = "3 swaps"),
                                     levels = names(condition_colors)))

## Damerau-Levenshtein distance from each version to its original, computed as
## in data/preprocess_stimuli.R and section 4.2.
sentence_letters <- materials %>%
  separate_rows(sentence, sep = " ") %>%
  mutate(word = str_remove_all(str_to_lower(sentence), "[^a-z0-9]")) %>%
  group_by(item, length_level) %>%
  mutate(letter = c(letters, LETTERS)[match(word, unique(word))]) %>%
  group_by(item, length_level, cond, condition_label) %>%
  summarize(letter_string = paste(letter, collapse = ""), n_words = n(), .groups = "drop")

## Sentence surprisal, the main analysis's information measure: the summed
## Goldfish surprisal (bits) of the original at each length, from
## data/slor_morcela_scores.csv, as in section 4.2.
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
         grammatical_numeric = ifelse(dl_distance == 0, 0.5, -0.5)) %>%
  select(item, cond, length_level, condition_label, n_words, dl_distance,
         weighted_dl_distance, sum_surprisal, weighted_dl_distance_surprisal,
         grammatical_numeric)
stopifnot(!anyNA(sentences))

cat("DL distance by condition (number of sentence versions):\n")
sentences %>% count(condition_label, dl_distance) %>% print()

## ---- Responses ----
responses <- read_csv("../data/en_artificial_acceptability_ratings/en_artificial_acceptability_ratings_data.csv",
                      show_col_types = FALSE)

cat("\n---- export integrity ----\n")
cat("participants:", n_distinct(responses$participant_id), "\n")
responses %>% count(participant_id, name = "trials") %>% count(trials, name = "participants") %>% print()
cat("(27 trials = complete; Surveyor writes no row for a skipped question)\n")
responses %>% distinct(participant_id, native_language) %>% count(native_language) %>% print()

## The study is of native English speakers: drop anyone who answered No.
## A blank answer is kept, as in sections 4.1 and 4.2.
non_native <- responses %>%
  distinct(participant_id, native_language) %>%
  filter(native_language %in% "No") %>%
  pull(participant_id)
cat("excluded as non-native:", length(non_native), "participant(s)\n")

trials <- responses %>%
  filter(!participant_id %in% non_native) %>%
  transmute(participant_id, item, cond,
            acceptability = as.numeric(substr(response, 1, 1))) %>%
  left_join(sentences, by = c("item", "cond"))

cat("participants analysed:", n_distinct(trials$participant_id), "| trials:", nrow(trials), "\n")

## ---- (1) Ratings by condition and length ----
cat("\n===== (1) Mean rating by condition and length =====\n")
trials %>%
  group_by(condition_label, length_level) %>%
  summarize(mean_rating = round(mean(acceptability), 2), .groups = "drop") %>%
  pivot_wider(names_from = length_level, values_from = mean_rating) %>%
  print()

## Figure: points are sentence means (one item at one length), jittered
## sideways; lines are linear fits to the trial-level ratings.
p_acceptability <- trials %>%
  group_by(condition_label, item, length_level, n_words) %>%
  summarize(value = mean(acceptability), .groups = "drop") %>%
  ggplot(aes(n_words, value, color = condition_label, fill = condition_label)) +
  geom_point(alpha = 0.4, position = position_jitter(width = 0.3, height = 0, seed = 2026)) +
  geom_smooth(data = trials, aes(y = acceptability), method = "lm", formula = y ~ x) +
  scale_color_manual(name = "Condition", values = condition_colors,
                     aesthetics = c("colour", "fill")) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Acceptability", breaks = seq(1, 7, 2)) +
  coord_cartesian(ylim = c(1, 7))
p_acceptability
ggsave("figures/artificial_acceptability_ratings_by_length.png", p_acceptability, width = 9, height = 5.5)

## ---- (2) Linear mixed model: condition x length ----
## Original is the reference level and length is centred, so the n_words_c
## term is the length slope of the originals, and each condition x n_words_c
## term is how much less steeply (if positive) that swap condition's ratings
## change with length. Length varies within item, so it takes by-item slopes.
## The fullest random-effect structure that converges is used.
cat("\n===== (2) acceptability ~ condition x length =====\n")
model_trials <- trials %>%
  mutate(participant_id = factor(participant_id), item = factor(item),
         condition = factor(recode(as.character(condition_label), "Original" = "original",
                                   "1 swap" = "swap_1", "3 swaps" = "swap_3"),
                            levels = c("original", "swap_1", "swap_3")),
         n_words_c = n_words - mean(n_words))
lmer_formulas <- c(
  "acceptability ~ condition * n_words_c + (condition * n_words_c | participant_id) + (condition * n_words_c | item)",
  "acceptability ~ condition * n_words_c + (condition + n_words_c | participant_id) + (condition + n_words_c | item)",
  "acceptability ~ condition * n_words_c + (condition + n_words_c || participant_id) + (condition + n_words_c || item)",
  "acceptability ~ condition * n_words_c + (1 | participant_id) + (1 | item)")
for (f in lmer_formulas) {
  m_lmer <- suppressWarnings(lmer(as.formula(f), data = model_trials, control = lmer_ctrl))
  problems <- grep("singular", m_lmer@optinfo$conv$lme4$messages,
                   ignore.case = TRUE, invert = TRUE, value = TRUE)
  if (length(problems) == 0) break
}
cat("random effects:", f, "\n")
print(round(summary(m_lmer)$coefficients, 4))

## ---- (3) The paper's three accounts: shape-constrained GAMs ----
## Fit to sentence means, with by-item random intercepts and slopes, compared
## on marginal log-likelihood (scam's own logLik() is conditional on the fitted
## random effects, so the random effects are integrated out directly, as in
## section 4.2 and marginal_logLik() in setup.R).
##
## The communication model is fitted with length and again with sentence
## surprisal in its place (same smooths, same random effects); binary and
## gradient grammaticality have no information term, so one fit of each serves
## under both, and the accounts are compared within each information measure.
cat("\n===== (3) Shape-constrained GAMs =====\n")
sentence_means <- trials %>%
  group_by(item, length_level, cond, n_words, dl_distance, weighted_dl_distance,
           sum_surprisal, weighted_dl_distance_surprisal, grammatical_numeric) %>%
  summarize(response_num = mean(acceptability), .groups = "drop") %>%
  mutate(item = factor(item))

binary_gam <- scam(
  response_num ~ grammatical_numeric + s(item, bs = "re") +
    s(item, grammatical_numeric, bs = "re"), data = sentence_means)
gradient_gam <- scam(
  response_num ~ s(dl_distance, bs = "mpd") + s(item, bs = "re") +
    s(item, dl_distance, bs = "re"), data = sentence_means)
scam_models <- list(
  "Length" = list(
    "Binary grammaticality" = binary_gam,
    "Gradient grammaticality" = gradient_gam,
    "Communication-based" = scam(
      response_num ~ s(weighted_dl_distance, bs = "mpd") + s(n_words, bs = "mpd") +
        s(item, bs = "re") + s(item, weighted_dl_distance, bs = "re"), data = sentence_means)),
  "Surprisal" = list(
    "Binary grammaticality" = binary_gam,
    "Gradient grammaticality" = gradient_gam,
    "Communication-based" = scam(
      response_num ~ s(weighted_dl_distance_surprisal, bs = "mpd") + s(sum_surprisal, bs = "mpd") +
        s(item, bs = "re") + s(item, weighted_dl_distance_surprisal, bs = "re"),
      data = sentence_means)))

likelihoods <- tibble()
for (information in names(scam_models)) {
  for (model_name in names(scam_models[[information]])) {
    m <- scam_models[[information]][[model_name]]
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
    likelihoods <- bind_rows(likelihoods, tibble(information = information, model = model_name,
                                                 converged = m$conv, log_likelihood = marginal))
  }
}
likelihoods <- likelihoods %>%
  mutate(information = factor(information, levels = c("Length", "Surprisal")),
         model = factor(model, levels = names(model_colors))) %>%
  group_by(information) %>%
  mutate(likelihood_diff = log_likelihood - max(log_likelihood)) %>%
  ungroup()
print(likelihoods %>% mutate(across(where(is.numeric), ~ round(.x, 1))))
write_csv(likelihoods, "outputs/artificial_acceptability_ratings_likelihoods.csv")

## Figure: likelihood differences from the best account, with length and with
## surprisal as the information measure.
p_likelihoods <- likelihoods %>%
  ggplot(aes(x = likelihood_diff, y = "Acceptability only", color = model, shape = model)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
  geom_point(size = 3) +
  facet_wrap(~ information) +
  scale_color_manual(name = "", values = model_colors) +
  scale_shape_manual(name = "", values = c(16, 17, 15)) +
  scale_x_continuous(name = "Log-likelihood difference from the best model") +
  scale_y_discrete(name = "") +
  guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  theme(legend.position = "top", legend.text = element_text(size = 15),
        axis.text.y = element_text(size = 13))
p_likelihoods
ggsave("figures/artificial_acceptability_ratings_likelihood_differences.png", p_likelihoods,
       width = 16, height = 4)
