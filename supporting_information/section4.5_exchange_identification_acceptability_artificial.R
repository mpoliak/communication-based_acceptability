## =============================================================================
## SI Section 4.5: Exchange identification task + acceptability on the artificial materials
##
## en_artificial_exchange_identification_acceptability. The 27 constructed items, each in a
## short, medium and long version, each original, with 1 swap or with 3 swaps:
## 9 conditions, every participant sees every item once. On each trial the
## participant first clicks the words they think are misplaced (the task of
## section 4.1) and then, on the next screen, rates the sentence (1-7, the
## prompt and scale of section 4.4).
##
##   A. the exchange identification task: outcomes by length, and length models (as section 4.1).
##   B. Acceptability: the figure, the condition x length lmer and the three
##      shape-constrained GAMs (as section 4.4), with sentence length and with
##      sentence surprisal as the information measure.
##   C. Acceptability conditional on the same trial's click: all trials vs
##      trials where the click was correct, and correct + wrongly flagged
##      originals (as section 4.2 did with the error question).
##
## A click is correct when an original has nothing clicked, or a swapped
## sentence has a swapped word clicked and nothing away from the swapped
## words (section 4.1's "located"; for 3 swaps this includes marking only some
## of the swapped regions).
##
## Mixed models use the fullest random-effect structure that converges: each
## is tried from the maximal structure down, stopping at the first without a
## convergence failure (singular fits are ignored, as throughout).
##
## Figures: figures/artificial_exchange_identification_acceptability_*.png and
## figures/artificial_exchange_identification_ratings*.png; numbers in
## outputs/artificial_exchange_identification_acceptability_likelihoods.csv
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
## The outcome palette is built with setNames: a "±" written as c("..." = colour)
## is mangled in a non-UTF-8 locale and the level silently falls back to grey.
condition_colors <- c("Original" = "#EB7AA3", "1 swap" = "#1E88E5", "3 swaps" = "#FFC107")
model_colors <- c("Binary grammaticality" = "#EB7AA3",
                  "Gradient grammaticality" = "#1E88E5",
                  "Communication-based" = "#004D40")
outcome_colors <- setNames(
  c("#004D40", "#5FA894", "#FFC107", "#EB7AA3", "#BBBBBB"),
  c("Correct\n(target \u00b1 adjacent words)", "Correct, partly\n(some swaps not marked)",
    "Mixed\n(target + distant words)", "Wrong\n(other words only)",
    "Wrong\n(clicked nothing)"))
outcome_linetypes <- setNames(c("solid", "dashed", "solid", "solid", "solid"),
                              names(outcome_colors))
original_colors <- c("Correct\n(clicked nothing)" = "#004D40",
                     "False alarm\n(1+ words marked)" = "#EB7AA3")
lmer_ctrl <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
glmer_ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

## ---- Materials ----
## The click rows carry the sentences (the rating rows have none).
materials <- read_csv("../data/en_artificial_exchange_identification_acceptability/materials_en_artificial_exchange_identification_acceptability.csv",
                      show_col_types = FALSE) %>%
  filter(experiment == "annot") %>%
  transmute(item, cond = condition, sentence,
            length_level = factor(str_extract(condition, "[^_]+$"),
                                  levels = c("short", "medium", "long")),
            condition_label = factor(recode(str_extract(condition, "^[^_]+"),
                                            "original" = "Original",
                                            "swap-1" = "1 swap",
                                            "swap-3" = "3 swaps"),
                                     levels = names(condition_colors)))

## One row per word. A word is a target (displaced) when it differs from the
## word at the same position in that item's original of the same length;
## consecutive targets form a run (one per swapped pair, or per chain where two
## swaps share a word).
words <- materials %>%
  separate_rows(sentence, sep = " ") %>%
  rename(word = sentence) %>%
  group_by(item, cond) %>%
  mutate(position = row_number(), n_words = n()) %>%
  ungroup()
words <- words %>%
  left_join(words %>%
              filter(condition_label == "Original") %>%
              select(item, length_level, position, original_word = word),
            by = c("item", "length_level", "position")) %>%
  mutate(is_target = word != original_word) %>%
  group_by(item, cond) %>%
  mutate(target_run = ifelse(is_target,
                             cumsum(is_target & !lag(is_target, default = FALSE)),
                             NA)) %>%
  ungroup()

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

## Damerau-Levenshtein distance from each version to its original, computed as
## in data/preprocess_stimuli.R: one letter per distinct word of the item at
## that length, then stringdist(method = "dl").
sentences <- words %>%
  mutate(bare = str_remove_all(str_to_lower(word), "[^a-z0-9]")) %>%
  group_by(item, length_level) %>%
  mutate(letter = c(letters, LETTERS)[match(bare, unique(bare))]) %>%
  group_by(item, length_level, cond, condition_label, n_words) %>%
  summarize(letter_string = paste(letter, collapse = ""), .groups = "drop")
sentences <- sentences %>%
  left_join(sentences %>%
              filter(condition_label == "Original") %>%
              select(item, length_level, original_string = letter_string),
            by = c("item", "length_level")) %>%
  left_join(sentence_surprisal, by = c("item", "length_level")) %>%
  mutate(dl_distance = stringdist(letter_string, original_string, method = "dl"),
         weighted_dl_distance = dl_distance / n_words,
         weighted_dl_distance_surprisal = dl_distance / sum_surprisal,
         grammatical_numeric = ifelse(dl_distance == 0, 0.5, -0.5)) %>%
  select(item, cond, n_words, dl_distance, weighted_dl_distance, sum_surprisal,
         weighted_dl_distance_surprisal, grammatical_numeric)
stopifnot(!anyNA(sentences))

cat("DL distance by condition (number of sentence versions):\n")
sentences %>% left_join(materials %>% select(item, cond, condition_label), by = c("item", "cond")) %>%
  count(condition_label, dl_distance) %>% print()

## ---- Responses ----
responses <- read_csv("../data/en_artificial_exchange_identification_acceptability/en_artificial_exchange_identification_acceptability_data.csv",
                      show_col_types = FALSE)

cat("\n---- export integrity ----\n")
cat("participants:", n_distinct(responses$participant_id), "\n")
responses %>% count(participant_id, experiment) %>%
  count(experiment, n, name = "participants") %>% print()
cat("(27 click rows each = complete: the exchange identification task always writes a row;",
    "a skipped rating writes none)\n")
responses %>% distinct(participant_id, native_language) %>% count(native_language) %>% print()

## The study is of native English speakers: drop anyone who answered No.
## A blank answer is kept, as in sections 4.1, 4.2 and 4.4.
non_native <- responses %>%
  distinct(participant_id, native_language) %>%
  filter(native_language %in% "No") %>%
  pull(participant_id)
cat("excluded as non-native:", length(non_native), "participant(s)\n")
responses <- responses %>% filter(!participant_id %in% non_native)

## Participants who also took the binary error detection + acceptability ratings
## study saw the same sentences twice; they are excluded from both.
took_both <- intersect(
  str_trim(responses$participant_id),
  str_trim(read_csv("../data/en_artificial_binary_error_detection/en_artificial_binary_error_detection_data.csv",
                    show_col_types = FALSE)$participant_id))
cat("excluded for also taking binary error detection + acceptability ratings:",
    length(took_both), "participant(s)\n")
responses <- responses %>% filter(!str_trim(participant_id) %in% took_both)

ratings <- responses %>%
  filter(experiment == "exp_b") %>%
  transmute(participant_id, item, cond,
            acceptability = as.numeric(substr(response, 1, 1)))

trials <- responses %>%
  filter(experiment == "annot") %>%
  mutate(trial = row_number()) %>%
  select(trial, participant_id, item, cond, response) %>%
  left_join(materials %>% select(item, cond, length_level, condition_label),
            by = c("item", "cond")) %>%
  left_join(ratings, by = c("participant_id", "item", "cond"))

cat("click trials:", nrow(trials), "| of which rated:", sum(!is.na(trials$acceptability)),
    "| ratings without a click trial:",
    nrow(anti_join(ratings, trials, by = c("participant_id", "item", "cond"))), "\n")

## One row per clicked word. Words are percent-escaped for % , and : so they
## survive the position:word encoding, hence URLdecode.
clicks <- trials %>%
  filter(response != "NONE") %>%
  separate_rows(response, sep = ",") %>%
  separate(response, into = c("position", "clicked_word"), sep = ":", extra = "merge") %>%
  transmute(trial, item, cond, position = as.integer(position),
            clicked_word = map_chr(clicked_word, URLdecode))

clicks %>%
  left_join(words %>% select(item, cond, position, word), by = c("item", "cond", "position")) %>%
  summarize(clicks = n(), matching_the_stimulus = sum(clicked_word == word, na.rm = TRUE),
            out_of_range = sum(is.na(word))) %>%
  print()

## Every word of the sentence each participant saw, marked clicked or not.
trial_words <- trials %>%
  select(trial, item, cond) %>%
  inner_join(words %>% select(item, cond, position, is_target, target_run),
             by = c("item", "cond"), relationship = "many-to-many") %>%
  left_join(clicks %>% transmute(trial, position, clicked = TRUE),
            by = c("trial", "position")) %>%
  mutate(clicked = replace_na(clicked, FALSE))

## ---- Scoring each click trial (as in section 4.1) ----
## Located: at least one displaced word clicked, and every other click in an
## unbroken run with the displaced words. Checked run by run: split the
## clicked-or-target positions into runs of consecutive positions and flag any
## run with no target in it as a stray click.
stray_clicks <- trial_words %>%
  filter(clicked | is_target) %>%
  arrange(trial, position) %>%
  group_by(trial) %>%
  mutate(run = cumsum(position != lag(position, default = -1L) + 1L)) %>%
  group_by(trial, run) %>%
  summarize(stray = !any(is_target), .groups = "drop") %>%
  group_by(trial) %>%
  summarize(any_stray = any(stray), .groups = "drop")

trials <- trials %>%
  left_join(trial_words %>%
              group_by(trial) %>%
              summarize(n_clicked = sum(clicked),
                        n_on_target = sum(clicked & is_target),
                        runs_total = n_distinct(target_run, na.rm = TRUE),
                        runs_found = n_distinct(target_run[clicked & is_target], na.rm = TRUE),
                        .groups = "drop"),
            by = "trial") %>%
  left_join(stray_clicks, by = "trial") %>%
  left_join(sentences, by = c("item", "cond")) %>%
  mutate(flagged = n_clicked > 0,
         located = n_on_target > 0 & !any_stray,
         click_correct = ifelse(condition_label == "Original", !flagged, located),
         outcome = case_when(
           !flagged         ~ "Wrong\n(clicked nothing)",
           located & runs_found < runs_total ~ "Correct, partly\n(some swaps not marked)",
           located          ~ "Correct\n(target \u00b1 adjacent words)",
           n_on_target == 0 ~ "Wrong\n(other words only)",
           TRUE             ~ "Mixed\n(target + distant words)"),
         outcome = factor(outcome, levels = names(outcome_colors)))

originals <- trials %>% filter(condition_label == "Original")
swapped <- trials %>% filter(condition_label != "Original")


## =============================================================================
## A. the exchange identification task
## =============================================================================
cat("\n\n===== A. Exchange identification task =====\n")
cat("false-alarm rate (any word clicked on an original):", round(mean(originals$flagged), 3), "\n")
trials %>%
  group_by(condition_label, length_level) %>%
  summarize(trials = n(), clicked_anything = round(mean(flagged), 3),
            click_correct = round(mean(click_correct), 3),
            all_swaps_marked = round(mean(runs_found == runs_total & located), 3),
            .groups = "drop") %>%
  print(n = Inf)

## Length models, as in section 4.1: locating a single swap (1-swap trials) and
## false alarms (original trials). Log-odds per SD of length.
click_models <- trials %>%
  mutate(participant_id = factor(participant_id), item = factor(item),
         n_words_z = as.numeric(scale(n_words)),
         located = as.integer(located), flagged = as.integer(flagged))
cat("\nn_words_z is per SD of length (SD =", round(sd(click_models$n_words), 1), "words)\n")

located_formulas <- c(
  "located ~ n_words_z + (n_words_z | participant_id) + (n_words_z | item)",
  "located ~ n_words_z + (n_words_z || participant_id) + (n_words_z || item)",
  "located ~ n_words_z + (1 | participant_id) + (1 | item)")
for (f in located_formulas) {
  m_located <- suppressWarnings(glmer(as.formula(f), family = binomial, control = glmer_ctrl,
                                      data = click_models %>% filter(condition_label == "1 swap")))
  problems <- grep("singular", m_located@optinfo$conv$lme4$messages,
                   ignore.case = TRUE, invert = TRUE, value = TRUE)
  if (length(problems) == 0) break
}
cat("\n---- located the swap ~ length (1-swap trials) ----\nrandom effects:", f, "\n")
print(round(summary(m_located)$coefficients, 4))

false_alarm_formulas <- c(
  "flagged ~ n_words_z + (n_words_z | participant_id) + (n_words_z | item)",
  "flagged ~ n_words_z + (n_words_z || participant_id) + (n_words_z || item)",
  "flagged ~ n_words_z + (1 | participant_id) + (1 | item)")
for (f in false_alarm_formulas) {
  m_false_alarm <- suppressWarnings(glmer(as.formula(f), family = binomial, control = glmer_ctrl,
                                          data = click_models %>% filter(condition_label == "Original")))
  problems <- grep("singular", m_false_alarm@optinfo$conv$lme4$messages,
                   ignore.case = TRUE, invert = TRUE, value = TRUE)
  if (length(problems) == 0) break
}
cat("\n---- false alarm ~ length (original trials) ----\nrandom effects:", f, "\n")
cat("(false alarms:", sum(click_models$flagged[click_models$condition_label == "Original"]),
    "-- check the estimates are sensible; with very few events a converged fit can be degenerate)\n")
print(round(summary(m_false_alarm)$coefficients, 4))

## ---- Figure A1: click outcomes on swapped sentences, by length ----
## Points are sentences (one item at one length), jittered; lines are linear
## fits to the trial-level 0/1 data. "Correct, partly" cannot happen with a
## single swap, so it is left out of the 1-swap panel.
outcomes_long <- swapped %>%
  select(condition_label, item, length_level, n_words, outcome) %>%
  expand_grid(level = names(outcome_colors)) %>%
  filter(!(condition_label == "1 swap" & level == "Correct, partly\n(some swaps not marked)")) %>%
  mutate(y = as.integer(as.character(outcome) == level),
         level = factor(level, levels = names(outcome_colors)))

p_outcomes_by_length <- ggplot(mapping = aes(x = n_words, color = level, fill = level,
                                             linetype = level)) +
  geom_point(data = outcomes_long %>%
               group_by(condition_label, level, item, length_level, n_words) %>%
               summarize(value = mean(y), .groups = "drop"),
             aes(y = value), alpha = 0.4,
             position = position_jitter(width = 0.3, height = 0.02, seed = 2026)) +
  geom_smooth(data = outcomes_long, aes(y = y), method = "lm", formula = y ~ x) +
  facet_wrap(~ condition_label) +
  scale_color_manual(name = "", values = outcome_colors, aesthetics = c("colour", "fill")) +
  scale_linetype_manual(name = "", values = outcome_linetypes) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Proportion of swapped trials", labels = scales::percent,
                     breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1)) +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1),
         linetype = guide_legend(nrow = 1)) +
  theme(legend.position = "top", legend.text = element_text(size = 11),
        legend.key.width = unit(1.4, "cm"))
p_outcomes_by_length
ggsave("figures/artificial_exchange_identification_acceptability_outcomes_by_length.png", p_outcomes_by_length,
       width = 13, height = 5.5)

## ---- Figure A2: click outcomes on original sentences, by length ----
original_long <- originals %>%
  mutate(outcome = ifelse(flagged, "False alarm\n(1+ words marked)", "Correct\n(clicked nothing)")) %>%
  select(item, length_level, n_words, outcome) %>%
  expand_grid(level = names(original_colors)) %>%
  mutate(y = as.integer(outcome == level),
         level = factor(level, levels = names(original_colors)))

p_outcomes_original_by_length <- ggplot(mapping = aes(x = n_words, color = level, fill = level)) +
  geom_point(data = original_long %>%
               group_by(level, item, length_level, n_words) %>%
               summarize(value = mean(y), .groups = "drop"),
             aes(y = value), alpha = 0.4,
             position = position_jitter(width = 0.3, height = 0.02, seed = 2026)) +
  geom_smooth(data = original_long, aes(y = y), method = "lm", formula = y ~ x) +
  scale_color_manual(name = "", values = original_colors, aesthetics = c("colour", "fill")) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Proportion of original trials", labels = scales::percent,
                     breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1)) +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
  theme(legend.position = "top", legend.text = element_text(size = 11))
p_outcomes_original_by_length
ggsave("figures/artificial_exchange_identification_acceptability_outcomes_original_by_length.png",
       p_outcomes_original_by_length, width = 11, height = 5.5)


## =============================================================================
## B. Acceptability (all rated trials)
## =============================================================================
cat("\n\n===== B. Acceptability =====\n")
rated <- trials %>%
  filter(!is.na(acceptability)) %>%
  mutate(participant_id = factor(participant_id), item = factor(item),
         condition = factor(recode(as.character(condition_label), "Original" = "original",
                                   "1 swap" = "swap_1", "3 swaps" = "swap_3"),
                            levels = c("original", "swap_1", "swap_3")),
         n_words_c = n_words - mean(n_words),
         click_label = factor(ifelse(click_correct, "Identification correct", "Identification incorrect"),
                              levels = c("Identification correct", "Identification incorrect")))

rated %>%
  group_by(condition_label, length_level) %>%
  summarize(mean_rating = round(mean(acceptability), 2), .groups = "drop") %>%
  pivot_wider(names_from = length_level, values_from = mean_rating) %>%
  print()

## ---- Figure B1: acceptability by length ----
p_acceptability <- rated %>%
  group_by(condition_label, item, length_level, n_words) %>%
  summarize(value = mean(acceptability), .groups = "drop") %>%
  ggplot(aes(n_words, value, color = condition_label, fill = condition_label)) +
  geom_point(alpha = 0.4, position = position_jitter(width = 0.3, height = 0, seed = 2026)) +
  geom_smooth(data = rated, aes(y = acceptability), method = "lm", formula = y ~ x) +
  scale_color_manual(name = "Condition", values = condition_colors,
                     aesthetics = c("colour", "fill")) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Acceptability", breaks = seq(1, 7, 2)) +
  coord_cartesian(ylim = c(1, 7))
p_acceptability
ggsave("figures/artificial_exchange_identification_ratings_by_length.png", p_acceptability,
       width = 9, height = 5.5)

## ---- B2: lmer, condition x length (all rated trials) ----
## Original is the reference and length is centred: n_words_c is the originals'
## length slope; a positive condition x n_words_c term means that swap
## condition's ratings change less steeply (or rise) with length.
lmer_formulas <- c(
  "acceptability ~ condition * n_words_c + (condition * n_words_c | participant_id) + (condition * n_words_c | item)",
  "acceptability ~ condition * n_words_c + (condition + n_words_c | participant_id) + (condition + n_words_c | item)",
  "acceptability ~ condition * n_words_c + (condition + n_words_c || participant_id) + (condition + n_words_c || item)",
  "acceptability ~ condition * n_words_c + (1 | participant_id) + (1 | item)")
for (f in lmer_formulas) {
  m_all <- suppressWarnings(lmer(as.formula(f), data = rated, control = lmer_ctrl))
  problems <- grep("singular", m_all@optinfo$conv$lme4$messages,
                   ignore.case = TRUE, invert = TRUE, value = TRUE)
  if (length(problems) == 0) break
}
cat("\n---- acceptability ~ condition x length, all rated trials ----\nrandom effects:", f, "\n")
print(round(summary(m_all)$coefficients, 4))


## =============================================================================
## C. Acceptability conditional on the click
## =============================================================================
cat("\n\n===== C. Acceptability conditional on the click =====\n")
rated %>%
  group_by(condition_label, length_level) %>%
  ## the rating means come first: naming a result click_correct would overwrite
  ## the column the later lines index by
  summarize(rated = n(),
            rating_click_correct = round(mean(acceptability[click_correct]), 2),
            rating_click_incorrect = round(mean(acceptability[!click_correct]), 2),
            n_click_incorrect = sum(!click_correct),
            share_click_correct = round(mean(click_correct), 3),
            .groups = "drop") %>%
  print(n = Inf)

## ---- Figure C1: acceptability by length, split by the click ----
p_by_click <- rated %>%
  group_by(click_label, condition_label, item, length_level, n_words) %>%
  summarize(value = mean(acceptability), .groups = "drop") %>%
  ggplot(aes(n_words, value, color = condition_label, fill = condition_label)) +
  geom_point(alpha = 0.4, position = position_jitter(width = 0.3, height = 0, seed = 2026)) +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(~ click_label, nrow = 1) +
  scale_color_manual(name = "Condition", values = condition_colors,
                     aesthetics = c("colour", "fill")) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Acceptability", breaks = seq(1, 7, 2)) +
  coord_cartesian(ylim = c(1, 7))
p_by_click
ggsave("figures/artificial_exchange_identification_ratings_by_identification.png", p_by_click,
       width = 11, height = 5.5)

## ---- Figure C2: all trials vs click-correct trials ----
panels <- bind_rows(
  rated %>% mutate(panel = "All trials"),
  rated %>% filter(click_correct) %>% mutate(panel = "Identification correct")) %>%
  mutate(panel = factor(panel, levels = c("All trials", "Identification correct")))
p_panels <- panels %>%
  group_by(panel, condition_label, item, length_level, n_words) %>%
  summarize(value = mean(acceptability), .groups = "drop") %>%
  ggplot(aes(n_words, value, color = condition_label, fill = condition_label)) +
  geom_point(alpha = 0.4, position = position_jitter(width = 0.3, height = 0, seed = 2026)) +
  geom_smooth(data = panels, aes(y = acceptability), method = "lm", formula = y ~ x) +
  facet_wrap(~ panel, nrow = 1) +
  scale_color_manual(name = "Condition", values = condition_colors,
                     aesthetics = c("colour", "fill")) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Acceptability", breaks = seq(1, 7, 2)) +
  coord_cartesian(ylim = c(1, 7))
p_panels
ggsave("figures/artificial_exchange_identification_ratings.png", p_panels, width = 11, height = 5.5)

## ---- C3: lmer, condition x length, click-correct trials only ----
for (f in lmer_formulas) {
  m_correct <- suppressWarnings(lmer(as.formula(f), control = lmer_ctrl,
                                     data = rated %>% filter(click_correct) %>% droplevels()))
  problems <- grep("singular", m_correct@optinfo$conv$lme4$messages,
                   ignore.case = TRUE, invert = TRUE, value = TRUE)
  if (length(problems) == 0) break
}
cat("\n---- acceptability ~ condition x length, click-correct trials ----\nrandom effects:", f, "\n")
print(round(summary(m_correct)$coefficients, 4))


## =============================================================================
## B3 / C4. The paper's three accounts: shape-constrained GAMs
## =============================================================================
## Fit to sentence means with by-item random intercepts and slopes, on three
## views of the rated trials, and compared on marginal log-likelihood within
## each view (scam's own logLik() is conditional on the fitted random effects,
## so they are integrated out directly, as in sections 4.2 and 4.4). The
## communication model is fitted with length and again with sentence surprisal
## in its place (same smooths, same random effects); binary and gradient
## grammaticality have no information term, so one fit of each serves under
## both, and the accounts are compared within each information measure.
datasets <- list(
  "All trials" = rated,
  "Identification correct" = rated %>% filter(click_correct),
  "Identification correct + flagged originals" = rated %>% filter(click_correct | condition_label == "Original"))

cat("\n\n===== Shape-constrained GAMs =====\n")
bind_rows(datasets, .id = "dataset") %>%
  group_by(dataset = factor(dataset, levels = names(datasets))) %>%
  summarize(trials = n(), originals = sum(condition_label == "Original"), .groups = "drop") %>%
  print()

likelihoods <- tibble()
for (dataset_name in names(datasets)) {
  sentence_means <- datasets[[dataset_name]] %>%
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
      likelihoods <- bind_rows(likelihoods, tibble(dataset = dataset_name, information = information,
                                                   model = model_name, converged = m$conv,
                                                   log_likelihood = marginal))
    }
  }
}
likelihoods <- likelihoods %>%
  mutate(dataset = factor(dataset, levels = names(datasets)),
         information = factor(information, levels = c("Length", "Surprisal")),
         model = factor(model, levels = names(model_colors))) %>%
  group_by(dataset, information) %>%
  mutate(likelihood_diff = log_likelihood - max(log_likelihood)) %>%
  ungroup()
likelihoods %>%
  arrange(dataset, information, desc(log_likelihood)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 1))) %>%
  as.data.frame() %>%
  print(right = FALSE)
write_csv(likelihoods, "outputs/artificial_exchange_identification_acceptability_likelihoods.csv")

## ---- Figure: likelihood differences, with length and with surprisal ----
p_likelihoods <- likelihoods %>%
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
        axis.text.y = element_text(size = 13))
p_likelihoods
ggsave("figures/artificial_exchange_identification_acceptability_likelihood_differences.png", p_likelihoods,
       width = 16, height = 5)
