## =============================================================================
## SI Section 4.1: Finding the misplaced words in the artificial materials
##
## en_artificial_exchange_identification. 27 constructed items -- written for the experiment
## rather than taken from a corpus -- each in a short, medium and long version,
## and each of those original, with 1 swap or with 3 swaps: 9 conditions, every
## participant sees every item once. Every word is clickable and participants
## click the words they think are misplaced. Surveyor's per_word_multiple_choice
## writes one row per trial, "NONE" when nothing was clicked.
##
## The same analysis as the en_exchange_identification part of section4_binary_error_detection.R, on
## these data only. Length varies WITHIN item here, so a length effect on
## noticing cannot be put down to long sentences simply being different ones.
##
## Figures: figures/artificial_exchange_identification_{by_length,outcomes,
## outcomes_by_length,outcomes_original,outcomes_original_by_length,agreement,
## agreement_by_length}.png
## =============================================================================

library(tidyverse)
library(lme4)

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

## ---- Colours ----
## The same colours as the other section 4 figures. The outcome palette is built
## with setNames: a "±" written as c("..." = colour) is mangled in a non-UTF-8
## locale, and the Correct level then silently falls back to grey.
condition_colors <- c("Original" = "#EB7AA3", "1 swap" = "#1E88E5", "3 swaps" = "#FFC107")
measure_colors <- c("Noticed the swap" = "#1E88E5",
                    "Located it correctly" = "#004D40",
                    "False alarm (no error present)" = "#EB7AA3")
## "Correct, partly" only arises with 3 swaps: the error was located, with no
## stray clicks, but at least one swapped region went unmarked. It is a lighter
## green, and a dashed curve in the by-length figure, to read as a kind of
## correct.
outcome_colors <- setNames(
  c("#004D40", "#5FA894", "#FFC107", "#EB7AA3", "#BBBBBB"),
  c("Correct\n(target \u00b1 adjacent words)", "Correct, partly\n(some swaps not marked)",
    "Mixed\n(target + distant words)", "Wrong\n(other words only)",
    "Wrong\n(clicked nothing)"))
outcome_linetypes <- setNames(c("solid", "dashed", "solid", "solid", "solid"),
                              names(outcome_colors))
original_colors <- c("Correct\n(clicked nothing)" = "#004D40",
                     "False alarm\n(1+ words marked)" = "#EB7AA3")

## ---- Materials: one row per word of every sentence version ----
materials <- read_csv("../data/en_artificial_exchange_identification/materials_en_artificial_exchange_identification.csv",
                      show_col_types = FALSE) %>%
  transmute(item, cond = condition, sentence,
            length_level = factor(str_extract(condition, "[^_]+$"),
                                  levels = c("short", "medium", "long")),
            condition_label = factor(recode(str_extract(condition, "^[^_]+"),
                                            "original" = "Original",
                                            "swap-1" = "1 swap",
                                            "swap-3" = "3 swaps"),
                                     levels = names(condition_colors)))

words <- materials %>%
  separate_rows(sentence, sep = " ") %>%
  rename(word = sentence) %>%
  group_by(item, cond) %>%
  mutate(position = row_number(), n_words = n()) %>%
  ungroup()

## A word is a target (displaced) when it differs from the word at the same
## position in that item's original of the same length. Consecutive targets form
## a run: one run per swapped pair, or per chain where two swaps share a word.
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

cat("displaced words and separate displaced runs, per sentence version:\n")
words %>%
  group_by(condition_label, item, cond) %>%
  summarize(targets = sum(is_target), runs = n_distinct(target_run, na.rm = TRUE),
            .groups = "drop") %>%
  count(condition_label, targets, runs, name = "sentences") %>%
  print()

## ---- Responses ----
responses <- read_csv("../data/en_artificial_exchange_identification/en_artificial_exchange_identification_data.csv",
                      show_col_types = FALSE)

cat("\n---- export integrity ----\n")
cat("participants:", n_distinct(responses$participant_id), "\n")
responses %>%
  count(participant_id, name = "trials") %>%
  count(trials, name = "participants") %>%
  print()
responses %>% count(cond, name = "trials") %>% print()

## Surveyor balances conditions WITHIN each participant (3 of each), but which
## item lands in which condition is pseudo-random, so item x condition cells
## hold uneven numbers of responses. Pooled rates and the models are
## unaffected; the per-sentence summaries below rest on these few each.
cells <- responses %>% count(item, cond)
cat("item x condition cells seen:", nrow(cells), "of", nrow(materials),
    "| responses per cell: median", median(cells$n),
    "range", min(cells$n), "-", max(cells$n), "\n")

trials <- responses %>%
  mutate(trial = row_number()) %>%
  select(trial, participant_id, item, cond, response) %>%
  left_join(materials %>% select(item, cond, length_level, condition_label),
            by = c("item", "cond"))

## One row per clicked word. Words are percent-escaped for % , and : so they
## survive the position:word encoding, hence URLdecode.
clicks <- trials %>%
  filter(response != "NONE") %>%
  separate_rows(response, sep = ",") %>%
  separate(response, into = c("position", "clicked_word"), sep = ":", extra = "merge") %>%
  transmute(trial, item, cond, position = as.integer(position),
            clicked_word = map_chr(clicked_word, URLdecode))

## Every stored word should be the word actually at that position.
clicks %>%
  left_join(words %>% select(item, cond, position, word),
            by = c("item", "cond", "position")) %>%
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

## ---- Scoring each trial ----
## Correctly located: at least one displaced word clicked, and every other
## click in an unbroken run with the displaced words -- marking one or both
## words of a swapped pair, or those plus the words right beside them, all
## point at the error. Checked run by run: take the positions that were clicked
## or are targets, split them into runs of consecutive positions, and flag any
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
              summarize(n_words = n(),
                        n_clicked = sum(clicked),
                        n_on_target = sum(clicked & is_target),
                        runs_total = n_distinct(target_run, na.rm = TRUE),
                        runs_found = n_distinct(target_run[clicked & is_target], na.rm = TRUE),
                        .groups = "drop"),
            by = "trial") %>%
  left_join(stray_clicks, by = "trial") %>%
  mutate(flagged = n_clicked > 0,
         located = n_on_target > 0 & !any_stray,
         ## "all swaps" is counted in swapped regions: where two swaps share a
         ## word they form one region, so a 3-swap sentence has 1 to 3 of them
         outcome = case_when(
           !flagged         ~ "Wrong\n(clicked nothing)",
           located & runs_found < runs_total ~ "Correct, partly\n(some swaps not marked)",
           located          ~ "Correct\n(target \u00b1 adjacent words)",
           n_on_target == 0 ~ "Wrong\n(other words only)",
           TRUE             ~ "Mixed\n(target + distant words)"),
         outcome = factor(outcome, levels = names(outcome_colors)))

originals <- trials %>% filter(condition_label == "Original")
swapped <- trials %>% filter(condition_label != "Original")

## ---- (1) Original sentences: false alarms ----
cat("\n---- (1) ORIGINAL sentences: false alarms ----\n")
cat("false-alarm rate (clicked at least one word in a sentence with no error):",
    round(mean(originals$flagged), 3), "\n")
cat("mean words clicked when a false alarm happened:",
    round(mean(originals$n_clicked[originals$flagged]), 2), "\n")
originals %>%
  group_by(length_level) %>%
  summarize(trials = n(), false_alarm = mean(flagged)) %>%
  print()

cat("\noriginal sentences that drew a false alarm:\n")
originals %>%
  group_by(item, length_level) %>%
  summarize(trials = n(), false_alarms = sum(flagged), .groups = "drop") %>%
  filter(false_alarms > 0) %>%
  left_join(materials %>% filter(condition_label == "Original") %>%
              select(item, length_level, sentence),
            by = c("item", "length_level")) %>%
  mutate(sentence = str_trunc(sentence, 70)) %>%
  print(n = Inf)

## ---- (2) Swapped sentences: noticed and located ----
cat("\n---- (2) SWAPPED sentences: noticed, located correctly ----\n")
swapped %>%
  group_by(condition_label, length_level) %>%
  summarize(trials = n(),
            noticed = mean(flagged),
            located_correctly = mean(located),
            located_given_noticed = mean(located[flagged]),
            words_clicked_when_noticed = mean(n_clicked[flagged]),
            all_errors_touched = mean(runs_found == runs_total),
            .groups = "drop") %>%
  print(width = Inf)

cat("\nswapped sentences never located by anyone:\n")
swapped %>%
  group_by(condition_label, item, length_level) %>%
  summarize(trials = n(), noticed = sum(flagged), located = sum(located),
            .groups = "drop") %>%
  filter(located == 0) %>%
  print(n = Inf)

## ---- (3) Agreement about where the error is ----
## Mean pairwise Jaccard overlap between the sets of words that different
## participants clicked on the same sentence, among those who clicked anything:
## 1 = identical sets, 0 = no overlap. "collapsed" makes each displaced run one
## unit, so clicking either word of a swapped pair, or both, is the same
## answer; "raw" keeps every word separate.
##
## Jaccard's chance level depends on sentence length -- two people clicking at
## random agree far more often in a short sentence -- so "chance" repeats the
## collapsed calculation 200 times with random clicks: every participant keeps
## their number of clicks, but which words is random.
set.seed(2026)
random_clicks <- trial_words %>%
  semi_join(swapped %>% filter(flagged), by = "trial") %>%
  group_by(trial) %>%
  mutate(n_clicks = sum(clicked)) %>%
  ungroup() %>%
  expand_grid(rep = 1:200) %>%
  mutate(u = runif(n())) %>%
  group_by(rep, trial) %>%
  filter(rank(u) <= n_clicks) %>%
  ungroup()

click_units <- bind_rows(
  trial_words %>%
    filter(clicked) %>%
    transmute(measure = "collapsed", rep = 0L, trial, item, cond,
              unit = ifelse(is_target, paste0("run", target_run), as.character(position))),
  trial_words %>%
    filter(clicked) %>%
    transmute(measure = "raw", rep = 0L, trial, item, cond,
              unit = as.character(position)),
  random_clicks %>%
    transmute(measure = "chance", rep, trial, item, cond,
              unit = ifelse(is_target, paste0("run", target_run), as.character(position)))) %>%
  distinct()

## For every pair of participants on the same sentence: the size of each set
## and how many units they share. Jaccard = shared / (size_a + size_b - shared).
set_sizes <- click_units %>%
  count(measure, rep, item, cond, trial, name = "size")

shared_units <- click_units %>%
  inner_join(click_units, by = c("measure", "rep", "item", "cond", "unit"),
             suffix = c("_a", "_b"), relationship = "many-to-many") %>%
  filter(trial_a < trial_b) %>%
  count(measure, rep, item, cond, trial_a, trial_b, name = "shared")

pairwise <- set_sizes %>%
  inner_join(set_sizes, by = c("measure", "rep", "item", "cond"),
             suffix = c("_a", "_b"), relationship = "many-to-many") %>%
  filter(trial_a < trial_b) %>%
  left_join(shared_units, by = c("measure", "rep", "item", "cond", "trial_a", "trial_b")) %>%
  mutate(shared = replace_na(shared, 0L),
         jaccard = shared / (size_a + size_b - shared))

## One row per sentence version with at least two participants who clicked
## something. Chance-corrected agreement is (observed - chance) / (1 - chance):
## 0 = chance, 1 = perfect. Undefined where random clicking already agrees
## perfectly (chance = 1).
sentence_agreement <- pairwise %>%
  group_by(measure, rep, item, cond) %>%
  summarize(agreement = mean(jaccard), .groups = "drop") %>%
  group_by(measure, item, cond) %>%
  summarize(agreement = mean(agreement), .groups = "drop") %>%
  pivot_wider(names_from = measure, values_from = agreement) %>%
  left_join(materials %>% select(item, cond, length_level, condition_label),
            by = c("item", "cond")) %>%
  left_join(words %>% distinct(item, cond, n_words), by = c("item", "cond")) %>%
  left_join(trials %>% group_by(item, cond) %>%
              summarize(n_flagged = sum(flagged), .groups = "drop"),
            by = c("item", "cond")) %>%
  mutate(corrected = ifelse(chance < 1, (collapsed - chance) / (1 - chance), NA))

cat("\n---- (3) agreement about which words are misplaced ----\n")
sentence_agreement %>%
  group_by(condition_label) %>%
  summarize(sentences = n(), collapsed = mean(collapsed), raw = mean(raw),
            .groups = "drop") %>%
  print()

cat("\nagreement vs length, swapped sentences:\n")
sentence_agreement %>%
  filter(condition_label != "Original") %>%
  group_by(condition_label, length_level) %>%
  summarize(sentences = n(), observed = mean(collapsed), chance = mean(chance),
            corrected = mean(corrected, na.rm = TRUE), .groups = "drop") %>%
  print()

## ---- (4) The length effect, within item ----
## Is a swap harder to notice in a longer version of the SAME sentence? Length
## varies within item, so it takes a by-item random slope. 1-swap trials only:
## 3 swaps are noticed on 99% of trials, and a model with a 3-swap term is
## quasi-separated -- it returns absurdly small standard errors while
## reporting convergence.
model_data <- trials %>%
  mutate(n_words_z = as.numeric(scale(n_words)),
         noticed = as.integer(flagged),
         located = as.integer(located))
cat("\nn_words_z is per SD of length (SD =", round(sd(model_data$n_words), 1),
    "words). Estimates are log-odds.\n")

m_noticed <- glmer(noticed ~ n_words_z + (n_words_z | participant_id) + (n_words_z | item),
                   family = binomial, data = model_data %>% filter(condition_label == "1 swap"),
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
cat("\n---- noticed the swap ~ length (1-swap trials) ----\n")
print(round(summary(m_noticed)$coefficients, 4))

m_located <- glmer(located ~ n_words_z + (n_words_z | participant_id) + (n_words_z | item),
                   family = binomial, data = model_data %>% filter(condition_label == "1 swap"),
                   control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
cat("\n---- located it correctly ~ length (1-swap trials) ----\n")
print(round(summary(m_located)$coefficients, 4))

## Random intercepts only: false alarms are rare (30 in 900 trials). The
## maximal structure does not converge, and the uncorrelated-slopes version
## reports convergence but is degenerate (intercept -72 log-odds, slope +11.6
## for a 3% event rate).
m_false_alarm <- glmer(noticed ~ n_words_z + (1 | participant_id) + (1 | item),
                       family = binomial, data = model_data %>% filter(condition_label == "Original"),
                       control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
cat("\n---- false alarm ~ length (original trials) ----\n")
print(round(summary(m_false_alarm)$coefficients, 4))

## ---- Figure 1: noticing and locating the swap, by length ----
## Points are sentence means (item x length); lines are linear fits (lm) to the
## trial-level 0/1 data, as in every figure in this script. The false-alarm line
## comes from the original sentences and is drawn in both panels as the
## reference.
by_length <- bind_rows(
  swapped %>% transmute(panel = as.character(condition_label), item, length_level, n_words,
                        measure = "Noticed the swap", y = as.integer(flagged)),
  swapped %>% transmute(panel = as.character(condition_label), item, length_level, n_words,
                        measure = "Located it correctly", y = as.integer(located)),
  originals %>% transmute(panel = "1 swap", item, length_level, n_words,
                          measure = "False alarm (no error present)", y = as.integer(flagged)),
  originals %>% transmute(panel = "3 swaps", item, length_level, n_words,
                          measure = "False alarm (no error present)", y = as.integer(flagged))) %>%
  mutate(measure = factor(measure, levels = names(measure_colors)))

p_by_length <- ggplot(mapping = aes(x = n_words, color = measure, fill = measure)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey80") +
  geom_point(data = by_length %>%
               group_by(panel, measure, item, length_level, n_words) %>%
               summarize(value = mean(y), .groups = "drop"),
             aes(y = value), alpha = 0.3) +
  geom_smooth(data = by_length, aes(y = y), method = "lm", formula = y ~ x) +
  facet_wrap(~ panel) +
  scale_color_manual(name = "", values = measure_colors, aesthetics = c("colour", "fill")) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Proportion of trials", labels = scales::percent,
                     breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1)) +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
  theme(legend.position = "top", legend.text = element_text(size = 13))
p_by_length
ggsave("figures/artificial_exchange_identification_by_length.png", p_by_length, width = 12, height = 5.5)

## ---- Figure 2: what got clicked on a swapped trial ----
p_outcomes <- swapped %>%
  count(condition_label, outcome) %>%
  group_by(condition_label) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  ggplot(aes(outcome, prop, fill = outcome)) +
  geom_col(alpha = 0.9, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", 100 * prop, n)), vjust = -0.3, size = 4) +
  ## free x, so the 1-swap panel has no empty "Correct, partly" slot
  facet_grid(~ condition_label, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = outcome_colors) +
  scale_y_continuous(name = "Proportion of swapped trials", labels = scales::percent,
                     expand = expansion(mult = c(0, 0.18))) +
  ## the long labels wrap once more to fit the bars
  scale_x_discrete(name = "", labels = setNames(c("Correct\n(target \u00b1\nadjacent words)",
                                                  "Correct, partly\n(some swaps\nnot marked)",
                                                  "Mixed\n(target +\ndistant words)",
                                                  "Wrong\n(other words only)",
                                                  "Wrong\n(clicked nothing)"),
                                                names(outcome_colors))) +
  theme(axis.text.x = element_text(size = 11))
p_outcomes
ggsave("figures/artificial_exchange_identification_outcomes.png", p_outcomes, width = 13, height = 5.5)

## ---- Figure 3: the same four outcomes by length ----
## Each outcome is its own 0/1 indicator with its own linear fit, so the lines
## sum to 1 at any length.
##
## Each point is one sentence -- one item at one length -- and its height is
## the share of that sentence's trials with that outcome. Several items share a
## length and each sentence has only a few trials, so points pile up on the same
## spot (above all at 0% and 100%). They are jittered -- 0.3 tokens sideways,
## 2 points of percentage up or down, fixed seed -- so every one stays visible.
## The lines use the unjittered trial-level data.
##
## "Correct, partly" cannot happen with a single swap, so it is left out of
## the 1-swap panel rather than drawn as a row of zeros.
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
ggsave("figures/artificial_exchange_identification_outcomes_by_length.png", p_outcomes_by_length,
       width = 13, height = 5.5)

## ---- Figures 4 and 5: original sentences ----
## No word is displaced, so the only correct response is to click nothing.
original_outcomes <- originals %>%
  mutate(outcome = factor(ifelse(flagged, "False alarm\n(1+ words marked)",
                                 "Correct\n(clicked nothing)"),
                          levels = names(original_colors)))

p_outcomes_original <- original_outcomes %>%
  count(outcome) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(outcome, prop, fill = outcome)) +
  geom_col(alpha = 0.9, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", 100 * prop, n)), vjust = -0.3, size = 4.5) +
  scale_fill_manual(values = original_colors) +
  scale_y_continuous(name = "Proportion of original trials", labels = scales::percent,
                     expand = expansion(mult = c(0, 0.18))) +
  scale_x_discrete(name = "") +
  theme(axis.text.x = element_text(size = 12))
p_outcomes_original
ggsave("figures/artificial_exchange_identification_outcomes_original.png", p_outcomes_original,
       width = 10, height = 5.5)

## Points are sentences (one item at one length), jittered as in Figure 3.
original_long <- original_outcomes %>%
  select(item, length_level, n_words, outcome) %>%
  expand_grid(level = names(original_colors)) %>%
  mutate(y = as.integer(as.character(outcome) == level),
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
ggsave("figures/artificial_exchange_identification_outcomes_original_by_length.png", p_outcomes_original_by_length,
       width = 11, height = 5.5)

## ---- Figure 6: agreement about where the error is ----
## One value per sentence version with at least two participants who clicked
## something. Few original sentences draw two false alarms, so the Original
## panel holds only a handful of sentences.
p_agreement <- ggplot(sentence_agreement, aes(collapsed, fill = condition_label)) +
  geom_histogram(binwidth = 0.1, boundary = 0, color = "white", alpha = 0.9) +
  facet_wrap(~ condition_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(name = "Condition", values = condition_colors) +
  ## coord_cartesian, not scale limits: a hard limit drops the bin at 1.0
  scale_x_continuous(name = "Agreement about which words are misplaced\n(mean pairwise Jaccard)",
                     breaks = seq(0, 1, 0.25)) +
  coord_cartesian(xlim = c(-0.05, 1.05)) +
  scale_y_continuous(name = "Sentences") +
  theme(legend.position = "none")
p_agreement
ggsave("figures/artificial_exchange_identification_agreement.png", p_agreement, width = 12, height = 5)

## ---- Figure 7: agreement against length ----
## Point size and line weight = participants who clicked something. Raw
## agreement only; the chance-corrected values are printed in (3). For 3 swaps
## the two diverge: raw agreement falls with length, but chance-corrected
## agreement does not, because short 3-swap sentences are mostly displaced
## and random clicks agree there anyway.
p_agreement_by_length <- sentence_agreement %>%
  filter(condition_label != "Original") %>%
  ggplot(aes(n_words, collapsed, color = condition_label, fill = condition_label)) +
  geom_point(aes(size = n_flagged), alpha = 0.45) +
  geom_smooth(aes(weight = n_flagged), method = "lm", formula = y ~ x) +
  scale_color_manual(name = "", values = condition_colors, aesthetics = c("colour", "fill")) +
  scale_size_continuous(name = "participants who\nmarked something", range = c(1, 6)) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Agreement about which words are misplaced\n(mean pairwise Jaccard)",
                     breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("Swapped sentences, trials where an error was marked") +
  theme(legend.position = "right", plot.title = element_text(size = 15))
p_agreement_by_length
ggsave("figures/artificial_exchange_identification_agreement_by_length.png", p_agreement_by_length,
       width = 10, height = 5.5)
