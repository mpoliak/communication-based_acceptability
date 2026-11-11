## =============================================================================
## SI Section 4: set-up shared by the three studies of constructed sentences
##
## Sourced by section4.1_acceptability_ratings_only.R,
## section4.2_binary_error_detection_acceptability_ratings.R,
## section4.3_exchange_identification_acceptability_ratings.R and
## section4_model_comparison_figure.R. It
##   - sources setup.R, for the libraries, the plotting theme, the main
##     analysis's marginal_logLik(), morcela_weights() and morcela_score(), the
##     five models' colours, shapes and legend labels, and the main text's
##     DL-distance colours and buckets (DL_COLORS, DL_BUCKET_LEVELS,
##     bucket_dl(), YLIM) for the acceptability figures (setup.R also reads the
##     main experiments' data, which section 4 does not use);
##   - sets the colours of the three conditions, for the correctness figures;
##   - builds the table of the 243 sentences that all three studies show.
##
## The sentences: 27 constructed (artificial) base items, each written in a
## short, medium and long version, and each version shown original, with 1 swap
## or with 3 swaps. That makes 9 conditions and 243 sentences, and every
## participant sees every base item once, in one of the 9 conditions. As in the
## controlled English experiment of the main text, every item comes in length
## versions of one base item, so the analyses group by the base item, with that
## experiment's by-item random effects. Surprisal then varies within item, which
## is what makes the by-item surprisal slopes identifiable.
## =============================================================================

## Libraries, plotting theme, helper functions, and the models' colours, shapes
## and labels
source("setup.R")
library(stringdist)

## ---- Colours ----
## The three conditions have DL distance 0, 1 and 3, and take the colours of
## those distances in the main text's figures
CONDITION_COLORS <- c("Original" = DL_COLORS[["0"]],
                      "1 swap"   = DL_COLORS[["1-2"]],
                      "3 swaps"  = DL_COLORS[["3"]])

## ---- The 243 sentences ----
## In each study's materials file, the sentence is on the row of the trial's
## first screen (a rating that follows another question has an empty sentence)
materials <- bind_rows(
  acceptability_ratings_only = read_csv(
    "../data/en_artificial_acceptability_ratings/materials_en_artificial_acceptability_ratings.csv",
    show_col_types = FALSE) %>%
    filter(experiment == "exp_b"),
  binary_error_detection_acceptability_ratings = read_csv(
    "../data/en_artificial_binary_error_detection/materials_en_artificial_binary_error_detection.csv",
    show_col_types = FALSE) %>%
    filter(experiment == "exp_b_comprehension"),
  exchange_identification_acceptability_ratings = read_csv(
    "../data/en_artificial_exchange_identification_acceptability/materials_en_artificial_exchange_identification_acceptability.csv",
    show_col_types = FALSE) %>%
    filter(experiment == "annot"),
  .id = "task"
) %>%
  transmute(task, item = as.integer(item), cond = condition, sentence)

## Every study shows each of the 243 sentences once, and all three show the
## same sentences
stopifnot(
  all(count(materials, task)$n == 27 * 9),
  nrow(distinct(materials, item, cond, sentence)) == 27 * 9
)

sentences <- materials %>%
  distinct(item, cond, sentence) %>%
  mutate(
    length_level = factor(str_extract(cond, "[^_]+$"),
                          levels = c("short", "medium", "long")),
    condition_label = factor(recode(str_extract(cond, "^[^_]+"),
                                    "original" = "Original",
                                    "swap-1" = "1 swap",
                                    "swap-3" = "3 swaps"),
                             levels = names(CONDITION_COLORS))
  )

## Damerau-Levenshtein distance from each version to its original, computed as
## in data/preprocess_stimuli.R: lowercase, drop punctuation, write each
## distinct word of the item (at that length) as one letter, and take
## stringdist(method = "dl") between the letter strings. Length is the number
## of space-separated words.
sentence_letters <- sentences %>%
  separate_rows(sentence, sep = " ") %>%
  mutate(word = str_remove_all(str_to_lower(sentence), "[^a-z0-9]")) %>%
  group_by(item, length_level) %>%
  mutate(letter = c(letters, LETTERS)[match(word, unique(word))]) %>%
  group_by(item, cond) %>%
  summarize(letter_string = paste(letter, collapse = ""), n_words = n(),
            .groups = "drop")

## The language-model scores, from the one source of the main analysis
## (data/slor_morcela_scores.csv: Goldfish 1000mb, log2, model-initial unigram)
artificial_scores <- read_csv("../data/slor_morcela_scores.csv",
                              col_types = cols(item = col_character()),
                              show_col_types = FALSE) %>%
  filter(language == "en_artificial", unigram_model == "model_initial") %>%
  transmute(item = as.integer(item), cond = condition, log_p_sentence, log_p_unigram)

## Information(s) in the communication account is the surprisal of the intended
## sentence, i.e. of the original at the same length: the original and swapped
## versions of an item at one length share it, and its short, medium and long
## versions differ. SLOR and MORCELA instead score the sentence that was
## actually read, so they take each condition's own log probabilities.
original_surprisal <- artificial_scores %>%
  filter(str_starts(cond, "original")) %>%
  transmute(item,
            length_level = factor(str_extract(cond, "[^_]+$"),
                                  levels = c("short", "medium", "long")),
            sum_surprisal = -log_p_sentence)

sentences <- sentences %>%
  left_join(sentence_letters, by = c("item", "cond")) %>%
  group_by(item, length_level) %>%
  mutate(original_string = letter_string[condition_label == "Original"]) %>%
  ungroup() %>%
  left_join(original_surprisal, by = c("item", "length_level")) %>%
  left_join(artificial_scores, by = c("item", "cond")) %>%
  mutate(
    dl_distance = stringdist(letter_string, original_string, method = "dl"),
    weighted_dl_distance = dl_distance / n_words,
    weighted_dl_distance_surprisal = dl_distance / sum_surprisal,
    grammatical_numeric = ifelse(dl_distance == 0, .5, -.5)
  ) %>%
  select(item, cond, sentence, length_level, condition_label, n_words,
         sum_surprisal, dl_distance, weighted_dl_distance,
         weighted_dl_distance_surprisal, grammatical_numeric,
         log_p_sentence, log_p_unigram)

## Every sentence found its scores, and every condition has the DL distance it
## was built to have: 0, 1 and 3. (So the main analysis's restriction to DL
## distances up to 3, which drops its shuffled condition, would remove nothing
## here, and is not repeated.)
stopifnot(
  nrow(sentences) == 27 * 9,
  !anyNA(sentences),
  all(sentences$dl_distance ==
        c("Original" = 0, "1 swap" = 1, "3 swaps" = 3)[as.character(sentences$condition_label)])
)

cat("DL distance by condition (number of sentences):\n")
sentences %>% count(condition_label, dl_distance) %>% print()

## The 27 original sentences of each length version: their length in words and
## their surprisal in bits (the Information(s) of all three conditions)
cat("Length (words) and surprisal (bits) of the original sentences, by length version:\n")
sentences %>%
  filter(condition_label == "Original") %>%
  group_by(length_level) %>%
  summarize(n_sentences = n(),
            min_words = min(n_words), median_words = median(n_words),
            max_words = max(n_words),
            min_surprisal = min(sum_surprisal), median_surprisal = median(sum_surprisal),
            max_surprisal = max(sum_surprisal)) %>%
  print(width = Inf)
