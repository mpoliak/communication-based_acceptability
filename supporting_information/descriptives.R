## =============================================================================
## Materials and participant descriptives (not an SI section)
##
## Descriptives of the materials and participants that no other script prints:
##   - sentence length in tokens, per experiment and pooled over all nine
##     experiments (Methods, Materials)
##   - comprehension accuracy by DL distance, per experiment and pooled (an
##     optional check that comprehension stays near ceiling despite errors)
##   - item x condition cells by DL distance and condition, written to
##     outputs/dl_distance_counts.txt, for the Methods (Materials), the Results
##     (cells excluded at DL >= 4), the Discussion (data with DL distance of 3
##     or less), the caption of the observed-vs-predicted figure (the conditions
##     behind the DL 0, 1 and 3 points; DL 2 cells not drawn) and SI Section 2
##     (shuffled cells kept at DL <= 3). Line references to this file depend on
##     its layout.
##   - each experiment's answers to the first-language question, and the number
##     of analysed participants whose first language is the experiment's
##     language, per experiment and in total (Methods, Participants)
## Only dl_distance_counts.txt is saved; the other tables are printed.
## =============================================================================

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

## Libraries, plotting theme, helper functions, `crit`, the by-item data frames
## and the shared plotting helpers all live in setup.R.
source("setup.R")

## ---- Materials Information ----
### ---- Sentence Length Across Languages ----
lengths <- crit %>%
  select(language, item, n_words) %>%
  distinct() %>%
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
      levels=c("Danish", "English", "French", "German", 
               "Hindi", "Korean", "Mandarin", "Russian", 
               "English\n(controlled)")
      )
  )


lengths %>%
  group_by(language) %>% 
  summarize(
    median_length = median(n_words),
    min_length = min(n_words),
    max_length = max(n_words),
    mean_length = mean(n_words),
    sd_length = sd(n_words),
  ) %>%
  print()

## Pooled over all nine experiments, over the same rows: the range and median
## of sentence length that the main text reports
lengths %>%
  summarize(
    language = "All nine experiments",
    median_length = median(n_words),
    min_length = min(n_words),
    max_length = max(n_words),
    mean_length = mean(n_words),
    sd_length = sd(n_words)
  ) %>%
  print()

### ---- Comprehension Accuracy by Edit Distance ----
## Descriptive check: how accurately did participants answer the comprehension
## question that followed each sentence, as a function of how corrupted that
## sentence was? Restricted to the three error conditions in the design (0, 1, and
## 3 adjacent-token exchanges); the rare intermediate DL distance of 2 is dropped
## as it is only sparsely populated. Scored the same way as the exclusion
## criterion in setup.R (correct = literal_response == response), on the same
## retained sample (participants below the 80% threshold already excluded).
comprehension_accuracy <- comprehension %>%
  filter(!participant_id %in% exclusions,
         dl_distance %in% c(0, 1, 3)) %>%
  mutate(correct = literal_response == response) %>%
  group_by(language, dl_distance) %>%
  summarize(n = n(), k = sum(correct, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    accuracy = k / n,
    ## Wilson 95% confidence interval for a binomial proportion (well behaved
    ## near the ceiling, where a normal-approximation interval would exceed 1)
    z = 1.96,
    wilson_center = (accuracy + z^2 / (2 * n)) / (1 + z^2 / n),
    wilson_half   = z * sqrt(accuracy * (1 - accuracy) / n + z^2 / (4 * n^2)) /
                    (1 + z^2 / n),
    lower = wilson_center - wilson_half,
    upper = wilson_center + wilson_half,
    language = label_language(language)
  )

## Accuracy per experiment and DL distance, with its Wilson 95% CI
comprehension_accuracy %>%
  select(language, dl_distance, n, k, accuracy, lower, upper) %>%
  print(n = Inf)

## Pooled accuracy per DL distance, printed for reference
comprehension %>%
  filter(!participant_id %in% exclusions, dl_distance %in% c(0, 1, 3)) %>%
  mutate(correct = literal_response == response) %>%
  group_by(dl_distance) %>%
  summarize(accuracy = mean(correct, na.rm = TRUE), n = n(), .groups = "drop") %>%
  print()

### ---- Count Damerau-Levenshtein distance for each language ----
distance_counts <- crit %>%
  select(language, item, condition, dl_distance, n_words) %>%
  distinct() %>% mutate(
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
      levels=c("Danish", "English", "French", "German", 
               "Hindi", "Korean", "Mandarin", "Russian", 
               "English\n(controlled)")
    )
  )

## Item x condition cells per experiment, over all DL distances
table(distance_counts$language) %>%
  print()

sink("outputs/dl_distance_counts.txt")
print("Original Sentences")
table(
  distance_counts$dl_distance[distance_counts$condition == "original"], 
  distance_counts$language[distance_counts$condition == "original"]
) %>%
  print(., n=10^4)

print("1 Exchange")
table(
  distance_counts$dl_distance[distance_counts$condition == "swap_1"], 
  distance_counts$language[distance_counts$condition == "swap_1"]
) %>%
  print(., n=10^4)

print("3 Exchanges")
table(
  distance_counts$dl_distance[distance_counts$condition == "swap_3"], 
  distance_counts$language[distance_counts$condition == "swap_3"]
) %>%
  print(., n=10^4)

print("Shuffle")
table(
  distance_counts$dl_distance[distance_counts$condition == "random"], 
  distance_counts$language[distance_counts$condition == "random"]
) %>%
  print(., n=10^4)

sink()

## ---- Demographic Information ----
## Each experiment's answers to its first-language question, analysed
## participants only. The question is first_language in da, de, hi, ko and ru,
## english_first_language in en, en_exp and zh, and french_first_language in
## fr. Most tables cross the answer with country of birth, which the paper does
## not report.
demographics <- dat %>%
  select(
    participant_id, country_birth, country_birth_text, first_language, 
    first_language_text,second_language , second_language_text, 
    language_contexts, english_first_language, country_of_birth, 
    french_first_language,  french_contexts, language
  ) %>%
  distinct

## Danish
danish_demographics <- demographics %>% filter(language == "da")
table(danish_demographics$first_language, danish_demographics$country_birth) %>%
  print()

## English
english_demographics <- demographics %>% filter(language == "en")

table(english_demographics$english_first_language, 
      english_demographics$country_of_birth) %>%
  print()

## English (Length-Controlled Experiment)
english_exp_demographics <- demographics %>% filter(language == "en_exp")

table(english_exp_demographics$english_first_language, 
      english_exp_demographics$country_of_birth,
      useNA = "always"
) %>%
  print()

## French
french_demographics <- demographics %>% filter(language == "fr")

table(
  french_demographics$french_first_language,
  french_demographics$country_of_birth
) %>%
  print()

## German
german_demographics <- demographics %>% filter(language == "de")

table(
  german_demographics$first_language,
  german_demographics$country_birth
) %>%
  print()

## Hindi
hindi_demographics <- demographics %>% filter(language == "hi")

table(
  hindi_demographics$first_language,
  hindi_demographics$country_birth
) %>%
  print()

table(hindi_demographics$first_language_text) %>%
  print()

## Korean
korean_demographics <- demographics %>% filter(language == "ko")

table(
  korean_demographics$first_language,
  korean_demographics$country_of_birth
) %>%
  print()

## Chinese
chinese_demographics <- demographics %>% filter(language == "zh")

table(chinese_demographics$english_first_language) %>%
  print()

## Russian
russian_demographics <- demographics %>% filter(language == "ru")

table(
  russian_demographics$first_language,
  russian_demographics$country_birth
) %>%
  print()

### ---- First Language Is the Experiment's Language ----
## Analysed participants whose first language is the experiment's language:
## first_language_is_target in full_data.csv (see ../data/DATA_DICTIONARY.md),
## which harmonizes the three first-language questions above. The analysed
## participants are those in `crit`; their rows are taken from `dat` with
## crit's filters, because crit does not keep first_language_is_target. One
## row per participant: the answer is the same on all of a participant's rows.
first_language_target <- dat %>%
  filter(experiment == "exp_b", !grepl("non_sense", condition)) %>%
  distinct(language, participant_id, first_language_is_target)

stopifnot(nrow(first_language_target) == n_distinct(crit$participant_id))

## Per experiment
first_language_target %>%
  mutate(language = label_language(language)) %>%
  group_by(language) %>%
  summarize(
    n_participants = n(),
    n_first_language_target = sum(first_language_is_target)
  ) %>%
  print()

## In total
first_language_target %>%
  summarize(
    n_participants = n(),
    n_first_language_target = sum(first_language_is_target)
  ) %>%
  print()
