## =============================================================================
## Build the Surveyor questionnaire for the binary error detection version of the
## artificial length-controlled English materials.
##
## Same sentences as data/en_artificial_exchange_identification -- read from its validated
## stimuli file, so the two tasks show byte-identical sentences -- presented as
## in data/en_binary_error_detection and data/en_experimental_binary_error_detection: an
## error-detection question, then an acceptability rating coordinated to follow
## it. Writes materials_en_artificial_binary_error_detection.csv (stimuli) and
## demographics_en_artificial_binary_error_detection.csv (survey properties).
## =============================================================================

library(tidyverse)

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

acceptability_options <-
  "1 - Definitely not acceptable_2_3_4_5_6_7 - Definitely acceptable"

## ---- Sentences ----
## 27 items x 9 conditions (original / swap-1 / swap-3 x short / medium / long,
## no shuffled condition), named as in en_experimental, so each participant
## rates 27 sentences.
sentences <- read_csv("../en_artificial_exchange_identification/materials_en_artificial_exchange_identification.csv",
                      show_col_types = FALSE) %>%
  select(item, condition, sentence)

## ---- Error-detection question ----
## This is the anchor trial: it carries the sentence and is what the participant
## sees first. literal_response is the correct answer -- "No" only when nothing
## was moved, "Yes" for every swap condition.
error_question <- sentences %>%
  transmute(
    experiment = "exp_b_comprehension",
    item, condition, sentence,
    prompt = "Are there any errors in this sentence?",
    question_type = "multiple_choice",
    options = "Yes_No",
    literal_response = ifelse(str_starts(condition, "original"), "No", "Yes"),
    config = ""
  )

## ---- Acceptability rating ----
## Coordinated follow-up, shown immediately after its anchor. `sentence` is left
## empty on purpose: Surveyor does not reprint the sentence on a coordinated
## trial, it renders the prompt alone directly beneath the anchor.
acceptability <- sentences %>%
  transmute(
    experiment = "exp_b",
    item, condition,
    sentence = "",
    prompt = "Please rate how acceptable the sentence is",
    question_type = "multiple_choice",
    options = acceptability_options,
    literal_response = "NA",
    config = paste0("coordinate;experiment%exp_b_comprehension",
                    ";item%", item, ";condition%", condition)
  )

## Interleave so each sentence's two trials sit together, as in en_binary_error_detection
materials <- bind_rows(error_question, acceptability) %>%
  arrange(as.integer(item), condition, desc(experiment))

write_csv(materials, "materials_en_artificial_binary_error_detection.csv", na = "")

## ---- Survey properties (demographics) ----
## Identical to data/en_experimental_binary_error_detection. Two thirds of the
## sentences contain a swap, so "Most of these materials have errors" still
## holds. Mirrors data/en_binary_error_detection/demographics_en_binary_error_detection.csv. The
## instructions describe the order the participant actually meets the two
## questions in: the error question comes first, the rating second.
demographics <- tribble(
  ~name, ~type, ~content, ~options, ~config,
  "consent", "text",
  paste0("<b>Consent Statement</b>:<br>By answering the following questions, ",
         "you are participating in a study being performed by cognitive ",
         "scientists in the MIT Department of Brain and Cognitive Sciences. If ",
         "you have questions about this research, please contact us. Your ",
         "participation in this research is voluntary. You may decline to ",
         "answer any or all of the following questions. You may decline ",
         "further participation, at any time, without adverse consequences. ",
         "Your anonymity is assured; the researchers who have requested your ",
         "participation will not receive any personal information about you."),
  "", "",

  "participant_ID", "question_text", "Please provide your Prolific ID", "", "",
  "native_language", "question_choices", "Is English your first language?",
  "Yes_No", "",
  "other_languages", "question_text",
  paste0("If you speak any languages other than English, please list them here ",
         "(write \"none\" if not)."), "", "",
  "country_of_birth", "question_choices", "In which country were you born?",
  "United States_United Kingdom_Canada_Australia_Ireland_New Zealand_Other", "",

  "instructions", "text",
  paste0("Instructions: In this study you will read a series of sentences, one ",
         "at a time.<br><br>For each sentence, you will first be asked whether ",
         "it contains any errors. Answer Yes or No. ",
         "<b>Most of these materials have errors.</b><br><br>",
         "You will then rate how acceptable the sentence sounds to you on a ",
         "scale from 1 (definitely not acceptable) to 7 (definitely acceptable). ",
         "Trust your intuition; there are no right or wrong answers, and you ",
         "should not overthink it.<br><br>Please read each sentence carefully ",
         "and answer every question"),
  "", ""
)

write_csv(demographics, "demographics_en_artificial_binary_error_detection.csv", na = "")

## ---- Checks ----
## Surveyor rejects a stimuli file unless, within each experiment, the number of
## items is a multiple of the number of conditions and every (item, condition)
## pair appears exactly once.
cat("rows written:", nrow(materials), "\n")
materials %>%
  group_by(experiment) %>%
  summarize(rows = n(), items = n_distinct(item), conditions = n_distinct(condition),
            latin_square_ok = n_distinct(item) %% n_distinct(condition) == 0,
            pairs_unique = n() == n_distinct(paste(item, condition)),
            .groups = "drop") %>%
  print()

anchors <- materials %>%
  filter(experiment == "exp_b_comprehension") %>%
  transmute(key = paste(experiment, item, condition))
links <- materials %>%
  filter(str_starts(config, "coordinate")) %>%
  transmute(key = paste0("exp_b_comprehension ", item, " ", condition))
cat("coordinate links resolving to a real anchor:",
    sum(links$key %in% anchors$key), "of", nrow(links), "\n")
cat("trials per participant:", n_distinct(materials$item), "sentences,",
    2 * n_distinct(materials$item), "screens\n")
