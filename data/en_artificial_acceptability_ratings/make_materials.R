## =============================================================================
## Build the Surveyor questionnaire for the acceptability-only version of the
## artificial length-controlled English materials.
##
## Same sentences as data/en_artificial_exchange_identification -- read from its validated
## stimuli file, so every artificial-materials study shows byte-identical
## sentences -- presented exactly as in the original English experiment
## (data/en): one acceptability rating per sentence, with the same prompt,
## the same 1-7 scale and the same instructions. There is no error question
## and no exchange identification task. The original experiment also asked a comprehension
## question after each sentence; this version has none. Writes materials_en_artificial_acceptability_ratings.csv (stimuli) and
## demographics_en_artificial_acceptability_ratings.csv (survey properties).
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

## Prompt and scale copied from the original English experiment's stimuli file.
original_en <- read_csv("../en/materials_en.csv", show_col_types = FALSE) %>%
  filter(experiment == "exp_b")
acceptability_prompt <- unique(original_en$prompt)
acceptability_options <- unique(original_en$options)
stopifnot(length(acceptability_prompt) == 1, length(acceptability_options) == 1)

## ---- Sentences ----
## 27 items x 9 conditions (original / swap-1 / swap-3 x short / medium / long,
## no shuffled condition), so each participant rates 27 sentences.
sentences <- read_csv("../en_artificial_exchange_identification/materials_en_artificial_exchange_identification.csv",
                      show_col_types = FALSE) %>%
  select(item, condition, sentence)

## ---- Acceptability rating ----
## Laid out as the original English exp_b rows: literal_response NA, no config.
materials <- sentences %>%
  transmute(
    experiment = "exp_b",
    item, condition, sentence,
    prompt = acceptability_prompt,
    question_type = "multiple_choice",
    options = acceptability_options,
    literal_response = "NA",
    config = ""
  ) %>%
  arrange(as.integer(item), condition)

write_csv(materials, "materials_en_artificial_acceptability_ratings.csv", na = "")

## ---- Survey properties (demographics) ----
## Consent and background questions as in the other artificial-materials studies.
## The instructions are the original English experiment's, minus its sentences
## about questions with correct answers and the 25% rejection rule: those were
## about the comprehension questions, which this version does not have.
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
  "Instructions: Please read each sentence and then answer the question immediately following.",
  "", ""
)

write_csv(demographics, "demographics_en_artificial_acceptability_ratings.csv", na = "")

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
cat("trials per participant:", n_distinct(materials$item), "sentences\n")
