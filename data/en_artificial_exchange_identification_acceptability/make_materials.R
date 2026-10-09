## =============================================================================
## Build the Surveyor questionnaire for the exchange identification + acceptability version of
## the artificial length-controlled English materials.
##
## Same sentences as data/en_artificial_exchange_identification -- read from its validated
## stimuli file, so every artificial-materials study shows byte-identical
## sentences. On each trial the participant first clicks the words they think
## are misplaced (exactly the en_artificial_exchange_identification task), then rates the
## sentence's acceptability (exactly the original English prompt and scale).
##
## The click comes first because it has to: Surveyor shows the sentence only on
## the anchor trial, and a coordinated follow-up renders its prompt alone. The
## exchange identification task needs the sentence to click on, so it is the anchor and the rating
## is its coordinated follow-up -- the same arrangement as the error question
## and rating in data/en_artificial_binary_error_detection. Writes
## materials_en_artificial_exchange_identification_acceptability.csv (stimuli) and
## demographics_en_artificial_exchange_identification_acceptability.csv (survey properties).
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

## The click rows are the en_artificial_exchange_identification stimuli, unchanged; the rating
## prompt and scale are the original English experiment's.
## Read as plain text: its options, literal_response and config hold the literal
## string "NA", which read_csv would otherwise turn into missing values.
click_task <- read_csv("../en_artificial_exchange_identification/materials_en_artificial_exchange_identification.csv",
                       col_types = cols(.default = col_character()), na = character())
original_en <- read_csv("../en/materials_en.csv", show_col_types = FALSE) %>%
  filter(experiment == "exp_b")
acceptability_prompt <- unique(original_en$prompt)
acceptability_options <- unique(original_en$options)
stopifnot(length(acceptability_prompt) == 1, length(acceptability_options) == 1)

## ---- Exchange identification task ----
## The anchor trial: it carries the sentence and is what the participant sees
## first. Every column as in en_artificial_exchange_identification (per_word_multiple_choice,
## options, literal_response and config all "NA").
clicks <- click_task %>%
  select(experiment, item, condition, sentence, prompt, question_type,
         options, literal_response, config)

## ---- Acceptability rating ----
## Coordinated follow-up, shown immediately after its anchor. `sentence` is left
## empty on purpose: Surveyor does not reprint the sentence on a coordinated
## trial, it renders the prompt alone directly beneath the anchor.
acceptability <- click_task %>%
  transmute(
    experiment = "exp_b",
    item, condition,
    sentence = "",
    prompt = acceptability_prompt,
    question_type = "multiple_choice",
    options = acceptability_options,
    literal_response = "NA",
    config = paste0("coordinate;experiment%annot",
                    ";item%", item, ";condition%", condition)
  )

## Interleave so each sentence's two trials sit together, click first.
materials <- bind_rows(clicks, acceptability) %>%
  arrange(as.integer(item), condition, experiment)

write_csv(materials, "materials_en_artificial_exchange_identification_acceptability.csv", na = "")

## ---- Survey properties (demographics) ----
## As in en_artificial_exchange_identification, with the instructions extended to the rating
## that now follows each click trial.
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
         "at a time.<br><br>In some of the sentences, one or more pairs of words ",
         "have been exchanged, so they appear in the wrong order. For each ",
         "sentence, you will first click on every word that you think is in the ",
         "wrong place. Not every sentence has something wrong with it. If the ",
         "sentence looks fine to you, do not click any word.<br><br>",
         "You will then rate how acceptable the sentence sounds to you on a ",
         "scale from 1 (definitely not acceptable) to 7 (definitely acceptable). ",
         "Trust your intuition; there are no trick sentences, and you should not ",
         "overthink it.<br><br>Please read each sentence carefully and answer ",
         "every question"),
  "", ""
)

write_csv(demographics, "demographics_en_artificial_exchange_identification_acceptability.csv", na = "")

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
  filter(experiment == "annot") %>%
  transmute(key = paste(item, condition))
links <- materials %>%
  filter(str_starts(config, "coordinate")) %>%
  transmute(key = paste(item, condition))
cat("coordinate links resolving to a real anchor:",
    sum(links$key %in% anchors$key), "of", nrow(links), "\n")
cat("trials per participant:", n_distinct(materials$item), "sentences,",
    2 * n_distinct(materials$item), "screens\n")
