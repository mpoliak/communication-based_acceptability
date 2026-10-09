## =============================================================================
## Build the Surveyor questionnaire for the "identify the misplaced words" task.
##
## Uses the ORIGINAL English materials (data/en), not the length-controlled ones,
## and only the original and swap_1 conditions -- so 72 items x 2 conditions.
## There is no acceptability rating and no Yes/No error question: every word in
## the sentence is clickable, and the participant clicks the ones they think are
## out of place. That is Surveyor's `per_word_multiple_choice` question type, so
## `options`, `literal_response` and `config` are all NA.
##
## `literal_response` is left NA because the answer is a set of word positions,
## which that column cannot express. The ground truth is recoverable at analysis
## time by diffing each swap_1 sentence against its own original: every swap_1
## sentence here displaces exactly 2 adjacent positions (verified below), and
## the correct response for an original sentence is to click nothing.
##
## Writes materials_en_exchange_identification.csv and demographics_en_exchange_identification.csv.
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

## ---- Sentences: original and swap_1 only ----
sentences <- read_csv("../en/materials_en.csv", show_col_types = FALSE) %>%
  filter(experiment == "exp_b", condition %in% c("original", "swap_1")) %>%
  select(item, condition, sentence)

materials <- sentences %>%
  transmute(
    experiment = "annot",
    item, condition,
    sentence,
    prompt = "Click on all the words that are misplaced",
    question_type = "per_word_multiple_choice",
    options = "NA",
    literal_response = "NA",
    config = "NA"
  ) %>%
  arrange(as.integer(item), condition)

write_csv(materials, "materials_en_exchange_identification.csv", na = "")

## ---- Survey properties (demographics) ----
## Mirrors the other questionnaires in this project. The instructions describe
## the click-the-words task and make clear that clicking nothing is a valid
## response, since half the items have nothing wrong with them.
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
         "at a time.<br><br>In some of the sentences, two words have been ",
         "exchanged, so they appear in the wrong order. Your task is to find ",
         "them: click on every word that you think is in the wrong place.<br><br>",
         "Not every sentence has something wrong with it. If the sentence looks ",
         "fine to you, move on without clicking any word.<br><br>Trust your ",
         "intuition; there are no trick sentences, and you should not overthink ",
         "it. Please read each sentence carefully"),
  "", ""
)

write_csv(demographics, "demographics_en_exchange_identification.csv", na = "")

## ---- Checks ----
## Surveyor rejects a stimuli file unless, within each experiment, the number of
## items is a multiple of the number of conditions and every (item, condition)
## pair appears exactly once.
cat("rows written:", nrow(materials), "\n")
materials %>%
  summarize(items = n_distinct(item), conditions = n_distinct(condition),
            latin_square_ok = n_distinct(item) %% n_distinct(condition) == 0,
            pairs_unique = n() == n_distinct(paste(item, condition))) %>%
  print()
cat("conditions:", paste(sort(unique(materials$condition)), collapse = ", "), "\n")
cat("question_type:", paste(unique(materials$question_type), collapse = ", "), "\n")

## Ground truth for scoring: how many positions each swap_1 sentence displaces.
displaced <- sentences %>%
  pivot_wider(names_from = condition, values_from = sentence) %>%
  rowwise() %>%
  mutate(n_displaced = sum(str_split(original, " ")[[1]] !=
                           str_split(swap_1, " ")[[1]])) %>%
  ungroup()
cat("displaced positions per swap_1 sentence:\n")
print(table(displaced$n_displaced))
