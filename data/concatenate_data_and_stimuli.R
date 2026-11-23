## =============================================================================
## Data preparation, step 2: run after preprocess_stimuli.R and
## preprocess_en_experimental.R. For each of the nine experiments, reads
## <experiment>/experiment_data.csv (the responses) and
## <experiment>/stimuli_data.csv (from step 1), joins them on item and
## condition, prefixes the item ids with the experiment (e.g. en_12), and
## stacks all experiments into full_data.csv, with the harmonized
## first_language_is_target and country_of_birth columns.
## get_slor_morcela_scores.py reads only the stimuli_data.csv files, so it does
## not depend on this step.
## =============================================================================

library(tidyverse)

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

all_languages <- c("da", "de", "en", "fr","hi", "ko", "ru", "zh", "en_experimental")

for (current_language in all_languages) {
  ## Read experiment data
  experiment_data <- read_csv(
    file.path(
      current_language, 
      "experiment_data.csv"
      )
    )
  
  ## Remove duplicated observations. A few sessions were submitted twice: the
  ## duplicate rows agree on every field (including the response) but their
  ## timestamps can differ by a second, so de-duplicating on the full row would
  ## miss them. De-duplicate ignoring the timestamp and keep the first
  ## submission. This affects participants hi_4 (144 rows, exact duplicates)
  ## and ko_34 (143 rows, timestamps one second apart).
  dedup_columns <- setdiff(names(experiment_data), "timestamp")
  experiment_data <- experiment_data %>%
    distinct(across(all_of(dedup_columns)), .keep_all = TRUE)

  ## The length-controlled English experiment stores its timestamp as text in
  ## another format (YYYY/MM/DD HH:MM:SS), which bind_rows() cannot stack with
  ## the other experiments' date-times, so it is dropped. The Mandarin timestamp
  ## is dropped too, although it is a normal date-time; Mandarin ids are already
  ## text, so as.character() is only a safeguard.
  if (current_language == "zh") {
    experiment_data <- experiment_data %>%
      select(-timestamp) %>%
      mutate(participant_id = as.character(participant_id))
  } else if (current_language == "en_experimental") {
    experiment_data <- experiment_data %>%
      select(-timestamp)
  }
  
  ## Read stimuli data
  stimuli_data <- read_csv(
    file.path(
      current_language, 
      "stimuli_data.csv"
    )
  )
  
  language_data <- experiment_data %>%
    rename(condition = cond) %>%
    left_join(stimuli_data,
              by=c("item", "condition"))
  
  language_data <- language_data %>%
    mutate(
      item = as.character(item),
      item = paste(current_language, item, sep="_")
    )
  
  ## Concatenate data
  if (current_language == all_languages[1]) {
    all_data <- language_data
  } else {
    all_data <- bind_rows(all_data, language_data)
  }
}

## Mark the sentences as UTF-8 (written with a byte-order mark below)
Encoding(all_data$original_sentence) <- "UTF-8"
Encoding(all_data$perceived_sentence) <- "UTF-8"

## Harmonize demographic columns that were named/asked differently across
## languages so that downstream code and external users have a single canonical
## column. The raw per-language columns are retained for provenance.
all_data <- all_data %>%
  mutate(
    ## Whether the participant reported the study's language as their first
    ## language. This was asked yes/no in en/en_exp/zh (english_first_language,
    ## incl. "是") and fr (french_first_language, "Oui"), and as a localized
    ## first-language dropdown in da/de/hi/ko/ru (first_language holds the
    ## language name). %in% is used so missing values resolve to FALSE, not NA.
    first_language_is_target =
      english_first_language %in% c("Yes", "是") |
      french_first_language  %in% "Oui"             |
      first_language %in% c("Dansk", "Deutsch",
                            "हिंदी",
                            "한국어",
                            "Русский"),
    ## Country of birth was stored as `country_of_birth` (en/en_exp/fr/ko/zh) and
    ## as `country_birth` (da/de/hi/ru); merge into the canonical column.
    country_of_birth = coalesce(country_of_birth, country_birth)
  )

## Write with a byte-order mark (write_excel_csv), so that Excel shows the
## non-Latin scripts correctly
write_excel_csv(all_data, "full_data.csv")
  
