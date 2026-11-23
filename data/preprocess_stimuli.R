## =============================================================================
## Data preparation, step 1 (with preprocess_en_experimental.R); then run
## concatenate_data_and_stimuli.R and get_slor_morcela_scores.py.
## Reads <language>/materials_<language>.csv for the 8 languages (Mandarin:
## materials_zh_for_analysis.csv) and writes <language>/stimuli_data.csv: one
## row per critical sentence (no comprehension questions, no nonsense
## sentences) with the normalized perceived and original sentence, the DL
## distance between them and the original's length in tokens (n_words).
## The length-controlled English experiment has its own script,
## preprocess_en_experimental.R, because its conditions also code length. The
## two scripts share the normalization and DL code line for line: a change to
## one must be copied to the other.
## =============================================================================

library(tidyverse)
library(stringdist)

## Note: the Mandarin materials were word-segmented with jiebaR ahead of time and
## are distributed pre-segmented as materials_zh_for_analysis.csv, which this
## script reads below. jiebaR is therefore not needed to run the pipeline; see
## "Provenance: Mandarin word segmentation" in DATA_DICTIONARY.md for the
## segmentation code that produced that file.

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

all_languages <- c("da", "de", "en", "fr","hi", "ko", "ru", "zh")

for (current_language in all_languages) {
  ## Read materials files
  if (current_language == "zh") {
    current_dat <- read_csv(
      file.path(
        current_language, 
        paste0("materials_", current_language, "_for_analysis.csv")
      )
    )
  } else {
    current_dat <- read_csv(
      file.path(
        current_language, 
        paste0("materials_", current_language, ".csv")
      )
    )
  }
  
  ## Remove comprehension questions
  current_dat <- current_dat %>%
    filter(experiment == "exp_b")
  
  ## Select relevant columns
  current_dat <- current_dat %>%
    select(item, condition, sentence)
  
  ## Remove nonsense conditions
  current_dat <- current_dat %>%
    filter(!grepl("non_sense", condition))
  
  ## Remove punctuation. Use the PCRE Unicode classes for punctuation (\p{P})
  ## and symbols (\p{S}) rather than [[:punct:]]: the latter also matches Unicode
  ## combining marks, which would strip the vowel signs and viramas out of the
  ## Devanagari (Hindi) materials and merge distinct words into one another.
  current_dat$sentence <- gsub("[\\p{P}\\p{S}]+", "", current_dat$sentence, perl = TRUE)
  
  ## Remove leading and trailing spaces
  current_dat$sentence <- gsub("^\\s+|\\s+$", "", current_dat$sentence)
  
  ## Collapse runs of whitespace into one space
  current_dat$sentence <- gsub("\\s+", " ", current_dat$sentence)
  
  ## Make everything lower case
  current_dat$sentence <- tolower(current_dat$sentence)
  
  ## Extract original sentences
  originals <- current_dat %>% 
    filter(condition == "original") %>%
    select(item, sentence) %>%
    rename(original_sentence = sentence)
  
  ## Join originals with full data frame
  current_dat <- current_dat %>%
    rename(perceived_sentence = sentence) %>%
    left_join(originals)
  
  ## Count swaps and other edits
  current_dat$dl_distance <- 0
  
  for (i in 1:nrow(current_dat)) {
    ## Create a superset of words in the perceived and intended sentences
    ## and then map each word to a unique character
    superset_of_words <- unique(c(
      unlist(strsplit(current_dat$perceived_sentence[i], " ")),
      unlist(strsplit(current_dat$original_sentence[i], " "))
    ))
    
    ## For each string (perceived and intended), map each word to a unique 
    ## character
    char_map <- setNames(
      c(letters, LETTERS)[1:length(superset_of_words)], 
      superset_of_words)
  
    char_string_perceived <- sapply(
      strsplit(current_dat$perceived_sentence[i], " ")[[1]], 
      function(word) char_map[[word]]
    ) %>%
      paste0(collapse = "")
    
    char_string_intended <- sapply(
      strsplit(current_dat$original_sentence[i], " ")[[1]], 
      function(word) char_map[[word]]
    ) %>%
      paste0(collapse = "")
    
    ## Compute Damerau-Levenshtein distance
    current_dat$dl_distance[i] <- stringdist(
      char_string_perceived, char_string_intended,
      method="dl",
      weight=c(1, 1, 1, 1)
    )
  }
  
  ## Get number of words from original_sentence
  current_dat$n_words <- sapply(
    strsplit(current_dat$original_sentence, " "), 
    length
  )
  
  ## Set language
  current_dat$language <- current_language
  
  ## Mark the sentences as UTF-8 and write with a byte-order mark
  ## (write_excel_csv), so that Excel shows the non-Latin scripts correctly
  Encoding(current_dat$original_sentence) <- "UTF-8"
  Encoding(current_dat$perceived_sentence) <- "UTF-8"
  
  write_excel_csv(current_dat, file.path(current_language, "stimuli_data.csv"))
  
}
  
