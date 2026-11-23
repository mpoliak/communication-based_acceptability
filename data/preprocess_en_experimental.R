## =============================================================================
## Data preparation, step 1 for the length-controlled English experiment
## (preprocess_stimuli.R does the 8 languages); then run
## concatenate_data_and_stimuli.R and get_slor_morcela_scores.py.
## Reads en_experimental/materials_en_experimental.csv and writes
## en_experimental/stimuli_data.csv, with the same columns as the other
## languages' files. It has its own script because its conditions also code
## length (e.g. swap-1_short), which goes into the item id instead. The
## normalization, the join to the originals and the DL code repeat
## preprocess_stimuli.R line for line: a change to one must be copied to the
## other.
## =============================================================================

library(tidyverse)
library(stringdist)

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

en_exp <- read_csv("en_experimental/materials_en_experimental.csv")

## Remove comprehension questions
en_exp <- en_exp %>%
  filter(experiment == "exp_b")

## Select relevant columns
en_exp <- en_exp %>%
  select(item, condition, sentence)

## Remove punctuation (see preprocess_stimuli.R: use the PCRE Unicode punctuation
## and symbol classes rather than [[:punct:]], which also strips combining marks)
en_exp$sentence <- gsub("[\\p{P}\\p{S}]+", "", en_exp$sentence, perl = TRUE)

## Remove leading and trailing spaces
en_exp$sentence <- gsub("^\\s+|\\s+$", "", en_exp$sentence)

## Collapse runs of whitespace into one space
en_exp$sentence <- gsub("\\s+", " ", en_exp$sentence)

## Make everything lower case
en_exp$sentence <- tolower(en_exp$sentence)

## Separate into length conditions
en_exp <- en_exp %>%
  separate(
    col="condition",
    into=c("condition", "length"),
    sep="_"
  ) %>%
  mutate(
    condition = gsub("-", "_", condition)
  )

## Item ids: the base item plus s, m or l for its short, medium or long version
## (e.g. 62m)
en_exp <- en_exp %>%
  mutate(
    item = paste0(item, substr(length, 1, 1))
  ) %>%
  select(-length)

## Extract original sentences
originals <- en_exp %>% 
  filter(condition == "original") %>%
  select(item, sentence) %>%
  rename(original_sentence = sentence)

## Join originals with full data frame
en_exp <- en_exp %>%
  rename(perceived_sentence = sentence) %>%
  left_join(originals)

## Count swaps and other edits
en_exp$dl_distance <- 0

for (i in 1:nrow(en_exp)) {
  ## Create a superset of words in the perceived and intended sentences
  ## and then map each word to a unique character
  superset_of_words <- unique(c(
    unlist(strsplit(en_exp$perceived_sentence[i], " ")),
    unlist(strsplit(en_exp$original_sentence[i], " "))
  ))
  
  ## For each string (perceived and intended), map each word to a unique 
  ## character
  char_map <- setNames(
    c(letters, LETTERS)[1:length(superset_of_words)], 
    superset_of_words)
  
  char_string_perceived <- sapply(
    strsplit(en_exp$perceived_sentence[i], " ")[[1]], 
    function(word) char_map[[word]]
  ) %>%
    paste0(collapse = "")
  
  char_string_intended <- sapply(
    strsplit(en_exp$original_sentence[i], " ")[[1]], 
    function(word) char_map[[word]]
  ) %>%
    paste0(collapse = "")
  
  ## Compute Damerau-Levenshtein distance
  en_exp$dl_distance[i] <- stringdist(
    char_string_perceived, char_string_intended,
    method="dl",
    weight=c(1, 1, 1, 1)
  )
}

## Get number of words from original_sentence
en_exp$n_words <- sapply(
  strsplit(en_exp$original_sentence, " "), 
  length
)

## Set language
en_exp$language <- "en_exp"

## Mark the sentences as UTF-8 and write with a byte-order mark
## (write_excel_csv), so that Excel reads the file as UTF-8, as for the other
## languages' files
Encoding(en_exp$original_sentence) <- "UTF-8"
Encoding(en_exp$perceived_sentence) <- "UTF-8"

write_excel_csv(en_exp, file.path("en_experimental", "stimuli_data.csv"))
