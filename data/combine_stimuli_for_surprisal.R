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

all_languages <- c("da", "de", "en", "fr","hi", "ko", "ru", "zh", "en_experimental")

materials <- read_csv("en/materials_en.csv") %>% 
  mutate(language = "en", item=as.character(item))
materials <- materials[0,]

for (current_language in all_languages) {
  ## Read materials files.
  ## NB: Mandarin deliberately uses the UNSEGMENTED materials_zh.csv here, not the
  ## jiebaR-segmented materials_zh_for_analysis.csv used elsewhere. The Goldfish
  ## Chinese model is a SentencePiece model trained on running (unsegmented) text,
  ## so feeding it the space-separated jieba output makes it score every inserted
  ## space as a separate token, inflating sum_surprisal for Mandarin relative to
  ## the other languages. Word segmentation is still what defines n_words and
  ## dl_distance -- only the language-model input is unsegmented.
  current_dat <- read_csv(
    file.path(
      current_language,
      paste0("materials_", current_language, ".csv")
    )
  )

  if (current_language == "en_experimental") {
    current_dat <- current_dat %>% 
      separate(condition, sep="_", into = c("condition", "type")) %>%
      mutate(
        type = substr(type, 1, 1),
        item = paste0(item, type)
      ) %>%
      select(-type)
  }
  
  materials <- materials %>%
    add_row(current_dat %>% mutate(language = current_language, item = as.character(item)))
  
}

materials <- materials %>%
  filter(experiment == "exp_b",
         grepl("original", condition)) %>%
  mutate(
    item = paste(language, item, sep = "_") 
  ) %>%
  select(
    language,
    item,
    sentence
  ) %>%
  distinct()
  
write_csv(materials, "materials_for_surprisal.csv")
