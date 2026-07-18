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
