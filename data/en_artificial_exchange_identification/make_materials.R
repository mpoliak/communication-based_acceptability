## =============================================================================
## Build the Surveyor questionnaire for the "identify the misplaced words" task
## on the ARTIFICIAL length-controlled English materials.
##
## Same task as en_exchange_identification: every word in the sentence is clickable, and the
## participant clicks the ones they think are out of place. That is Surveyor's
## `per_word_multiple_choice` question type, so `options`, `literal_response`
## and `config` are all NA, exactly as in en_exchange_identification.
##
## The design differs from en_exchange_identification in one way: length is manipulated within
## item. Each item exists in a short, medium and long version, each of which is
## original, 1-exchange or 3-exchange -- 9 conditions, no shuffled condition.
## Conditions are named version_length the way en_experimental names them
## (original_short, swap-1_medium, swap-3_long, ...), so they split cleanly on
## "_" at analysis time.
##
## `literal_response` is left NA because the answer is a set of word positions,
## which that column cannot express. The ground truth lives in
## fillers_exchanges.csv (`swap_slots`), which this script leaves untouched and
## which joins back to the data on item + condition.
##
## Surveyor requires the item count to be a multiple of the condition count, so
## this needs a multiple of 9 items. The script stops, and writes no stimuli
## file, until that holds.
##
## Writes materials_en_artificial_exchange_identification.csv and
## demographics_en_artificial_exchange_identification.csv.
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

## ---- Survey properties (demographics) ----
## Identical to en_exchange_identification except for one phrase in the instructions: "two
## words have been exchanged" becomes "one or more pairs of words have been
## exchanged", because the 3-exchange sentences contain three swapped pairs.
## Written first because it does not depend on the item count.
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
         "have been exchanged, so they appear in the wrong order. Your task is to ",
         "find them: click on every word that you think is in the wrong place.<br><br>",
         "Not every sentence has something wrong with it. If the sentence looks ",
         "fine to you, move on without clicking any word.<br><br>Trust your ",
         "intuition; there are no trick sentences, and you should not overthink ",
         "it. Please read each sentence carefully"),
  "", ""
)

write_csv(demographics, "demographics_en_artificial_exchange_identification.csv", na = "")

## ---- Sentences ----
## trim_ws = FALSE so the whitespace check below sees the file as written
sentences <- read_csv("fillers_exchanges.csv", show_col_types = FALSE, trim_ws = FALSE,
                      col_types = cols(swap_slots = col_character())) %>%
  mutate(version = recode(version, "original" = "original",
                          "1-exchange" = "swap-1", "3-exchange" = "swap-3"),
         condition = paste(version, length, sep = "_"))

## ---- Checks on the source materials ----
## Run before anything is written, so newly added items get the same scrutiny
## as the first 24. Every check returns the offending rows; all should be empty.
tokens <- function(s) str_split(s, " ")[[1]]
bare <- function(t) str_to_lower(str_remove_all(t, "[^[:alnum:]'-]"))
punct <- function(t) str_remove_all(t, "[[:alnum:]'-]")

paired <- sentences %>%
  filter(version != "original") %>%
  left_join(sentences %>% filter(version == "original") %>%
              select(item_id, length, original = sentence),
            by = c("item_id", "length")) %>%
  rowwise() %>%
  mutate(
    o = list(tokens(original)), x = list(tokens(sentence)),
    permutation = identical(sort(bare(o)), sort(bare(x))),
    punctuation_in_place = identical(punct(o), punct(x)),
    replay = list(reduce(as.integer(str_split(swap_slots, " ")[[1]]),
                         function(w, i) { w[c(i + 1, i + 2)] <- w[c(i + 2, i + 1)]; w },
                         .init = bare(o))),
    slots_replay = identical(replay, bare(x)),
    displaced = list(which(bare(o) != bare(x))),
    one_adjacent_swap = version != "swap-1" ||
      (length(displaced) == 2 && diff(displaced) == 1),
    ## a capital anywhere but the first word must be one the original also has
    ## there, i.e. a proper noun -- catches a moved sentence-initial "The"
    capitals_ok = all(x[-1][str_detect(x[-1], "^[[:upper:]]")] %in%
                        o[-1][str_detect(o[-1], "^[[:upper:]]")])) %>%
  ungroup()

problems <- list(
  "cell missing or duplicated" =
    sentences %>% count(item_id, length, version) %>%
    complete(item_id, length, version, fill = list(n = 0)) %>% filter(n != 1),
  "n_words disagrees with the sentence" =
    sentences %>% filter(n_words != map_int(sentence, ~ length(tokens(.x)))),
  "stray whitespace" =
    sentences %>% filter(str_detect(sentence, "  ") | sentence != str_trim(sentence)),
  "% or : in a sentence (Surveyor escapes these)" =
    sentences %>% filter(str_detect(sentence, "[%:]")),
  "first word not capitalised" =
    sentences %>% filter(!str_detect(sentence, "^[[:upper:]]")),
  "no sentence-final punctuation" =
    sentences %>% filter(!str_detect(sentence, "[.!?]$")),
  "not short < medium < long" =
    sentences %>% filter(version == "original") %>%
    pivot_wider(id_cols = item_id, names_from = length, values_from = n_words) %>%
    filter(!(short < medium & medium < long)),
  "exchange is not a permutation of the original" = filter(paired, !permutation),
  "punctuation moved with a word" = filter(paired, !punctuation_in_place),
  "swap_slots does not reproduce the sentence" = filter(paired, !slots_replay),
  "1-exchange is not a single adjacent swap" = filter(paired, !one_adjacent_swap),
  "capital letter moved mid-sentence" = filter(paired, !capitals_ok),
  "duplicate sentence text" =
    sentences %>% add_count(sentence) %>% filter(n > 1)
)
cat("---- source checks ----\n")
iwalk(problems, ~ cat(sprintf("  %-50s %d\n", .y, nrow(.x))))
bad <- keep(problems, ~ nrow(.x) > 0)
if (length(bad)) {
  iwalk(bad, ~ { cat("\n##", .y, "\n"); print(.x %>% select(any_of(
    c("item_id", "length", "version", "n_words", "sentence")))) })
  stop("fix the source materials before building the survey")
}

## ---- Stimuli ----
materials <- sentences %>%
  transmute(
    experiment = "annot",
    item = item_id, condition,
    sentence,
    prompt = "Click on all the words that are misplaced",
    question_type = "per_word_multiple_choice",
    options = "NA",
    literal_response = "NA",
    config = "NA"
  ) %>%
  arrange(item,
          factor(str_extract(condition, "[^_]+$"), c("short", "medium", "long")),
          factor(str_extract(condition, "^[^_]+"), c("original", "swap-1", "swap-3")))

## Surveyor rejects a stimuli file unless, within each experiment, the number of
## items is a multiple of the number of conditions and every (item, condition)
## pair appears exactly once.
latin <- materials %>%
  summarize(items = n_distinct(item), conditions = n_distinct(condition),
            latin_square_ok = n_distinct(item) %% n_distinct(condition) == 0,
            pairs_unique = n() == n_distinct(paste(item, condition)))
cat("\n---- Latin square ----\n")
print(latin)
if (!latin$latin_square_ok || !latin$pairs_unique)
  stop(sprintf("%d items x %d conditions: Surveyor needs the item count to be a multiple of %d",
               latin$items, latin$conditions, latin$conditions))

write_csv(materials, "materials_en_artificial_exchange_identification.csv", na = "")
