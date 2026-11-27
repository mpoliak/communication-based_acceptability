## =============================================================================
## Build the Surveyor questionnaire for the highlighted errors version of the
## English controlled experiment.
##
## Same sentences as data/en_experimental, restricted to original/swap-1/swap-3
## (no shuffled condition). There is no comprehension or binary error detection
## question at all -- participants only rate acceptability. Wherever a swap
## condition displaces a word, that word is wrapped in bold red HTML in the
## displayed sentence, so the error (when present) is marked for the reader
## rather than something they have to find themselves. Writes
## materials_en_highlighted_errors.csv (stimuli) and
## demographics_en_highlighted_errors.csv (survey properties).
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

## ---- Highlighting ----
## Number the tokens of the original sentence 1..n, and do the same for the
## perturbed sentence (they are permutations of one another, so both have n
## tokens). Wherever position i holds a different token than the original's
## position i, that token is a displaced word -- wrap it in the highlight span.
## This is a positional diff, not a per-word dictionary lookup, so it is
## exactly as good as the person who built the swap conditions: it flags every
## position the permutation touched, and matches the "3 swaps" design so it
## flags 6 positions (3 pairwise exchanges) when both members of a pair land
## somewhere new. The comparison is case-sensitive, so an unrelated
## pre-existing quirk in these materials, where a moved word loses its
## capitalization, also flags such a word where it stands at its original
## position (4 words, in the short 3-swap sentences of items 26, 27, 39 and
## 54).

## Inline style, not a CSS class: Surveyor renders `sentence` as raw HTML with
## no guarantee a stylesheet is attached, so the styling has to travel with
## the tag itself.
highlight_word <- function(word) {
  paste0('<b style="color:red;">', word, '</b>')
}

## Wrap every token whose position differs from the original in red bold.
## `original` and `perturbed` must already be equal-length token vectors (true
## for every swap-1/swap-3 sentence here -- they are permutations of the
## original, verified below).
tokenize <- function(sentence) str_split(str_squish(sentence), " ")[[1]]

## The source materials end every perturbed sentence in a period, even when the
## original sentence in that (item, length) pair does not have one (a corpus
## sentence ending "?" or simply missing a full stop -- 21 cases, verified by
## direct inspection: in every one, the perturbed sentence's final token is
## EXACTLY the original's final token with one "." appended, and the word
## itself was not touched by the swap). Comparing raw tokens would flag that
## trailing period as a displaced word. This does not touch multi-word swaps
## that legitimately land a different word in the final slot -- 73 such cases
## exist and stay highlighted, since their final tokens differ by more than a
## trailing period.
final_token_is_punctuation_only <- function(original_last, perturbed_last) {
  identical(perturbed_last, paste0(original_last, "."))
}

highlight_diff <- function(original_tokens, perturbed_tokens) {
  changed <- original_tokens != perturbed_tokens
  n <- length(original_tokens)
  if (changed[n] && final_token_is_punctuation_only(original_tokens[n], perturbed_tokens[n])) {
    changed[n] <- FALSE
  }
  out <- perturbed_tokens
  out[changed] <- vapply(perturbed_tokens[changed], highlight_word, character(1))
  paste(out, collapse = " ")
}

## ---- Sentences: original, swap-1, swap-3 only -- no shuffled condition ----
sentences <- read_csv("../en_experimental/materials_en_experimental.csv",
                      show_col_types = FALSE) %>%
  filter(experiment == "exp_b", !str_starts(condition, "random")) %>%
  select(item, condition, sentence) %>%
  mutate(family = str_remove(condition, "_(short|medium|long)$"),
         length = str_extract(condition, "(short|medium|long)$"))

## One original sentence per item x length, to diff every swap condition
## against and to verify token counts line up before touching any text.
originals <- sentences %>%
  filter(family == "original") %>%
  transmute(item, length, original_sentence = sentence)

sentences <- sentences %>% left_join(originals, by = c("item", "length"))

## Item 12's original sentence (medium/long) has a stray double space in the
## source materials -- a pre-existing quirk, not introduced here. `tokenize()`
## collapses repeated whitespace before splitting so the diff still lines up;
## the double space is untouched in the sentence actually shown to participants
## (the `original` rows below use `sentence` verbatim, never
## `original_sentence`).
stopifnot(all(lengths(lapply(sentences$sentence, tokenize)) ==
              lengths(lapply(sentences$original_sentence, tokenize))))

materials <- sentences %>%
  rowwise() %>%
  mutate(sentence = if (family == "original") sentence else
    highlight_diff(tokenize(original_sentence), tokenize(sentence))) %>%
  ungroup() %>%
  transmute(
    experiment = "exp_b",
    item, condition,
    sentence,
    prompt = "Please rate how acceptable the sentence is (note that errors consisting of exchanged words are shown in red)",
    question_type = "multiple_choice",
    options = acceptability_options,
    literal_response = "NA",
    config = ""
  ) %>%
  arrange(as.integer(item), condition)

write_csv(materials, "materials_en_highlighted_errors.csv", na = "")

## ---- Survey properties (demographics) ----
## Mirrors the demographics files of the binary error detection studies. The
## instructions drop every mention of an error question -- there isn't one --
## and explain what the red bold text means instead.
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
         "at a time.<br><br>Some sentences contain an error created by ",
         "exchanging the positions of two or more words. When a sentence does, ",
         "the exchanged words will be shown in ",
         "<b style=\"color:red;\">bold red</b>. Not every sentence contains an ",
         "error.<br><br>For each sentence, rate how acceptable it sounds to you ",
         "on a scale from 1 (definitely not acceptable) to 7 (definitely ",
         "acceptable). Trust your intuition; there are no right or wrong ",
         "answers, and you should not overthink it.<br><br>Please read each ",
         "sentence carefully and answer every question"),
  "", ""
)

write_csv(demographics, "demographics_en_highlighted_errors.csv", na = "")

## ---- Checks ----
cat("rows written:", nrow(materials), "\n")
materials %>%
  summarize(items = n_distinct(item), conditions = n_distinct(condition),
            latin_square_ok = n_distinct(item) %% n_distinct(condition) == 0,
            pairs_unique = n() == n_distinct(paste(item, condition))) %>%
  print()
cat("condition labels:", paste(sort(unique(materials$condition)), collapse = ", "), "\n")
cat("rows with a highlight span (should be 0 for original, >0 for swap-1/swap-3):\n")
materials %>%
  mutate(family = str_remove(condition, "_(short|medium|long)$"),
         highlighted = str_detect(sentence, "color:red")) %>%
  count(family, highlighted) %>%
  print()
