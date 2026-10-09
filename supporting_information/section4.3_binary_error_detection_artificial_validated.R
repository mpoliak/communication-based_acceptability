## =============================================================================
## SI Section 4.3: Section 4.2 on the sentences validated in the exchange identification task
##
## Re-runs every analysis of section 4.2 (error question + acceptability on the
## artificial materials) on only the sentences whose swap people reliably find
## and whose original people reliably leave alone. The validation comes from
## the exchange identification task of section 4.1, run on a different sample of 100 people: a
## sentence -- one item at one length -- is kept when at least 80% of its
## 1-swap trials located the swap AND at least 80% of its original trials
## clicked nothing. All three versions of a kept sentence (original, 1 swap,
## 3 swaps) are kept, so every analysed sentence still has all its conditions.
##
## Figures: figures/artificial_binary_error_detection_validated_*.png
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

## ---- Accuracy of each sentence in the exchange identification task ----
## Section 4.1 builds `trials` for the exchange identification task: one row per participant x
## sentence version, with `flagged` (clicked anything) and `located` (clicked
## the swapped words and nothing away from them). Its own printed output is
## suppressed here; run section 4.1 directly to see it.
invisible(capture.output(source("section4.1_exchange_identification_artificial.R")))

sentence_accuracy <- trials %>%
  filter(condition_label %in% c("Original", "1 swap")) %>%
  mutate(correct = ifelse(condition_label == "1 swap", located, !flagged)) %>%
  group_by(item, length_level, condition_label) %>%
  summarize(trials = n(), accuracy = mean(correct), .groups = "drop")

keep_sentences <- sentence_accuracy %>%
  group_by(item, length_level) %>%
  filter(all(accuracy >= 0.8)) %>%
  ungroup() %>%
  distinct(item, length_level)

cat("---- sentences validated in the exchange identification task ----\n")
cat("kept", nrow(keep_sentences), "of",
    n_distinct(paste(sentence_accuracy$item, sentence_accuracy$length_level)),
    "item x length combinations\n")
keep_sentences %>% count(length_level, name = "kept") %>% print()
cat("cells sitting exactly on 80% (kept, under the at-least-80% rule):",
    sum(abs(sentence_accuracy$accuracy - 0.8) < 1e-9), "\n")
cat("\ndropped (accuracy in the exchange identification task):\n")
sentence_accuracy %>%
  anti_join(keep_sentences, by = c("item", "length_level")) %>%
  mutate(accuracy = round(accuracy, 2)) %>%
  pivot_wider(id_cols = c(item, length_level), names_from = condition_label,
              values_from = accuracy) %>%
  arrange(item, length_level) %>%
  print(n = Inf)

## ---- Section 4.2 on those sentences ----
## Section 4.2 checks for keep_sentences and figure_prefix, and restricts every
## analysis to these sentences when they are set.
figure_prefix <- "figures/artificial_binary_error_detection_validated_"
source("section4.2_binary_error_detection_artificial.R")

## Clear the switches, so running section 4.2 afterwards in the same session
## analyses every sentence again.
rm(keep_sentences, figure_prefix)
