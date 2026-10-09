## =============================================================================
## One-off, already run: replaced the Prolific IDs in each experiment_data.csv
## with <language>_<n> before release. Do NOT re-run: the raw IDs are not
## distributed, re-running would renumber the released IDs, and it loops over
## every data/ subfolder, so it fails on the folders that have no
## experiment_data.csv (e.g. the SI follow-up experiments).
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
options(ggplot2.discrete.colour= c("#A31F34", "#8A8B8C"))

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

directories <- list.dirs(full.names=F, recursive=F)
directories

for (language in directories) {
  print(paste0("Anonymizing ", language))
  dat <- read_csv(paste(language, "experiment_data.csv", sep="/"))
  dat$participant_id <- as.factor(dat$participant_id)
  levels(dat$participant_id) <- paste(language, 1:length(levels(dat$participant_id)), sep="_")
  write_csv(dat, paste(language, "experiment_data.csv", sep="/"))
}
