## =============================================================================
## Run every Supporting Information analysis, in document order: the six
## scripts below.
##
## Each section is also runnable on its own -- that is the point of the split --
## so this is only for reproducing the whole SI in one go. Each script runs in a
## separate R process (Rscript), exactly as when it is run alone, and writes its
## output to its own log, ../logs/<script name>.log. Each script sources setup.R
## itself, so the data are read once per section rather than shared.
##   - Sections 1-3 repeat the main analyses with a surprisal term in every
##     model (1), with the shuffled condition included (2), and with
##     information measured as length in tokens instead of surprisal (3).
##   - Section 4 refits SLOR, MORCELA and the communication model for the
##     English experiment under other language models and unigram models.
##   - Section 5 (one script) has the six follow-up experiments: binary error
##     detection + acceptability ratings on the original materials of the
##     English experiment (5.1, ../data/en_binary_error_detection); on the
##     length-controlled materials of the English (controlled) experiment,
##     binary error detection + acceptability ratings (5.2,
##     ../data/en_experimental_binary_error_detection) and acceptability
##     ratings with the errors highlighted (5.3, ../data/en_highlighted_errors);
##     and the three studies of constructed sentences, acceptability ratings
##     only (5.4), binary error detection + acceptability ratings (5.5) and
##     exchange identification + acceptability ratings (5.6), read from the
##     ../data/en_artificial_* folders. It ends with the figure of the five
##     models' log-likelihoods in the six studies, with the main text's
##     English and English (controlled) experiments.
##   - descriptives.R is not an SI section: it prints the sentence-length range
##     and median that the Methods report, plus unreported checks, and writes
##     outputs/dl_distance_counts.txt.
##
## Sections 1-3 cache their model fits under models/ (delete a cache to refit
## it). A fresh clone has no caches, so its first run fits every model, which
## takes long: run it plugged in, from inside supporting_information/
## (caffeinate keeps a Mac awake; leave it out elsewhere):
##   mkdir -p ../logs
##   LC_ALL=en_US.UTF-8 caffeinate -i Rscript run_all.R > ../logs/run_all.log 2>&1
##
## The figure of the five models' log-likelihoods
## (../analysis/figures/likelihood_differences.png), now in the main text, is
## not drawn here: ../analysis/main_analysis.R draws it from the main models.
## =============================================================================

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  ## Only while the path still resolves: once a script has moved into its
  ## folder, a relative --file path no longer does, and getwd() is the folder
  if (length(file_arg) && file.exists(sub("^--file=", "", file_arg))) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    return(dirname(rstudioapi::getActiveDocumentContext()$path))
  }
  getwd()
}
setwd(get_script_dir())

scripts <- c(
  "section1_number_of_variables.R",
  "section2_shuffled_condition.R",
  "section3_length.R",
  "section4_slor_morcela_settings.R",
  "section5_awareness_of_grammatical_errors.R",
  "descriptives.R"
)

## Each script runs in its own R process, so no script can change or remove
## another's objects or this script's. A script that fails does not stop the
## others: its failure is reported (the error is at the end of its log), the
## remaining scripts run, and the run stops at the end with the list of failures.
dir.create("../logs", showWarnings = FALSE)
failed <- character(0)
for (script in scripts) {
  log_file <- file.path("..", "logs", sub("\\.R$", ".log", script))
  message("\n########## ", script, " (log: ", log_file, ") ##########")
  ok <- system2(
    file.path(R.home("bin"), "Rscript"), script,
    stdout = log_file, stderr = log_file
  ) == 0
  if (!ok) {
    message("!!!!! ", script, " failed; see ", log_file)
    failed <- c(failed, script)
  }
}

if (length(failed)) {
  stop(length(failed), " of ", length(scripts), " scripts failed: ",
       paste(failed, collapse = ", "))
}
message("\nAll ", length(scripts), " scripts ran.")
