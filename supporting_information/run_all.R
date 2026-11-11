## =============================================================================
## Run every Supporting Information analysis, in document order: the nine
## scripts below.
##
## Each section is also runnable on its own -- that is the point of the split --
## so this is only for reproducing the whole SI in one go. Each script sources
## setup.R itself, so the data are read once per section rather than shared.
##   - Sections 1-3 repeat the main analyses with a surprisal term in every
##     model (1), with the shuffled condition included (2), and with
##     information measured as length in words instead of surprisal (3).
##   - Section 4 has three studies of constructed sentences: acceptability
##     ratings only (4.1), binary error detection + acceptability ratings (4.2)
##     and exchange identification + acceptability ratings (4.3). They source
##     setup.R through section4_setup.R and read the ../data/en_artificial_*
##     folders (4.2 and 4.3 each also read the other's responses, to exclude
##     the people who took both). section4_model_comparison_figure.R draws the
##     three studies' log-likelihoods from the CSVs that 4.1-4.3 write, so it
##     runs after them.
##   - Section 5 refits SLOR and MORCELA for the English experiment under other
##     language models and unigram models.
##   - descriptives.R is not an SI section: it prints the materials and
##     participant descriptives the paper cites and writes
##     outputs/dl_distance_counts.txt.
##
## Runtime: sections 1-3 cache their model fits under models/ (delete a
## cache to refit it). Sections 1 and 2 took 23 and 21 minutes on 2026-11-02,
## refitting SLOR and MORCELA and loading the other three models from their
## caches; from empty caches they take longer (not timed). Run it plugged in,
## from inside supporting_information/:
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
  "section4.1_acceptability_ratings_only.R",
  "section4.2_binary_error_detection_acceptability_ratings.R",
  "section4.3_exchange_identification_acceptability_ratings.R",
  "section4_model_comparison_figure.R",
  "section5_slor_morcela_settings.R",
  "descriptives.R"
)

## A script that fails does not stop the others: its error is reported, the
## remaining scripts run, and the run stops at the end with the list of failures.
failed <- character(0)
for (script in scripts) {
  message("\n########## ", script, " ##########")
  ok <- tryCatch({
    source(script, echo = FALSE)
    TRUE
  }, error = function(e) {
    message("!!!!! ", script, " failed: ", conditionMessage(e))
    FALSE
  })
  if (!ok) failed <- c(failed, script)
  ## A failed script may have stopped in another working directory
  setwd(get_script_dir())
}

if (length(failed)) {
  stop(length(failed), " of ", length(scripts), " scripts failed: ",
       paste(failed, collapse = ", "))
}
message("\nAll ", length(scripts), " scripts ran.")
