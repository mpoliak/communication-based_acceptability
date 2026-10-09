## =============================================================================
## Run every Supporting Information analysis, in document order: the ten
## scripts below.
##
## Each section is also runnable on its own -- that is the point of the split --
## so this is only for reproducing the whole SI in one go. Each script sources
## setup.R itself, so the data are read once per section rather than shared --
## except sections 4.1-4.5, which need none of setup.R. They read the
## ../data/en_artificial_* folders (4.2 and 4.5 each also read the other's
## responses, to exclude the people who took both) and, for 4.2-4.5, the
## sentence surprisals in ../data/slor_morcela_scores.csv. 4.3 sources 4.1
## (silently) and 4.2 itself, so 4.1 runs twice here; it is deterministic, so
## the second run rewrites identical figures.
##
## The SI curve figure over surprisal (figures/curves_surprisal.png/.pdf) is not
## drawn here: ../analysis/main_analysis.R draws it from the main models.
## =============================================================================

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

scripts <- c(
  "section1_number_of_variables.R",
  "section2_shuffled_condition.R",
  "section3_length.R",
  "section4_binary_error_detection.R",
  "section4.1_exchange_identification_artificial.R",
  "section4.2_binary_error_detection_artificial.R",
  "section4.3_binary_error_detection_artificial_validated.R",
  "section4.4_acceptability_ratings_artificial.R",
  "section4.5_exchange_identification_acceptability_artificial.R",
  "diagnostics_and_extras.R"
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
