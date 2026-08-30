## =============================================================================
## Run every Supporting Information analysis, in document order.
##
## Each section is also runnable on its own -- that is the point of the split --
## so this is only for reproducing the whole SI in one go. Each script sources
## setup.R itself, so the data are read once per section rather than shared.
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
  "section3_surprisal.R",
  "section4_error_awareness.R",
  "diagnostics_and_extras.R"
)

for (script in scripts) {
  message("\n########## ", script, " ##########")
  source(script, echo = FALSE)
}
