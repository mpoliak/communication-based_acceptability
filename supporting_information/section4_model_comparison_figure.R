## =============================================================================
## SI Section 4: the five models' log-likelihoods in the three studies of
## constructed sentences
##
## One row per study, drawn like the by-experiment panel of the main text's
## log-likelihood figure (p_language_likelihoods in analysis/main_analysis.R):
## each model's marginal log-likelihood minus that of the best model in the same
## study, with sentence surprisal as the information measure. The binary error
## detection + acceptability ratings study is shown on the trials judged
## correctly, and the exchange identification + acceptability ratings study on
## the trials identified correctly; the acceptability ratings only study has no
## second task, so all its trials are used.
##
## Reads the log-likelihoods written by section4.1_acceptability_ratings_only.R,
## section4.2_binary_error_detection_acceptability_ratings.R and
## section4.3_exchange_identification_acceptability_ratings.R, so run those
## first.
##
## Figure: figures/s4_likelihood_differences.png/.pdf
## Numbers drawn: outputs/s4_likelihood_differences.csv
## =============================================================================

## Print warnings as they happen, so that they land next to the step that raised
## them in the log
options(warn = 1)

## Set working directory
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

## Libraries, plotting theme, and the models' colours, shapes and legend labels
## (MODEL_LEVELS, MODEL_LABELS, MODEL_COLORS, MODEL_SHAPES, from setup.R)
source("section4_setup.R")

## ---- The rows of the figure ----
## The trials each study contributes, and its label, from top to bottom. The
## labels are broken over lines, as the main figure breaks its long row label
## ("English\n(controlled)"), with the trials used in parentheses on the last
## line.
figure_rows <- tribble(
  ~task,                                           ~subset,                ~dataset,
  "acceptability_ratings_only",                    "All trials",           "Acceptability ratings only\n(all trials)",
  "binary_error_detection_acceptability_ratings",  "Judged correctly",     "Binary error detection +\nacceptability ratings\n(trials judged correctly)",
  "exchange_identification_acceptability_ratings", "Identified correctly", "Exchange identification +\nacceptability ratings\n(trials identified correctly)"
)

likelihood_files <- paste0("outputs/s4_", figure_rows$task, "_likelihoods.csv")
if (!all(file.exists(likelihood_files))) {
  stop("Run the section 4.1-4.3 scripts first; missing: ",
       paste(likelihood_files[!file.exists(likelihood_files)], collapse = ", "))
}

## Surprisal as the information measure; the length check stays in the CSVs
s4_likelihoods <- bind_rows(lapply(likelihood_files, read_csv, show_col_types = FALSE)) %>%
  filter(information == "Surprisal") %>%
  inner_join(figure_rows, by = c("task", "subset")) %>%
  mutate(
    ## The first study on top: a discrete y axis puts the first level at the
    ## bottom
    dataset = factor(dataset, levels = rev(figure_rows$dataset)),
    model = factor(model, levels = MODEL_LEVELS)
  )

stopifnot(
  nrow(s4_likelihoods) == nrow(figure_rows) * length(MODEL_LEVELS),
  all(s4_likelihoods$converged)
)

s4_likelihoods %>%
  select(task, subset, model, log_likelihood, likelihood_diff) %>%
  mutate(log_likelihood = round(log_likelihood, 1),
         likelihood_diff = round(likelihood_diff, 1)) %>%
  print(n = Inf, width = Inf)

write_csv(
  s4_likelihoods %>%
    arrange(desc(dataset), model) %>%
    select(task, subset, model, log_likelihood, max_likelihood, likelihood_diff,
           n_participants, n_trials, n_cells),
  "outputs/s4_likelihood_differences.csv"
)

## ---- All five sets of trials ----
## The figure shows one set of trials per study; the studies with a second task
## were also analysed on all trials. In each of the five sets: the models from
## best to worst (surprisal as the information measure), and the length check:
## Communication's log-likelihood with length minus that with surprisal
s4_all_likelihoods <- bind_rows(lapply(likelihood_files, read_csv, show_col_types = FALSE))

s4_all_likelihoods %>%
  filter(information == "Surprisal") %>%
  arrange(task, subset, desc(log_likelihood)) %>%
  group_by(task, subset) %>%
  summarize(order = paste(model, collapse = " > "), .groups = "drop") %>%
  print(width = Inf)

s4_all_likelihoods %>%
  filter(model == "Communication") %>%
  select(task, subset, information, log_likelihood) %>%
  pivot_wider(names_from = information, values_from = log_likelihood) %>%
  mutate(length_minus_surprisal = Length - Surprisal) %>%
  print(width = Inf)

## ---- Figure ----
## The by-experiment panel of the main figure, with the legend of that figure's
## top panel, since this figure has no pooled panel to carry it. Five entries
## do not fit on one row at this width, so the paper's three models go on the
## first row and SLOR and MORCELA on the second. The figure is as wide as the
## main one (8 in), so that printed one PNAS column (8.7 cm) wide its smallest
## text (tick labels and legend) comes out at about 6 pt, as in the main
## figure. The long row labels push the panel to the right, and a legend
## centred over the panel would run off the right edge, so the legend is
## centred on the whole plot instead (legend.location = "plot").
p_dataset_likelihoods <- ggplot(s4_likelihoods) +
  geom_point(
    aes(
      y=dataset,
      x=likelihood_diff,
      color=model,
      fill=model,
      shape=model
    ),
    size=3
  ) +
  scale_color_manual(
    name="",
    values=MODEL_COLORS,
    labels=MODEL_LABELS
  ) +
  scale_fill_manual(values=MODEL_COLORS) +
  scale_shape_manual(values=MODEL_SHAPES) +
  scale_x_continuous(
    name = "Log-Likelihood Difference"
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE), shape = "none", fill = "none") +
  theme(legend.position="top", legend.location = "plot") +
  scale_y_discrete(name="")

if (interactive()) print(p_dataset_likelihoods)

ggsave("figures/s4_likelihood_differences.png",
       plot = p_dataset_likelihoods, width = 8, height = 5.5)

## Save as pdf, with the fonts embedded (cairo_pdf)
ggsave("figures/s4_likelihood_differences.pdf",
       plot = p_dataset_likelihoods, width = 8, height = 5.5,
       device = cairo_pdf)

print("Done.")
