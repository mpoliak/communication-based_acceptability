## =============================================================================
## SI Section 1: Controlling for Number of Variables
##
## Adds sentence surprisal as a predictor to the binary- and gradient-
## grammaticality models so every model has the same number of predictors (the
## communication model already has it: its Information(s) is sentence surprisal).
## The main text's two language-model baselines, SLOR and MORCELA, get the same
## surprisal term, so all five models of the main text are compared.
## Produces with_surprisal_term_figure.png and
## likelihood_differences_with_surprisal_term.png (and .pdf).
## =============================================================================

## Print warnings as they happen
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

## Libraries, plotting theme, helper functions, `crit`, the by-item data frames
## and the shared plotting helpers all live in setup.R.
source("setup.R")

## ---- SLOR and MORCELA ----
## The main text's two language-model baselines. SLOR is in the by-item data
## already (setup.R). MORCELA's two weights are learned from the ratings exactly
## as in the main text (../analysis/main_analysis.R): once, on all nine
## experiments, for the pooled fit (here), and separately for each experiment
## for the by-experiment fits (below). Surprisal is not added to that
## regression -- the score would then no longer be MORCELA -- so the weights are
## the main text's.
crit_by_items_with_morcela <- crit_by_items

crit_by_items_with_morcela$morcela <- morcela_score(
  crit_by_items_with_morcela, morcela_weights(crit_by_items_with_morcela)
)

## ---- Main Analysis with Surprisal in all models ----
## Adds sentence surprisal as a predictor to the binary- and gradient-
## grammaticality models, matching the communication model's predictor count,
## and likewise to SLOR and MORCELA (MORCELA then also has its two weights).
## As in the communication model, the surprisal effect is monotone decreasing
## ("mpd"): surprisal is strongly correlated with length, and higher surprisal
## means more potential processing difficulty, never less. The SLOR and MORCELA
## smooths themselves are monotone increasing ("mpi"), as in the main text.
## Produces with_surprisal_term_figure.png and
## likelihood_differences_with_surprisal_term.png.
##
## As in the main analysis, the pooled fit pools all nine experiments (the eight
## languages and controlled English). The controlled English item ids have lost
## their s/m/l length suffix in setup.R, so the three length versions of a base
## item share one item intercept.

if (!file.exists("models/s1_surprisal_term_decreasing_pooled.RData")) {
  fit_overall_surprisal_term_binary <- scam(
    response_num ~ grammatical_numeric + s(sum_surprisal, bs="mpd") + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, grammatical_numeric, bs="re") +
      s(language, sum_surprisal, bs="re"), # RE
    data=crit_by_items
  )

  fit_overall_surprisal_term_gradient <- scam(
    response_num ~ s(dl_distance, bs="mpd") + s(sum_surprisal, bs="mpd") + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, dl_distance, bs="re") +
      s(language, sum_surprisal, bs="re"), # RE
    data=crit_by_items
  )

  fit_overall_surprisal_term_comm <- scam(
    response_num ~ s(weighted_dl_distance_surprisal, bs="mpd") + s(sum_surprisal, bs="mpd") + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, weighted_dl_distance_surprisal, bs="re") +
      s(language, sum_surprisal, bs="re"), # RE
    data=crit_by_items
  )

  save(
    fit_overall_surprisal_term_binary,
    fit_overall_surprisal_term_gradient,
    fit_overall_surprisal_term_comm,
    file="models/s1_surprisal_term_decreasing_pooled.RData"
  )

} else {
  load("models/s1_surprisal_term_decreasing_pooled.RData")
}
## The surprisal smooth (the first smooth of the binary model, the second of the
## gradient model) must be monotone decreasing. Checked after a cache is loaded,
## so a model fit under an older specification is refused.
stopifnot(
  fit_overall_surprisal_term_binary$smooth[[1]]$label == "s(sum_surprisal)",
  inherits(fit_overall_surprisal_term_binary$smooth[[1]], "mpd.smooth"),
  fit_overall_surprisal_term_gradient$smooth[[2]]$label == "s(sum_surprisal)",
  inherits(fit_overall_surprisal_term_gradient$smooth[[2]], "mpd.smooth")
)

## SLOR and MORCELA with the same surprisal term. They were added after the three
## models above, so they are cached in a file of their own: the three models'
## cache, and every number reported from it, stay as they were.
if (!file.exists("models/s1_surprisal_term_slor_morcela_pooled.RData")) {
  print(paste(Sys.time(), "Fitting the pooled SLOR model"))
  fit_overall_surprisal_term_slor <- scam(
    response_num ~ s(slor, bs="mpi") + s(sum_surprisal, bs="mpd") + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, slor, bs="re") +
      s(language, sum_surprisal, bs="re"), # RE
    data=crit_by_items_with_morcela,
    control = scam.control(trace = TRUE)
  )

  print(paste(Sys.time(), "Fitting the pooled MORCELA model"))
  fit_overall_surprisal_term_morcela <- scam(
    response_num ~ s(morcela, bs="mpi") + s(sum_surprisal, bs="mpd") + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, morcela, bs="re") +
      s(language, sum_surprisal, bs="re"), # RE
    data=crit_by_items_with_morcela,
    control = scam.control(trace = TRUE)
  )

  save(
    fit_overall_surprisal_term_slor,
    fit_overall_surprisal_term_morcela,
    file="models/s1_surprisal_term_slor_morcela_pooled.RData"
  )

} else {
  load("models/s1_surprisal_term_slor_morcela_pooled.RData")
}
## The surprisal smooth (the second smooth of each model) must be monotone
## decreasing, and the SLOR and MORCELA smooths (the first) monotone increasing,
## as in the main text
stopifnot(
  fit_overall_surprisal_term_slor$smooth[[2]]$label == "s(sum_surprisal)",
  inherits(fit_overall_surprisal_term_slor$smooth[[2]], "mpd.smooth"),
  fit_overall_surprisal_term_morcela$smooth[[2]]$label == "s(sum_surprisal)",
  inherits(fit_overall_surprisal_term_morcela$smooth[[2]], "mpd.smooth"),
  fit_overall_surprisal_term_slor$smooth[[1]]$label == "s(slor)",
  inherits(fit_overall_surprisal_term_slor$smooth[[1]], "mpi.smooth"),
  fit_overall_surprisal_term_morcela$smooth[[1]]$label == "s(morcela)",
  inherits(fit_overall_surprisal_term_morcela$smooth[[1]], "mpi.smooth")
)

## Evaluate Convergence
fit_overall_surprisal_term_binary$conv
fit_overall_surprisal_term_gradient$conv
fit_overall_surprisal_term_comm$conv
fit_overall_surprisal_term_slor$conv
fit_overall_surprisal_term_morcela$conv

## Extract the marginal log-likelihood of each model. scam's own logLik()/AIC()/
## BIC() are conditional on the fitted random effects, so they do not match the
## values reported in the SI; use marginal_logLik() throughout.
marginal_logLik(fit_overall_surprisal_term_binary)
marginal_logLik(fit_overall_surprisal_term_gradient)
marginal_logLik(fit_overall_surprisal_term_comm)
marginal_logLik(fit_overall_surprisal_term_slor)
marginal_logLik(fit_overall_surprisal_term_morcela)

## How much the surprisal term moves each pooled model's predictions: its
## fitted contribution's maximum minus its minimum over the data, in rating
## points (0: the term was estimated flat)
tibble(
  model = c("Binary", "Gradient", "Communication", "SLOR", "MORCELA"),
  surprisal_term_range = c(
    diff(range(predict(fit_overall_surprisal_term_binary, type = "terms")[, "s(sum_surprisal)"])),
    diff(range(predict(fit_overall_surprisal_term_gradient, type = "terms")[, "s(sum_surprisal)"])),
    diff(range(predict(fit_overall_surprisal_term_comm, type = "terms")[, "s(sum_surprisal)"])),
    diff(range(predict(fit_overall_surprisal_term_slor, type = "terms")[, "s(sum_surprisal)"])),
    diff(range(predict(fit_overall_surprisal_term_morcela, type = "terms")[, "s(sum_surprisal)"]))
  )
) %>%
  print()

## Each model's convergence is recorded next to its log-likelihood: `converged`
## is scam's own flag ($conv, whether its penalized IRLS converged), and
## `bfgs_convergence` how its search for the smoothing parameters ended ("Full
## convergence", or e.g. "Iteration limit reached")
main_likelihoods <- tibble(
  model = c("Binary", "Gradient", "Communication", "SLOR", "MORCELA"),
  log_likelihood = c(
    marginal_logLik(fit_overall_surprisal_term_binary),
    marginal_logLik(fit_overall_surprisal_term_gradient),
    marginal_logLik(fit_overall_surprisal_term_comm),
    marginal_logLik(fit_overall_surprisal_term_slor),
    marginal_logLik(fit_overall_surprisal_term_morcela)
  ),
  converged = c(
    isTRUE(fit_overall_surprisal_term_binary$conv),
    isTRUE(fit_overall_surprisal_term_gradient$conv),
    isTRUE(fit_overall_surprisal_term_comm$conv),
    isTRUE(fit_overall_surprisal_term_slor$conv),
    isTRUE(fit_overall_surprisal_term_morcela$conv)
  ),
  bfgs_convergence = c(
    fit_overall_surprisal_term_binary$bfgs.info$conv,
    fit_overall_surprisal_term_gradient$bfgs.info$conv,
    fit_overall_surprisal_term_comm$bfgs.info$conv,
    fit_overall_surprisal_term_slor$bfgs.info$conv,
    fit_overall_surprisal_term_morcela$bfgs.info$conv
  )
) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood,
    model = factor(model, levels=MODEL_LEVELS)
  )

write_csv(main_likelihoods, "outputs/s1_surprisal_term_likelihoods_pooled.csv")

## The change in log-likelihood from adding the surprisal term: these fits minus
## the main text's (../analysis/outputs, written by ../analysis/main_analysis.R).
## Communication has the term already, so it is the main text's model.
main_likelihoods %>%
  mutate(model = as.character(model)) %>%
  left_join(
    read_csv("../analysis/outputs/likelihoods_all_experiments.csv", show_col_types = FALSE) %>%
      select(model, main_text_log_likelihood = log_likelihood),
    by = "model"
  ) %>%
  mutate(change_from_surprisal_term = log_likelihood - main_text_log_likelihood) %>%
  select(model, log_likelihood, main_text_log_likelihood, change_from_surprisal_term) %>%
  print(width = Inf)

## Drawn like the main text's five-model figure (../analysis/main_analysis.R),
## with MODEL_COLORS, MODEL_SHAPES and MODEL_LABELS from setup.R
p_main_likelihoods <- ggplot(main_likelihoods) +
  geom_point(
    aes(
      y="All Experiments",
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
  scale_fill_manual(name="", values=MODEL_COLORS, labels=MODEL_LABELS) +
  scale_shape_manual(name="", values=MODEL_SHAPES, labels=MODEL_LABELS) +
  scale_x_continuous(
    name = ""
  ) +
  ## Five entries do not fit on one row at this width, so the paper's three
  ## models go on the first row and SLOR and MORCELA on the second. The colour,
  ## fill and shape scales share a name, labels and guide, so ggplot merges them
  ## into one legend whose keys show each model's marker.
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    fill = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme(legend.position="top") +
  scale_y_discrete(name="")

if (interactive()) print(p_main_likelihoods)

### ---- Analysis By Language ----
if (file.exists("models/s1_surprisal_term_decreasing_by_language.RData")) {
  load("models/s1_surprisal_term_decreasing_by_language.RData")
} else {

  ## Fit a model for each language
  ## NB: this in-loop tibble records the CONDITIONAL log-likelihood (scam's own
  ## logLik(), which plugs in the fitted random effects). It is kept only as a
  ## diagnostic; the figures below use the MARGINAL values rebuilt with
  ## marginal_stats_bylang() after this cache block.
  fits_byLang_surprisal_term_stats <- tibble(
    language = character(),
    model = character(),
    log_likelihood_conditional = double()
  )

  fits_byLang_surprisal_term <- list()
  for (language_iter in levels(crit_by_items$language)) {
    print(
      paste0("Fitting models for ", language_iter)
    )
    ## Get specific data
    crit_by_items_language_specific <- crit_by_items %>%
      filter(language == language_iter)

    if (language_iter == "en_exp") {
      ## Fit models. Length, and so surprisal, varies within item here, so the
      ## by-item surprisal slopes are identifiable.
      binary <- scam(
        response_num ~ grammatical_numeric + s(sum_surprisal, bs="mpd") +
          s(item, bs="re") + s(item, grammatical_numeric, bs="re") + s(item, sum_surprisal, bs="re"),
        data=crit_by_items_language_specific
      )

      gradient <- scam(
        response_num ~ s(dl_distance, bs="mpd") + s(sum_surprisal, bs="mpd") +
          s(item, bs="re") + s(item, dl_distance, bs="re") + s(item, sum_surprisal, bs="re"),
        data=crit_by_items_language_specific
      )

      comm <- scam(
        response_num ~ s(weighted_dl_distance_surprisal, bs="mpd") + s(sum_surprisal, bs="mpd") +
          s(item, bs="re") + s(item, weighted_dl_distance_surprisal, bs="re") + s(item, sum_surprisal, bs="re"),
        data=crit_by_items_language_specific
      )
    } else {
      ## Fit models
      binary <- scam(
        response_num ~ grammatical_numeric + s(sum_surprisal, bs="mpd") +
          s(item, bs="re") ,
        data=crit_by_items_language_specific
      )

      gradient <- scam(
        response_num ~ s(dl_distance, bs="mpd") + s(sum_surprisal, bs="mpd") +
          s(item, bs="re"),
        data=crit_by_items_language_specific
      )

      comm <- scam(
        response_num ~ s(weighted_dl_distance_surprisal, bs="mpd") + s(sum_surprisal, bs="mpd") +
          s(item, bs="re"),
        data=crit_by_items_language_specific
      )
    }



    fits_byLang_surprisal_term[[language_iter]] <- list(
      binary = binary,
      gradient = gradient,
      comm = comm
    )

    fits_byLang_surprisal_term_stats <- fits_byLang_surprisal_term_stats %>%
      add_row(
        language = language_iter,
        model = "binary",
        log_likelihood_conditional = binary %>% logLik() %>% as.numeric()
      ) %>%
      add_row(
        language = language_iter,
        model = "gradient",
        log_likelihood_conditional = gradient %>% logLik() %>% as.numeric()
      ) %>%
      add_row(
        language = language_iter,
        model = "comm",
        log_likelihood_conditional = comm %>% logLik() %>% as.numeric()
      )

  }

  fits_byLang_surprisal_term_stats <- fits_byLang_surprisal_term_stats %>%
    mutate(
      language = label_language(language),
      model = case_when(
        model == "binary" ~ "Binary\nGrammaticality",
        model == "gradient" ~ "Gradient\nGrammaticality",
        model == "comm" ~ "Communication-based"
      ) %>%
        factor(levels=c("Binary\nGrammaticality", "Gradient\nGrammaticality", "Communication-based"))
    )

  ## The saved *_stats tibble holds CONDITIONAL log-likelihoods
  ## (log_likelihood_conditional); the plotted values are the marginal ones
  ## recomputed immediately after this block.
  save(
    fits_byLang_surprisal_term,
    fits_byLang_surprisal_term_stats,
    file="models/s1_surprisal_term_decreasing_by_language.RData"
  )
}
## The same check of the surprisal smooth in every experiment
stopifnot(all(sapply(fits_byLang_surprisal_term, function(M) {
  M$binary$smooth[[1]]$label == "s(sum_surprisal)" && inherits(M$binary$smooth[[1]], "mpd.smooth") &&
    M$gradient$smooth[[2]]$label == "s(sum_surprisal)" && inherits(M$gradient$smooth[[2]], "mpd.smooth")
})))

### ---- SLOR and MORCELA By Language ----
## The same surprisal term added to SLOR and MORCELA in each experiment, cached
## in a file of their own (see the pooled fits above). Only the fitted models
## are cached; the log-likelihoods are computed below with marginal_logLik().
if (file.exists("models/s1_surprisal_term_slor_morcela_by_language.RData")) {
  load("models/s1_surprisal_term_slor_morcela_by_language.RData")
} else {

  fits_byLang_surprisal_term_slor_morcela <- list()
  for (language_iter in levels(crit_by_items$language)) {
    print(
      paste0("Fitting SLOR and MORCELA for ", language_iter)
    )
    ## Get specific data
    crit_by_items_language_specific <- crit_by_items %>%
      filter(language == language_iter)

    ## Score MORCELA with this experiment's own weights, as in the main text
    crit_by_items_language_specific$morcela <- morcela_score(
      crit_by_items_language_specific, morcela_weights(crit_by_items_language_specific)
    )

    if (language_iter == "en_exp") {
      ## Fit models. Length, and so surprisal, varies within item here, so the
      ## by-item surprisal slopes are identifiable.
      slor <- scam(
        response_num ~ s(slor, bs="mpi") + s(sum_surprisal, bs="mpd") +
          s(item, bs="re") + s(item, slor, bs="re") + s(item, sum_surprisal, bs="re"),
        data=crit_by_items_language_specific
      )

      morcela <- scam(
        response_num ~ s(morcela, bs="mpi") + s(sum_surprisal, bs="mpd") +
          s(item, bs="re") + s(item, morcela, bs="re") + s(item, sum_surprisal, bs="re"),
        data=crit_by_items_language_specific
      )
    } else {
      ## Fit models
      slor <- scam(
        response_num ~ s(slor, bs="mpi") + s(sum_surprisal, bs="mpd") +
          s(item, bs="re"),
        data=crit_by_items_language_specific
      )

      morcela <- scam(
        response_num ~ s(morcela, bs="mpi") + s(sum_surprisal, bs="mpd") +
          s(item, bs="re"),
        data=crit_by_items_language_specific
      )
    }

    fits_byLang_surprisal_term_slor_morcela[[language_iter]] <- list(
      slor = slor,
      morcela = morcela
    )
  }

  save(
    fits_byLang_surprisal_term_slor_morcela,
    file="models/s1_surprisal_term_slor_morcela_by_language.RData"
  )
}
## All nine experiments are there, with the decreasing surprisal smooth (the
## second smooth of each model) and an increasing SLOR or MORCELA smooth (the
## first)
stopifnot(
  setequal(names(fits_byLang_surprisal_term_slor_morcela), names(fits_byLang_surprisal_term)),
  all(sapply(fits_byLang_surprisal_term_slor_morcela, function(M) {
    M$slor$smooth[[2]]$label == "s(sum_surprisal)" && inherits(M$slor$smooth[[2]], "mpd.smooth") &&
      M$morcela$smooth[[2]]$label == "s(sum_surprisal)" && inherits(M$morcela$smooth[[2]], "mpd.smooth") &&
      M$slor$smooth[[1]]$label == "s(slor)" && inherits(M$slor$smooth[[1]], "mpi.smooth") &&
      M$morcela$smooth[[1]]$label == "s(morcela)" && inherits(M$morcela$smooth[[1]], "mpi.smooth")
  }))
)

## All five models of each experiment, in the order of the figure legend
fits_byLang_surprisal_term_5models <- Map(
  c, fits_byLang_surprisal_term,
  fits_byLang_surprisal_term_slor_morcela[names(fits_byLang_surprisal_term)]
)

## Evaluate Convergence
sapply(fits_byLang_surprisal_term_5models, function(M) sapply(M, function(m) m$conv))

### ---- Figure based on individual-language likelihoods
## Each fit's convergence is recorded next to its log-likelihood, as for the
## pooled fits. marginal_stats_bylang() lists the fits in the order of the model
## list, experiment by experiment, and so do the two columns added here.
fits_byLang_surprisal_term_stats <- marginal_stats_bylang(
  fits_byLang_surprisal_term_5models,
  c(binary = "Binary", gradient = "Gradient", comm = "Communication",
    slor = "SLOR", morcela = "MORCELA")) %>%
  mutate(
    converged = unlist(lapply(fits_byLang_surprisal_term_5models, function(M) sapply(M, function(m) isTRUE(m$conv)))),
    bfgs_convergence = unlist(lapply(fits_byLang_surprisal_term_5models, function(M) sapply(M, function(m) m$bfgs.info$conv)))
  )

fits_byLang_surprisal_term_stats <- fits_byLang_surprisal_term_stats %>%
  group_by(language) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood
  )

write_csv(
  fits_byLang_surprisal_term_stats %>%
    ungroup() %>%
    mutate(across(c(language, model), ~ gsub("\n", " ", .x))),
  "outputs/s1_surprisal_term_likelihoods_by_language.csv"
)

## The same change in log-likelihood in each experiment
fits_byLang_surprisal_term_stats %>%
  ungroup() %>%
  mutate(language = as.character(language), model = as.character(model)) %>%
  left_join(
    read_csv("../analysis/outputs/likelihoods_by_language.csv", show_col_types = FALSE) %>%
      mutate(language = as.character(label_language(experiment))) %>%
      select(language, model, main_text_log_likelihood = log_likelihood),
    by = c("language", "model")
  ) %>%
  mutate(change_from_surprisal_term = log_likelihood - main_text_log_likelihood) %>%
  select(language, model, change_from_surprisal_term) %>%
  pivot_wider(names_from = model, values_from = change_from_surprisal_term) %>%
  print(width = Inf)

p_language_likelihoods <- ggplot(fits_byLang_surprisal_term_stats) +
  geom_point(
    aes(
      y=language,
      x=likelihood_diff,
      color=model,
      fill=model,
      shape=model
    ),
    size=3
  ) +
  scale_color_manual(
    name="",
    values=MODEL_COLORS
  ) +
  scale_fill_manual(values=MODEL_COLORS) +
  scale_shape_manual(values=MODEL_SHAPES) +
  scale_x_continuous(
    name = "Log-Likelihood Difference"
  ) +
  theme(legend.position="none") +
  scale_y_discrete(name="",
                   limits = rev)

if (interactive()) print(p_language_likelihoods)

## Stack the combined "All Experiments" panel (top) over the by-language panel
## (bottom), matching the main-text likelihood figure.
likelihood_differences_with_surprisal_term_figure <- p_main_likelihoods / p_language_likelihoods +
  plot_layout(heights = c(1, 3))

ggsave("figures/likelihood_differences_with_surprisal_term.png",
       plot = likelihood_differences_with_surprisal_term_figure, width = 8, height = 8)

## Save as pdf, with the fonts embedded (cairo_pdf)
ggsave("figures/likelihood_differences_with_surprisal_term.pdf",
       plot = likelihood_differences_with_surprisal_term_figure, width = 8, height = 8,
       device = cairo_pdf)

## ---- Model Predictions Figure ----
## Predictions are population-level: random_effect_labels() drops every item
## effect, including the by-item slopes of the controlled English models.

### ---- Baseline scatterplot (item means) ----
p_baseline <- ggplot(
  plot_data_by_item,
  aes(x = sum_surprisal, y = response, color = dl_distance_bucket)
) +
  geom_point(alpha = 0.3) +
  facet_wrap(~language, ncol = 1, strip.position = "left") +
  labs(x="Sentence Surprisal (bits)", y="Rating") +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),

    # put facet strips on the outside + rotate text
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90)
  ) +
  scale_y_continuous(
    limits = YLIM,
    breaks = NULL,
    labels = NULL
  ) +
  scale_dl_color()

### ---- GAM smooth over item means ----
p_gam <- p_baseline +
  geom_smooth(
    aes(x = sum_surprisal, y = response, color = dl_distance_bucket),
    method  = "gam",
    formula = y ~ s(x, bs = "tp"),
    se = TRUE,
    linewidth = 1.5,
    show.legend = FALSE
  ) +
  theme(
    legend.position = "none"
  ) +
  ggtitle("Human Judgments") +
  scale_y_continuous(
    limits = YLIM,
    breaks = seq(1, 7, 2),
    labels = seq(1, 7, 2),
    name = "Acceptability Ratings (1-7 scale)"
  ) +
  ggplot2::guides(color = "none")

### ---- Model-based plots ----
#### ---- Grammatical model ----
#### ---- Predictions from grammatical model ----
grammatical_predictions <- setNames(
  lapply(names(fits_byLang_surprisal_term), function(lang) {
    predict(
      fits_byLang_surprisal_term[[lang]]$binary,
      newdata = get_prediction_grid_surprisal(fits_byLang_surprisal_term, lang, "binary", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal)),
      type = "response",
      se.fit = TRUE,
      exclude = random_effect_labels(fits_byLang_surprisal_term[[lang]]$binary)
    )
  }
  ),
  names(fits_byLang_surprisal_term)
)

## Combine list to tibble
grammatical_plot_data <- bind_rows(
  lapply(names(grammatical_predictions), function(lang) {
    preds <- grammatical_predictions[[lang]]
    grid  <- get_prediction_grid_surprisal(fits_byLang_surprisal_term, lang, "binary", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal))

    tibble(
      language = label_language(lang),
      sum_surprisal = grid$sum_surprisal,
      dl_distance_bucket = bucket_dl(grid$dl_distance),
      response = preds$fit,
      se = preds$se.fit
    )
  })
)

##### ---- Plot grammatical model ----
p_grammatical <- p_baseline +
  ggnewscale::new_scale_color() +   # <- reset color scale AFTER baseline
  geom_ribbon(
    data=grammatical_plot_data,
    aes(
      x = sum_surprisal,
      ymin = response - 2 * se,
      ymax = response + 2 * se,
      group = dl_distance_bucket
    ),
    alpha = 0.25
  ) +
  geom_line(
    data=grammatical_plot_data,
    aes(
      x = sum_surprisal,
      y = response,
      color = dl_distance_bucket
    ),
    linewidth = 1.5
  ) +
  facet_wrap(~language, ncol = 1) +
  labs(
    x = "Sentence Surprisal (bits)",
    y = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  ggtitle("Binary Grammaticality Model") +
  scale_y_continuous(
    limits = YLIM,
    breaks = NULL,
    labels = NULL
  ) +
  scale_color_manual(
    values = c(
      "0" = "#EB7AA3",
      "1-2" = "#06402B",
      "3" = "#06402B",
      "4+" = "#06402B"
    )
  ) +
  ggplot2::guides(color = "none")

if (interactive()) print(p_grammatical)

#### ---- Distance model ----
#### ---- Predictions from distance model ----
distance_predictions <- setNames(
  lapply(names(fits_byLang_surprisal_term), function(lang) {
    predict(
      fits_byLang_surprisal_term[[lang]]$gradient,
      newdata = get_prediction_grid_surprisal(fits_byLang_surprisal_term, lang, "gradient", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal)),
      type = "response",
      se.fit = TRUE,
      exclude = random_effect_labels(fits_byLang_surprisal_term[[lang]]$gradient)
    )
  }
  ),
  names(fits_byLang_surprisal_term)
)

## Combine list to tibble
distance_plot_data <- bind_rows(
  lapply(names(distance_predictions), function(lang) {
    preds <- distance_predictions[[lang]]
    grid  <- get_prediction_grid_surprisal(fits_byLang_surprisal_term, lang, "gradient", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal))

    tibble(
      language = label_language(lang),
      sum_surprisal = grid$sum_surprisal,
      dl_distance_bucket = bucket_dl(grid$dl_distance),
      response = preds$fit,
      se = preds$se.fit
    )
  })
)

##### ---- Plot distance model ----
p_distance <- p_baseline +
  geom_ribbon(
    data=distance_plot_data,
    aes(
      x = sum_surprisal,
      ymin = response - 2 * se,
      ymax = response + 2 * se,
      group = dl_distance_bucket
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(
    data=distance_plot_data,
    aes(
      x = sum_surprisal,
      y = response,
      color = dl_distance_bucket
    ),
    linewidth = 1.5
  ) +
  facet_wrap(~language, ncol = 1) +
  labs(
    x = "Sentence Surprisal (bits)",
    y = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  ggtitle("Gradient Grammaticality Model") +
  scale_y_continuous(
    limits = YLIM,
    breaks = NULL,
    labels = NULL
  ) +
  ggplot2::guides(color = "none")

if (interactive()) print(p_distance)

#### ---- Weighted distance model ----
#### ---- Predictions from weighted distance model ----
weighted_distance_predictions <- setNames(
  lapply(names(fits_byLang_surprisal_term), function(lang) {
    predict(
      fits_byLang_surprisal_term[[lang]]$comm,
      newdata = get_prediction_grid_surprisal(fits_byLang_surprisal_term, lang, "comm", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal)),
      type = "response",
      se.fit = TRUE,
      exclude = random_effect_labels(fits_byLang_surprisal_term[[lang]]$comm)
    )
  }
  ),
  names(fits_byLang_surprisal_term)
)

## Combine list to tibble
weighted_distance_plot_data <- bind_rows(
  lapply(names(weighted_distance_predictions), function(lang) {
    preds <- weighted_distance_predictions[[lang]]
    grid  <- get_prediction_grid_surprisal(fits_byLang_surprisal_term, lang, "comm", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal))

    tibble(
      language = label_language(lang),
      sum_surprisal = grid$sum_surprisal,
      dl_distance_bucket = bucket_dl(grid$dl_distance),
      response = preds$fit,
      se = preds$se.fit
    )
  })
)

##### ---- Plot distance model ----
p_weighted_distance <- p_baseline +
  geom_ribbon(
    data=weighted_distance_plot_data %>%
      mutate(
        lower = response - 2*se,
        upper = response + 2*se,
        upper = ifelse(upper > 7.5, 7.5, upper),
      ),
    aes(
      x = sum_surprisal,
      ymin = lower,
      ymax = upper,
      group = dl_distance_bucket
    ),
    alpha = 0.25,
    color = NA
  ) +
  geom_line(
    data=weighted_distance_plot_data %>%
      mutate(
        response = ifelse(response > 7.5, 7.5, response)
      ),
    aes(
      x = sum_surprisal,
      y = response,
      color = dl_distance_bucket
    ),
    linewidth = 1.5
  ) +
  facet_wrap(~language, ncol = 1) +
  labs(
    x = "Sentence Surprisal (bits)",
    y = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  ggtitle("Communication Model") +
  scale_y_continuous(
    limits = YLIM,
    breaks = NULL,
    labels = NULL
  ) +
  ggplot2::guides(color = "none")

if (interactive()) print(p_weighted_distance)

### ---- Combine plots into main figure ----
main_figure = (p_gam |  p_weighted_distance | p_distance |p_grammatical) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
  )

if (interactive()) print(main_figure)

ggsave(
  filename = "figures/with_surprisal_term_figure.png",
  plot     = main_figure,
  width    = 25,
  height   = 20,
  units    = "in",
  dpi      = 300
)

## Save as pdf
ggsave(
  filename = "figures/with_surprisal_term_figure.pdf",
  plot     = main_figure,
  width    = 25,
  height   = 20,
  units    = "in"
)
