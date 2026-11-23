## =============================================================================
## SI Section 2: Analyses including the shuffled condition
##
## Repeats the main comparison of all five models (Binary, Gradient,
## Communication, SLOR and MORCELA) with the fully shuffled sentences included.
## Information(s) is sentence surprisal, as in the main analysis. MORCELA's
## weights are learned on these data (written to outputs/s2_morcela_weights.csv).
## Produces Figure S3 (figures/with_shuffle_figure.png) and Figure S4
## (figures/likelihood_differences_with_shuffled.png), and writes
## outputs/s2_shuffled_surprisal_likelihoods_pooled.csv and _by_language.csv.
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
## already (setup.R). MORCELA's two weights are learned by the main text's
## recipe (../analysis/main_analysis.R), which fits them to the very ratings the
## score is then used to predict -- here all of them, shuffled condition
## included. So the weights are learned anew on these data rather than taken
## from the main text, just as the other models' smooths are refit on them. As
## in the main text, they are learned once, on the eight language experiments
## (controlled English left out of the pool), for the pooled fit (here), and
## separately for each experiment for the by-experiment fits (below).
## crit_by_items_full_with_morcela is the data of all five pooled fits;
## crit_by_items_full keeps all nine experiments for the by-experiment fits.
crit_by_items_full_with_morcela <- crit_by_items_full %>%
  filter(language != "en_exp") %>%
  droplevels()

morcela_weights_full <- morcela_weights(crit_by_items_full_with_morcela)

crit_by_items_full_with_morcela$morcela <- morcela_score(
  crit_by_items_full_with_morcela, morcela_weights_full
)

## MORCELA's weights, pooled and for each experiment
morcela_weights_full_table <- bind_rows(
  tibble(
    experiment = "all_experiments",
    beta = morcela_weights_full$beta,
    gamma = morcela_weights_full$gamma
  ),
  bind_rows(lapply(levels(crit_by_items_full$language), function(language_iter) {
    weights <- morcela_weights(crit_by_items_full %>% filter(language == language_iter))
    tibble(
      experiment = language_iter,
      beta = weights$beta,
      gamma = weights$gamma
    )
  }))
)

print(morcela_weights_full_table)

write_csv(morcela_weights_full_table, "outputs/s2_morcela_weights.csv")

## ---- Model Comparison with the Shuffled Condition ----
## Repeats the main model comparison on the full data, including the sentences
## with DL distance >= 4 that the main text excludes (nearly all of the shuffled
## condition).
## Produces with_shuffle_figure.png and likelihood_differences_with_shuffled.png.
##
## As in the main analysis, the pooled fit pools the eight language experiments
## (crit_by_items_full_with_morcela, above); controlled English is left out of
## the pool and fit only on its own (below).
if (!file.exists("models/s2_shuffled_surprisal_pooled.RData")) {
  model_grammatical_incl_shuffled <- scam(
    response_num ~ grammatical_numeric + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, grammatical_numeric, bs="re"), # RE
    data=crit_by_items_full_with_morcela
  )
  
  model_distance_incl_shuffled <- scam(
    response_num ~ s(dl_distance, bs="mpd") + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, dl_distance, bs="re"), # RE
    data=crit_by_items_full_with_morcela
  )
  
  model_weighted_distance_incl_shuffled <- scam(
    response_num ~ s(weighted_dl_distance_surprisal, bs="mpd") + s(sum_surprisal, bs="mpd") + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, weighted_dl_distance_surprisal, bs="re") + 
      s(language, sum_surprisal, bs="re"), # RE
    data=crit_by_items_full_with_morcela
  )
  
  
  save(
    model_grammatical_incl_shuffled, 
    model_distance_incl_shuffled, 
    model_weighted_distance_incl_shuffled,
    file="models/s2_shuffled_surprisal_pooled.RData"
  )
  
} else {
  load("models/s2_shuffled_surprisal_pooled.RData")
}
## The three models (cached together) pool the eight language experiments, so an
## old cache of all nine experiments stops the script here
stopifnot(
  nlevels(model_grammatical_incl_shuffled$model$language) == 8,
  !"en_exp" %in% levels(model_grammatical_incl_shuffled$model$language)
)

## SLOR and MORCELA, specified as in the main text. They were added after the
## three models above, so they are cached in a file of their own.
if (!file.exists("models/s2_shuffled_slor_morcela_pooled.RData")) {
  print(paste(Sys.time(), "Fitting the pooled SLOR model"))
  model_slor_incl_shuffled <- scam(
    response_num ~ s(slor, bs="mpi") + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, slor, bs="re"), # RE
    data=crit_by_items_full_with_morcela
  )

  print(paste(Sys.time(), "Fitting the pooled MORCELA model"))
  model_morcela_incl_shuffled <- scam(
    response_num ~ s(morcela, bs="mpi") + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, morcela, bs="re"), # RE
    data=crit_by_items_full_with_morcela
  )

  save(
    model_slor_incl_shuffled,
    model_morcela_incl_shuffled,
    file="models/s2_shuffled_slor_morcela_pooled.RData"
  )

} else {
  load("models/s2_shuffled_slor_morcela_pooled.RData")
}
## The SLOR and MORCELA smooths (the first smooth of each model) must be
## monotone increasing, as in the main text, and the two fits pool the eight
## language experiments (no old cache of all nine)
stopifnot(
  model_slor_incl_shuffled$smooth[[1]]$label == "s(slor)",
  inherits(model_slor_incl_shuffled$smooth[[1]], "mpi.smooth"),
  model_morcela_incl_shuffled$smooth[[1]]$label == "s(morcela)",
  inherits(model_morcela_incl_shuffled$smooth[[1]], "mpi.smooth"),
  nlevels(model_slor_incl_shuffled$model$language) == 8,
  !"en_exp" %in% levels(model_slor_incl_shuffled$model$language),
  nlevels(model_morcela_incl_shuffled$model$language) == 8,
  !"en_exp" %in% levels(model_morcela_incl_shuffled$model$language)
)


## Inspect models
print(model_grammatical_incl_shuffled$conv)
print(model_distance_incl_shuffled$conv)
print(model_weighted_distance_incl_shuffled$conv)
print(model_slor_incl_shuffled$conv)
print(model_morcela_incl_shuffled$conv)

## Extract the marginal log-likelihood of each model (see the note on
## marginal_logLik() in setup.R: scam's own logLik()/AIC()/BIC() are conditional
## on the fitted random effects).
print(marginal_logLik(model_grammatical_incl_shuffled))
print(marginal_logLik(model_distance_incl_shuffled))
print(marginal_logLik(model_weighted_distance_incl_shuffled))
print(marginal_logLik(model_slor_incl_shuffled))
print(marginal_logLik(model_morcela_incl_shuffled))

## Each model's convergence is recorded next to its log-likelihood: `converged`
## is scam's own flag ($conv, whether its penalized IRLS converged), and
## `bfgs_convergence` how its search for the smoothing parameters ended ("Full
## convergence", or e.g. "Iteration limit reached")
pooled_likelihoods <- tibble(
  model = c("Binary", "Gradient", "Communication", "SLOR", "MORCELA"),
  log_likelihood = c(
    marginal_logLik(model_grammatical_incl_shuffled),
    marginal_logLik(model_distance_incl_shuffled),
    marginal_logLik(model_weighted_distance_incl_shuffled),
    marginal_logLik(model_slor_incl_shuffled),
    marginal_logLik(model_morcela_incl_shuffled)
  ),
  converged = c(
    isTRUE(model_grammatical_incl_shuffled$conv),
    isTRUE(model_distance_incl_shuffled$conv),
    isTRUE(model_weighted_distance_incl_shuffled$conv),
    isTRUE(model_slor_incl_shuffled$conv),
    isTRUE(model_morcela_incl_shuffled$conv)
  ),
  bfgs_convergence = c(
    model_grammatical_incl_shuffled$bfgs.info$conv,
    model_distance_incl_shuffled$bfgs.info$conv,
    model_weighted_distance_incl_shuffled$bfgs.info$conv,
    model_slor_incl_shuffled$bfgs.info$conv,
    model_morcela_incl_shuffled$bfgs.info$conv
  )
) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood,
    model = factor(model, levels=MODEL_LEVELS)
  )

write_csv(pooled_likelihoods, "outputs/s2_shuffled_surprisal_likelihoods_pooled.csv")

## Drawn like the main text's five-model figure (../analysis/main_analysis.R),
## with MODEL_COLORS, MODEL_SHAPES and MODEL_LABELS from setup.R
p_pooled_likelihoods <- ggplot(pooled_likelihoods) +
  geom_point(
    aes(
      y="All Languages",
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
  ## Five entries do not fit on one row at this width, so Binary, Gradient and
  ## Communication go on the first row and SLOR and MORCELA on the second. The
  ## colour, fill and shape scales share a name, labels and guide, so ggplot
  ## merges them into one legend whose keys show each model's marker.
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    fill = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme(legend.position="top") +
  scale_y_discrete(name="")

### ---- Analysis By Language ----
if (file.exists("models/s2_shuffled_surprisal_by_language.RData")) {
  load("models/s2_shuffled_surprisal_by_language.RData")
} else {
  
  ## Fit a model for each language. Only the fitted models are cached; the
  ## log-likelihoods are computed below with marginal_logLik().
  fits_byLang_full_with_shuffle <- list()
  for (language_iter in levels(crit_by_items_full$language)) {
    print(
      paste0("Fitting models for ", language_iter)
    )
    ## Get specific data
    crit_by_items_language_specific <- crit_by_items_full %>%
      filter(language == language_iter)
    
    if (language_iter == "en_exp") {
      ## The controlled English item ids have lost their s/m/l length suffix in
      ## setup.R, so the three length versions of a base item share one item
      ## intercept.
      ## Fit models
      binary <- scam(
        response_num ~ grammatical_numeric +
          s(item, bs="re") + s(item, grammatical_numeric, bs="re"),
        data=crit_by_items_language_specific
      )
      
      gradient <- scam(
        response_num ~ s(dl_distance, bs="mpd") + 
          s(item, bs="re") + s(item, dl_distance, bs="re"),
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
        response_num ~ grammatical_numeric +
          s(item, bs="re") ,
        data=crit_by_items_language_specific
      )
      
      gradient <- scam(
        response_num ~ s(dl_distance, bs="mpd") +
          s(item, bs="re"),
        data=crit_by_items_language_specific
      )
      
      comm <- scam(
        response_num ~ s(weighted_dl_distance_surprisal, bs="mpd") + s(sum_surprisal, bs="mpd") +
          s(item, bs="re"),
        data=crit_by_items_language_specific
      )
    }
    
    
    
    fits_byLang_full_with_shuffle[[language_iter]] <- list(
      binary = binary,
      gradient = gradient,
      comm = comm
    )
    
  }

  save(
    fits_byLang_full_with_shuffle,
    file="models/s2_shuffled_surprisal_by_language.RData"
  )
}

### ---- SLOR and MORCELA By Language ----
## SLOR and MORCELA in each experiment, specified as in the main text and cached
## in a file of their own (see the pooled fits above). Only the fitted models are
## cached; the log-likelihoods are computed below with marginal_logLik().
if (file.exists("models/s2_shuffled_slor_morcela_by_language.RData")) {
  load("models/s2_shuffled_slor_morcela_by_language.RData")
} else {

  fits_byLang_full_with_shuffle_slor_morcela <- list()
  for (language_iter in levels(crit_by_items_full$language)) {
    print(
      paste0("Fitting SLOR and MORCELA for ", language_iter)
    )
    ## Get specific data
    crit_by_items_language_specific <- crit_by_items_full %>%
      filter(language == language_iter)

    ## Score MORCELA with this experiment's own weights, learned on these data
    crit_by_items_language_specific$morcela <- morcela_score(
      crit_by_items_language_specific, morcela_weights(crit_by_items_language_specific)
    )

    if (language_iter == "en_exp") {
      ## Fit models
      slor <- scam(
        response_num ~ s(slor, bs="mpi") +
          s(item, bs="re") + s(item, slor, bs="re"),
        data=crit_by_items_language_specific
      )

      morcela <- scam(
        response_num ~ s(morcela, bs="mpi") +
          s(item, bs="re") + s(item, morcela, bs="re"),
        data=crit_by_items_language_specific
      )
    } else {
      ## Fit models
      slor <- scam(
        response_num ~ s(slor, bs="mpi") +
          s(item, bs="re"),
        data=crit_by_items_language_specific
      )

      morcela <- scam(
        response_num ~ s(morcela, bs="mpi") +
          s(item, bs="re"),
        data=crit_by_items_language_specific
      )
    }

    fits_byLang_full_with_shuffle_slor_morcela[[language_iter]] <- list(
      slor = slor,
      morcela = morcela
    )
  }

  save(
    fits_byLang_full_with_shuffle_slor_morcela,
    file="models/s2_shuffled_slor_morcela_by_language.RData"
  )
}
## The two caches hold the same experiments, with an increasing SLOR or MORCELA
## smooth (the first smooth of each model)
stopifnot(
  setequal(names(fits_byLang_full_with_shuffle_slor_morcela), names(fits_byLang_full_with_shuffle)),
  all(sapply(fits_byLang_full_with_shuffle_slor_morcela, function(M) {
    M$slor$smooth[[1]]$label == "s(slor)" && inherits(M$slor$smooth[[1]], "mpi.smooth") &&
      M$morcela$smooth[[1]]$label == "s(morcela)" && inherits(M$morcela$smooth[[1]], "mpi.smooth")
  }))
)

## All five models of each experiment, in the order of the figure legend
fits_byLang_full_with_shuffle_5models <- Map(
  c, fits_byLang_full_with_shuffle,
  fits_byLang_full_with_shuffle_slor_morcela[names(fits_byLang_full_with_shuffle)]
)

## Evaluate Convergence
print(sapply(fits_byLang_full_with_shuffle_5models, function(M) sapply(M, function(m) m$conv)))

## Each fit's convergence is recorded next to its log-likelihood, as for the
## pooled fits. marginal_stats_bylang() lists the fits in the order of the model
## list, experiment by experiment, and so do the two columns added here.
fits_byLang_full_with_shuffle_stats <- marginal_stats_bylang(
  fits_byLang_full_with_shuffle_5models,
  c(binary = "Binary", gradient = "Gradient", comm = "Communication",
    slor = "SLOR", morcela = "MORCELA")) %>%
  mutate(
    converged = unlist(lapply(fits_byLang_full_with_shuffle_5models, function(M) sapply(M, function(m) isTRUE(m$conv)))),
    bfgs_convergence = unlist(lapply(fits_byLang_full_with_shuffle_5models, function(M) sapply(M, function(m) m$bfgs.info$conv)))
  )

fits_byLang_full_with_shuffle_stats <- fits_byLang_full_with_shuffle_stats %>%
  group_by(language) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood
  )

write_csv(
  fits_byLang_full_with_shuffle_stats %>%
    ungroup() %>%
    mutate(across(c(language, model), ~ gsub("\n", " ", .x))),
  "outputs/s2_shuffled_surprisal_likelihoods_by_language.csv"
)

p_language_likelihoods <- ggplot(fits_byLang_full_with_shuffle_stats) +
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

## Stack the combined "All Languages" panel (top) over the by-language panel
## (bottom), matching the main-text likelihood figure.
likelihood_differences_with_shuffled_figure <- p_pooled_likelihoods / p_language_likelihoods +
  plot_layout(heights = c(1, 3))

ggsave("figures/likelihood_differences_with_shuffled.png",
       plot = likelihood_differences_with_shuffled_figure, width = 8, height = 8)

## ---- Model Predictions Figure (Figure S3) ----
## Drawn like the main text's Figure 1 (../analysis/main_analysis.R), at its
## print size (theme_print_8pt, setup.R). The columns are the human judgments
## and then the three models in the order of the other figures: Binary, Gradient,
## Communication. Predictions are population-level: random_effect_labels() drops
## every item effect, including the by-item slopes of the controlled English
## models. The models here are fit to all DL distances, so each panel also has
## a curve for the items at DL >= 4 (see "Predictions for the items at DL >= 4").

### ---- Baseline scatterplot (item means) ----
## plot_data_by_item (setup.R) keeps the items at DL >= 4, in the grey "4+"
## bucket
p_baseline <- ggplot(
  plot_data_by_item,
  aes(x = sum_surprisal, y = response, color = dl_distance_bucket)
) +
  geom_point(alpha = 0.3, size = 0.6, stroke = 0) +
  facet_wrap(~language, ncol = 1, strip.position = "left") +
  labs(x = "Sentence surprisal (bits)", y = NULL) +
  theme_print_8pt +
  theme(
    legend.position = "top",
    ## Facet strips outside, with the language names horizontal. Only this
    ## first column shows them; the other columns blank their strips.
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    panel.spacing = unit(4, "pt")
  ) +
  scale_dl_color()

### ---- GAM smooth over item means ----
p_gam <- p_baseline +
  geom_smooth(
    aes(x = sum_surprisal, y = response, color = dl_distance_bucket),
    method  = "gam",
    formula = y ~ s(x, bs = "tp"),
    se = TRUE,
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  theme(
    legend.position = "none"
  ) +
  ggtitle("Human Judgments") +
  scale_y_continuous(
    limits = YLIM,
    breaks = seq(1, 7, 2),
    name = "Acceptability Ratings (1-7 scale)"
  ) +
  ggplot2::guides(color = "none")

### ---- Model predictions ----
## Each model's predictions over the experiment's range of surprisal, at DL
## distances 0, 1 and 3, one experiment at a time
binary_plot_data <- bind_rows(lapply(names(fits_byLang_full_with_shuffle), function(lang) {
  fit <- fits_byLang_full_with_shuffle[[lang]]$binary
  grid <- get_prediction_grid_surprisal(
    fits_byLang_full_with_shuffle, lang, "binary",
    sum_surprisal = range(crit$sum_surprisal[crit$language == lang])
  )
  preds <- predict(fit, newdata = grid, type = "response", se.fit = TRUE,
                   exclude = random_effect_labels(fit))
  tibble(
    language = label_language(lang),
    sum_surprisal = grid$sum_surprisal,
    dl_distance_bucket = bucket_dl(grid$dl_distance),
    response = preds$fit,
    se = preds$se.fit
  )
}))

gradient_plot_data <- bind_rows(lapply(names(fits_byLang_full_with_shuffle), function(lang) {
  fit <- fits_byLang_full_with_shuffle[[lang]]$gradient
  grid <- get_prediction_grid_surprisal(
    fits_byLang_full_with_shuffle, lang, "gradient",
    sum_surprisal = range(crit$sum_surprisal[crit$language == lang])
  )
  preds <- predict(fit, newdata = grid, type = "response", se.fit = TRUE,
                   exclude = random_effect_labels(fit))
  tibble(
    language = label_language(lang),
    sum_surprisal = grid$sum_surprisal,
    dl_distance_bucket = bucket_dl(grid$dl_distance),
    response = preds$fit,
    se = preds$se.fit
  )
}))

communication_plot_data <- bind_rows(lapply(names(fits_byLang_full_with_shuffle), function(lang) {
  fit <- fits_byLang_full_with_shuffle[[lang]]$comm
  grid <- get_prediction_grid_surprisal(
    fits_byLang_full_with_shuffle, lang, "comm",
    sum_surprisal = range(crit$sum_surprisal[crit$language == lang])
  )
  preds <- predict(fit, newdata = grid, type = "response", se.fit = TRUE,
                   exclude = random_effect_labels(fit))
  tibble(
    language = label_language(lang),
    sum_surprisal = grid$sum_surprisal,
    dl_distance_bucket = bucket_dl(grid$dl_distance),
    response = preds$fit,
    se = preds$se.fit
  )
}))

### ---- Predictions for the items at DL >= 4 ----
## The items at DL >= 4 (nearly all of the shuffled condition) have DL
## distances from 4 to 41, so no single DL distance on the grid stands for them.
## Instead, each model predicts every item x condition cell at DL >= 4 at its
## own DL distance and surprisal, and the "4+" curve smooths these predictions
## over surprisal (with the smooth of the human-judgment panel). The binary model
## needs no such curve: it predicts the same rating at every DL distance above
## 0, so its single error line stands for the items at DL >= 4 too.
gradient_dl4_data <- bind_rows(lapply(names(fits_byLang_full_with_shuffle), function(lang) {
  fit <- fits_byLang_full_with_shuffle[[lang]]$gradient
  ## The item level is a placeholder, as in the grids (item effects are excluded)
  dl4_cells <- crit_by_items_full %>%
    filter(language == lang, dl_distance >= 4) %>%
    mutate(item = factor(fit$model$item[1], levels = levels(fit$model$item)))
  tibble(
    language = label_language(lang),
    sum_surprisal = dl4_cells$sum_surprisal,
    dl_distance_bucket = bucket_dl(dl4_cells$dl_distance),
    response = as.numeric(predict(fit, newdata = dl4_cells, type = "response",
                                  exclude = random_effect_labels(fit)))
  )
}))

communication_dl4_data <- bind_rows(lapply(names(fits_byLang_full_with_shuffle), function(lang) {
  fit <- fits_byLang_full_with_shuffle[[lang]]$comm
  dl4_cells <- crit_by_items_full %>%
    filter(language == lang, dl_distance >= 4) %>%
    mutate(item = factor(fit$model$item[1], levels = levels(fit$model$item)))
  tibble(
    language = label_language(lang),
    sum_surprisal = dl4_cells$sum_surprisal,
    dl_distance_bucket = bucket_dl(dl4_cells$dl_distance),
    response = as.numeric(predict(fit, newdata = dl4_cells, type = "response",
                                  exclude = random_effect_labels(fit)))
  )
}))

### ---- Plot binary grammaticality model ----
p_binary <- p_baseline +
  ggnewscale::new_scale_color() +   # <- reset color scale AFTER baseline
  geom_ribbon(
    data = binary_plot_data,
    aes(
      x = sum_surprisal,
      ymin = response - 2 * se,
      ymax = response + 2 * se,
      group = dl_distance_bucket
    ),
    alpha = 0.25
  ) +
  geom_line(
    data = binary_plot_data,
    aes(
      x = sum_surprisal,
      y = response,
      color = dl_distance_bucket
    ),
    linewidth = 0.5
  ) +
  facet_wrap(~language, ncol = 1) +
  theme(
    legend.position = "none",
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  ggtitle(MODEL_LABELS[1]) +
  scale_y_continuous(limits = YLIM, breaks = NULL) +
  scale_color_manual(values = DL_COLORS_BINARY) +
  ggplot2::guides(color = "none")

### ---- Plot gradient grammaticality model ----
p_gradient <- p_baseline +
  geom_ribbon(
    data = gradient_plot_data,
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
    data = gradient_plot_data,
    aes(
      x = sum_surprisal,
      y = response,
      color = dl_distance_bucket
    ),
    linewidth = 0.5
  ) +
  geom_smooth(
    data = gradient_dl4_data,
    aes(x = sum_surprisal, y = response, color = dl_distance_bucket),
    method  = "gam",
    formula = y ~ s(x, bs = "tp"),
    se = FALSE,
    linewidth = 0.5
  ) +
  facet_wrap(~language, ncol = 1) +
  theme(
    legend.position = "none",
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  ggtitle(MODEL_LABELS[2]) +
  scale_y_continuous(limits = YLIM, breaks = NULL) +
  ggplot2::guides(color = "none")

### ---- Plot communication-based model ----
## Predictions above 7.5 are drawn at 7.5, the top of the panel
p_communication <- p_baseline +
  geom_ribbon(
    data = communication_plot_data %>%
      mutate(
        lower = response - 2 * se,
        upper = pmin(response + 2 * se, 7.5)
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
    data = communication_plot_data %>%
      mutate(response = pmin(response, 7.5)),
    aes(
      x = sum_surprisal,
      y = response,
      color = dl_distance_bucket
    ),
    linewidth = 0.5
  ) +
  geom_smooth(
    data = communication_dl4_data,
    aes(x = sum_surprisal, y = response, color = dl_distance_bucket),
    method  = "gam",
    formula = y ~ s(x, bs = "tp"),
    se = FALSE,
    linewidth = 0.5
  ) +
  facet_wrap(~language, ncol = 1) +
  theme(
    legend.position = "none",
    strip.text = element_blank(),
    strip.background = element_blank()
  ) +
  ggtitle(MODEL_LABELS[3]) +
  scale_y_continuous(limits = YLIM, breaks = NULL) +
  ggplot2::guides(color = "none")

### ---- Combine the panels into Figure S3 ----
## As in the main text's Figure 1: one legend and one x-axis title for the four
## columns, and the print theme for the figure as a whole
prediction_figure <- (p_gam | p_binary | p_gradient | p_communication) +
  plot_layout(guides = "collect", axis_titles = "collect") +
  plot_annotation(theme = theme_print_8pt) &
  theme(legend.position = "bottom",
        legend.margin = margin(0, 0, 0, 0),
        legend.box.spacing = unit(2, "pt"),
        plot.title = element_text(hjust = 0.5, size = rel(1))
  )

if (interactive()) print(prediction_figure)

ggsave(
  filename = "figures/with_shuffle_figure.png",
  plot     = prediction_figure,
  width    = 7,
  height   = 7.5,
  units    = "in",
  dpi      = 300
)
