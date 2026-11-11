## =============================================================================
## SI Section 2: Analyses including the shuffled condition
##
## Repeats the main comparison of all five models (Binary, Gradient,
## Communication, SLOR and MORCELA) with the fully shuffled sentences included.
## Information(s) is sentence surprisal, as in the main analysis. MORCELA's
## weights are learned on these data (written to s2_morcela_weights.csv).
## Produces with_shuffle_figure.png and likelihood_differences_with_shuffled.png
## (and .pdf).
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
## in the main text, they are learned once, on all nine experiments, for the
## pooled fit (here), and separately for each experiment for the by-experiment
## fits (below).
crit_by_items_full_with_morcela <- crit_by_items_full

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

morcela_weights_full_table

write_csv(morcela_weights_full_table, "outputs/s2_morcela_weights.csv")

## ---- Model Comparison with the Shuffled Condition ----
## Repeats the main model comparison on the full data, including the shuffled
## (DL distance >= 4) condition that the main text excludes.
## Produces with_shuffle_figure.png and likelihood_differences_with_shuffled.png.
##
## As in the main analysis, the pooled fit pools all nine experiments (the eight
## languages and controlled English). The controlled English item ids have lost
## their s/m/l length suffix in setup.R, so the three length versions of a base
## item share one item intercept.
if (!file.exists("models/s2_shuffled_surprisal_pooled.RData")) {
  model_grammatical_incl_shuffled <- scam(
    response_num ~ grammatical_numeric + # + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, grammatical_numeric, bs="re"), # RE
    data=crit_by_items_full
  )
  
  model_distance_incl_shuffled <- scam(
    response_num ~ s(dl_distance, bs="mpd") + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, dl_distance, bs="re"), # RE
    data=crit_by_items_full
  )
  
  model_weighted_distance_incl_shuffled <- scam(
    response_num ~ s(weighted_dl_distance_surprisal, bs="mpd") + s(sum_surprisal, bs="mpd") + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, weighted_dl_distance_surprisal, bs="re") + 
      s(language, sum_surprisal, bs="re"), # RE
    data=crit_by_items_full
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

## SLOR and MORCELA, specified as in the main text. They were added after the
## three models above, so they are cached in a file of their own: the three
## models' cache, and every number reported from it, stay as they were.
if (!file.exists("models/s2_shuffled_slor_morcela_pooled.RData")) {
  print(paste(Sys.time(), "Fitting the pooled SLOR model"))
  model_slor_incl_shuffled <- scam(
    response_num ~ s(slor, bs="mpi") + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, slor, bs="re"), # RE
    data=crit_by_items_full_with_morcela,
    control = scam.control(trace = TRUE)
  )

  print(paste(Sys.time(), "Fitting the pooled MORCELA model"))
  model_morcela_incl_shuffled <- scam(
    response_num ~ s(morcela, bs="mpi") + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, morcela, bs="re"), # RE
    data=crit_by_items_full_with_morcela,
    control = scam.control(trace = TRUE)
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
## monotone increasing, as in the main text
stopifnot(
  model_slor_incl_shuffled$smooth[[1]]$label == "s(slor)",
  inherits(model_slor_incl_shuffled$smooth[[1]], "mpi.smooth"),
  model_morcela_incl_shuffled$smooth[[1]]$label == "s(morcela)",
  inherits(model_morcela_incl_shuffled$smooth[[1]], "mpi.smooth")
)


## Inspect models
model_grammatical_incl_shuffled$conv
model_distance_incl_shuffled$conv
model_weighted_distance_incl_shuffled$conv
model_slor_incl_shuffled$conv
model_morcela_incl_shuffled$conv

## Extract the marginal log-likelihood of each model (see the note above: scam's
## own logLik()/AIC()/BIC() are conditional on the fitted random effects).
marginal_logLik(model_grammatical_incl_shuffled)
marginal_logLik(model_distance_incl_shuffled)
marginal_logLik(model_weighted_distance_incl_shuffled)
marginal_logLik(model_slor_incl_shuffled)
marginal_logLik(model_morcela_incl_shuffled)

## Each model's convergence is recorded next to its log-likelihood: `converged`
## is scam's own flag ($conv, whether its penalized IRLS converged), and
## `bfgs_convergence` how its search for the smoothing parameters ended ("Full
## convergence", or e.g. "Iteration limit reached")
incl_shuffled_main_likelihoods <- tibble(
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

write_csv(incl_shuffled_main_likelihoods, "outputs/s2_shuffled_surprisal_likelihoods_pooled.csv")

## Drawn like the main text's five-model figure (../analysis/main_analysis.R),
## with MODEL_COLORS, MODEL_SHAPES and MODEL_LABELS from setup.R
p_incl_shuffled_main_likelihoods <- ggplot(incl_shuffled_main_likelihoods) +
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
  scale_fill_manual(values=MODEL_COLORS) +
  scale_shape_manual(values=MODEL_SHAPES) +
  scale_x_continuous(
    name = ""
  ) +
  ## Five entries do not fit on one row at this width, so the paper's three
  ## models go on the first row and SLOR and MORCELA on the second.
  guides(color = guide_legend(nrow = 2, byrow = TRUE), shape = "none", fill = "none") +
  theme(legend.position="top") +
  scale_y_discrete(name="")

### ---- Analysis By Language ----
if (file.exists("models/s2_shuffled_surprisal_by_language.RData")) {
  load("models/s2_shuffled_surprisal_by_language.RData")
} else {
  
  ## Fit a model for each language
  ## Conditional log-likelihoods, diagnostic only -- see the note in Section 1.
  ## The figures use the marginal values rebuilt after this cache block.
  fits_byLang_full_with_shuffle_stats <- tibble(
    language = character(),
    model = character(),
    log_likelihood_conditional = double()
  )
  
  fits_byLang_full_with_shuffle <- list()
  for (language_iter in levels(crit_by_items$language)) {
    print(
      paste0("Fitting models for ", language_iter)
    )
    ## Get specific data
    crit_by_items_language_specific <- crit_by_items_full %>%
      filter(language == language_iter)
    
    if (language_iter == "en_exp") {
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
    
    fits_byLang_full_with_shuffle_stats <- fits_byLang_full_with_shuffle_stats %>%
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
  
  fits_byLang_full_with_shuffle_stats <- fits_byLang_full_with_shuffle_stats %>%
    mutate(
      language = label_language(language),
      model = case_when(
        model == "binary" ~ "Binary\nGrammaticality",
        model == "gradient" ~ "Gradient\nGrammaticality",
        model == "comm" ~ "Communication-based"
      ) %>%
        factor(levels=c("Binary\nGrammaticality", "Gradient\nGrammaticality", "Communication-based"))
    )
  
  ## Saved *_stats holds CONDITIONAL log-likelihoods; plotted values are marginal.
  save(
    fits_byLang_full_with_shuffle,
    fits_byLang_full_with_shuffle_stats,
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
  for (language_iter in levels(crit_by_items$language)) {
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
## All nine experiments are there, with an increasing SLOR or MORCELA smooth
## (the first smooth of each model)
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
sapply(fits_byLang_full_with_shuffle_5models, function(M) sapply(M, function(m) m$conv))

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

## Stack the combined "All Experiments" panel (top) over the by-language panel
## (bottom), matching the main-text likelihood figure.
likelihood_differences_with_shuffled_figure <- p_incl_shuffled_main_likelihoods / p_language_likelihoods +
  plot_layout(heights = c(1, 3))

ggsave("figures/likelihood_differences_with_shuffled.png",
       plot = likelihood_differences_with_shuffled_figure, width = 8, height = 8)

## Save as pdf, with the fonts embedded (cairo_pdf)
ggsave("figures/likelihood_differences_with_shuffled.pdf",
       plot = likelihood_differences_with_shuffled_figure, width = 8, height = 8,
       device = cairo_pdf)


## ---- Model Predictions Figure ----
## Predictions are population-level: random_effect_labels() drops every item
## effect, including the by-item slopes of the controlled English models.

### ---- Baseline scatterplot (item means) ----
DL_COLORS <- c("0"="#EB7AA3", "1-2"="#1E88E5", "3"="#FFC107", "4+"="#004D40")

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
  lapply(names(fits_byLang_full_with_shuffle), function(lang) {
    predict(
      fits_byLang_full_with_shuffle[[lang]]$binary,
      newdata = get_prediction_grid_surprisal(fits_byLang_full_with_shuffle, lang, "binary", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal),
                                              dl_distance=c(0, 1, 3, 10)),
      type = "response",
      se.fit = TRUE,
      exclude = random_effect_labels(fits_byLang_full_with_shuffle[[lang]]$binary)
    )
  }
  ),
  names(fits_byLang_full_with_shuffle)
)

## Combine list to tibble
grammatical_plot_data <- bind_rows(
  lapply(names(grammatical_predictions), function(lang) {
    preds <- grammatical_predictions[[lang]]
    grid  <- get_prediction_grid_surprisal(fits_byLang_full_with_shuffle, lang, "binary", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal),
                                           dl_distance=c(0, 1, 3, 10))
    
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
      "1-2" = "darkgrey",
      "3" = "darkgrey",
      "4+" = "darkgrey"
    )
  ) +
  ggplot2::guides(color = "none") 

if (interactive()) print(p_grammatical)

#### ---- Distance model ----
#### ---- Predictions from distance model ----
distance_predictions <- setNames(
  lapply(names(fits_byLang_full_with_shuffle), function(lang) {
    predict(
      fits_byLang_full_with_shuffle[[lang]]$gradient,
      newdata = get_prediction_grid_surprisal(fits_byLang_full_with_shuffle, lang, "gradient", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal),
                                              dl_distance=c(0, 1, 3, 10)),
      type = "response",
      se.fit = TRUE,
      exclude = random_effect_labels(fits_byLang_full_with_shuffle[[lang]]$gradient)
    )
  }
  ),
  names(fits_byLang_full_with_shuffle)
)

## Combine list to tibble
distance_plot_data <- bind_rows(
  lapply(names(distance_predictions), function(lang) {
    preds <- distance_predictions[[lang]]
    grid  <- get_prediction_grid_surprisal(fits_byLang_full_with_shuffle, lang, "gradient", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal),
                                           dl_distance=c(0, 1, 3, 10))
    
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
  lapply(names(fits_byLang_full_with_shuffle), function(lang) {
    predict(
      fits_byLang_full_with_shuffle[[lang]]$comm,
      newdata = get_prediction_grid_surprisal(fits_byLang_full_with_shuffle, lang, "comm", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal),
                                              dl_distance=c(0, 1, 3, 10)),
      type = "response",
      se.fit = TRUE,
      exclude = random_effect_labels(fits_byLang_full_with_shuffle[[lang]]$comm)
    )
  }
  ),
  names(fits_byLang_full_with_shuffle)
)

## Combine list to tibble
weighted_distance_plot_data <- bind_rows(
  lapply(names(weighted_distance_predictions), function(lang) {
    preds <- weighted_distance_predictions[[lang]]
    grid  <- get_prediction_grid_surprisal(fits_byLang_full_with_shuffle, lang, "comm", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal),
                                           dl_distance=c(0, 1, 3, 10))
    
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
  filename = "figures/with_shuffle_figure.png",
  plot     = main_figure,
  width    = 25,
  height   = 20,
  units    = "in",
  dpi      = 300
)

## Save as pdf
ggsave(
  filename = "figures/with_shuffle_figure.pdf",
  plot     = main_figure,
  width    = 25,
  height   = 20,
  units    = "in"
)
