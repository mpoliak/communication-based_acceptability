## =============================================================================
## SI Section 2: Analyses including the shuffled condition
##
## Repeats the main comparison with the fully shuffled sentences included.
## Information(s) is sentence surprisal, as in the main analysis.
## Produces with_shuffle_figure.png and likelihood_differences_with_shuffled.png.
## =============================================================================

## Print warnings as they happen
options(warn = 1)

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

## Libraries, plotting theme, helper functions, `crit`, the by-item data frames
## and the shared plotting helpers all live in setup.R.
source("setup.R")

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


## Inspect models
model_grammatical_incl_shuffled$conv
model_distance_incl_shuffled$conv
model_weighted_distance_incl_shuffled$conv

## Extract the marginal log-likelihood of each model (see the note above: scam's
## own logLik()/AIC()/BIC() are conditional on the fitted random effects).
marginal_logLik(model_grammatical_incl_shuffled)
marginal_logLik(model_distance_incl_shuffled)
marginal_logLik(model_weighted_distance_incl_shuffled)

incl_shuffled_main_likelihoods <- tibble(
  model = c("grammatical", "distance", "weighted_distance"),
  log_likelihood = c(
    marginal_logLik(model_grammatical_incl_shuffled),
    marginal_logLik(model_distance_incl_shuffled),
    marginal_logLik(model_weighted_distance_incl_shuffled)
  )
) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood,
    model = factor(model, levels=c("grammatical", "distance", "weighted_distance"))
  )

write_csv(incl_shuffled_main_likelihoods, "outputs/s2_shuffled_surprisal_likelihoods_pooled.csv")

p_incl_shuffled_main_likelihoods <- ggplot(incl_shuffled_main_likelihoods) +
  geom_point(
    aes(
      y="All Experiments",
      x=likelihood_diff,
      color=model,
      shape=model
    ),
    size=3
  ) +
  scale_color_manual(
    name="",
    values=c("#EB7AA3", "#1E88E5", "#004D40"),
    labels=c("Binary\nGrammaticality", "Gradient\nGrammaticality", "Communication")
  ) +
  scale_x_continuous(
    name = ""
  ) +
  guides(shape = "none") +
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

## Evaluate Convergence
sapply(fits_byLang_full_with_shuffle, function(M) sapply(M, function(m) m$conv))

fits_byLang_full_with_shuffle_stats <- marginal_stats_bylang(
  fits_byLang_full_with_shuffle,
  c(binary = "Binary\nGrammaticality", gradient = "Gradient\nGrammaticality", comm = "Communication-based"))

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
      shape=model
    ),
    size=3
  ) +
  scale_color_manual(
    name="",
    values=c("#EB7AA3", "#1E88E5", "#004D40")
  ) +
  scale_x_continuous(
    name = "Log-Likelihood Difference"
  ) +
  theme(legend.position="none") +
  scale_y_discrete(name="",
                   limits = rev)

p_language_likelihoods

## Stack the combined "All Experiments" panel (top) over the by-language panel
## (bottom), matching the main-text likelihood figure.
likelihood_differences_with_shuffled_figure <- p_incl_shuffled_main_likelihoods / p_language_likelihoods +
  plot_layout(heights = c(1, 3))

ggsave("figures/likelihood_differences_with_shuffled.png",
       plot = likelihood_differences_with_shuffled_figure, width = 8, height = 8)


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

p_grammatical

ggsave(
  p_grammatical,
  filename = "figures/with_shuffle_figure_grammatical.png",
  width    = 6,
  height   = 16,
  units    = "in",
  dpi      = 300
)

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

p_distance

ggsave(
  p_distance,
  filename = "figures/with_shuffle_figure_distance.png",
  width    = 6,
  height   = 16,
  units    = "in",
  dpi      = 300
)

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

p_weighted_distance

ggsave(
  p_weighted_distance,
  filename = "figures/with_shuffle_figure_weighted_distance.png",
  width    = 6,
  height   = 16,
  units    = "in",
  dpi      = 300
)

### ---- Combine plots into main figure ----
main_figure = (p_gam |  p_weighted_distance | p_distance |p_grammatical) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
  )

main_figure

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
