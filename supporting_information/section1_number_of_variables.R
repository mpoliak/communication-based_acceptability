## =============================================================================
## SI Section 1: Controlling for Number of Variables
##
## Adds utterance length as a predictor to the binary- and gradient-
## grammaticality models so every model has the same number of predictors.
## Produces withLength_figure.png and likelihood_differences_with_length.png.
## =============================================================================

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

## ---- Main Analysis with Length in all models ----
## Adds utterance length (P-spline) as a predictor to the binary- and gradient-
## grammaticality models, matching the communication model's predictor count.
## Produces withLength_figure.png and likelihood_differences_with_length.png.

if (!file.exists("models/fits_overall_full_withLength.RData")) {
  fit_overall_full_withLength_binary <- scam(
    ## length enters as an unconstrained P-spline ("ps"), not a monotone one, so
    ## the grammar-based models are free to fit length in either direction
    response_num ~ grammatical_numeric + s(n_words, bs="ps") + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, grammatical_numeric, bs="re") +
      s(language, n_words, bs="re"), # RE
    data=crit_by_items %>% 
      filter(language != "en_exp")
  )
  
  fit_overall_full_withLength_gradient <- scam(
    ## length enters as an unconstrained P-spline ("ps"), not a monotone one, so
    ## the grammar-based models are free to fit length in either direction
    response_num ~ s(dl_distance, bs="mpd") + s(n_words, bs="ps") + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, dl_distance, bs="re") +
      s(language, n_words, bs="re"), # RE
    data=crit_by_items %>% 
      filter(language != "en_exp")
  )
  
  fit_overall_full_withLength_comm <- scam(
    response_num ~ s(weighted_dl_distance, bs="mpd") + s(n_words, bs="mpd") + # Fixed effects 
      s(item, bs="re") + s(language, bs="re") + s(language, weighted_dl_distance, bs="re") +
      s(language, n_words, bs="re"), # RE
    data=crit_by_items %>% 
      filter(language != "en_exp")
  )
  
  save(
    fit_overall_full_withLength_binary, 
    fit_overall_full_withLength_gradient, 
    fit_overall_full_withLength_comm,
    file="models/fits_overall_full_withLength.RData"
  )
  
} else {
  load("models/fits_overall_full_withLength.RData")
}

## Evaluate Convergence
fit_overall_full_withLength_binary$conv
fit_overall_full_withLength_gradient$conv
fit_overall_full_withLength_comm$conv

## Extract the marginal log-likelihood of each model. scam's own logLik()/AIC()/
## BIC() are conditional on the fitted random effects, so they do not match the
## values reported in the SI; use marginal_logLik() throughout.
marginal_logLik(fit_overall_full_withLength_binary)
marginal_logLik(fit_overall_full_withLength_gradient)
marginal_logLik(fit_overall_full_withLength_comm)

## Save the figures
png("figures/mod_binary_nwords.png", width=6, height=6, units="in", res=300)
plot(fit_overall_full_withLength_binary, select = 1,
     shade = TRUE,            # optional: shaded CI instead of dashed lines
     shade.col = "lightblue", # color for shading
     main = "Binary Grammaticality", # title
     xlab = "Length",    # x-axis label
     ylab = "Acceptability (centered)", # y-axis label
     rug = TRUE,              # keep the rug
     jit = TRUE               # jitter the rug ticks
)
dev.off()

png("figures/mod_gradient_dl.png", width=6, height=6, units="in", res=300)
plot(fit_overall_full_withLength_gradient, select = 1,
     shade = TRUE,            # optional: shaded CI instead of dashed lines
     shade.col = "lightblue", # color for shading
     main = "Gradient Grammaticality", # title
     xlab = "DL distance",    # x-axis label
     ylab = "Acceptability (centered)", # y-axis label
     rug = TRUE,              # keep the rug
     jit = TRUE               # jitter the rug ticks
)
dev.off()

png("figures/mod_gradient_nwords.png", width=6, height=6, units="in", res=300)
plot(fit_overall_full_withLength_gradient, select = 2,
     shade = TRUE,            # optional: shaded CI instead of dashed lines
     shade.col = "lightblue", # color for shading
     main = "Gradient Grammaticality", # title
     xlab = "Length",    # x-axis label
     ylab = "Acceptability (centered)", # y-axis label
     rug = TRUE,              # keep the rug
     jit = TRUE               # jitter the rug ticks
)
dev.off()

png("figures/mod_communication_dlweighted.png", width=6, height=6, units="in", res=300)
plot(fit_overall_full_withLength_comm, select = 1,
     shade = TRUE,            # optional: shaded CI instead of dashed lines
     shade.col = "lightblue", # color for shading
     main = "Communication", # title
     xlab = "DL distance / length",    # x-axis label
     ylab = "Acceptability (centered)", # y-axis label
     rug = TRUE,              # keep the rug
     jit = TRUE               # jitter the rug ticks
)
dev.off()

png("figures/mod_communication_nwords.png", width=6, height=6, units="in", res=300)
plot(fit_overall_full_withLength_comm, select = 2,
     shade = TRUE,            # optional: shaded CI instead of dashed lines
     shade.col = "lightblue", # color for shading
     main = "Communication", # title
     xlab = "Length",    # x-axis label
     ylab = "Acceptability (centered)", # y-axis label
     rug = TRUE,              # keep the rug
     jit = TRUE               # jitter the rug ticks
)
dev.off()

main_likelihoods <- tibble(
  model = c("binary", "gradient", "comm"),
  log_likelihood = c(
    marginal_logLik(fit_overall_full_withLength_binary),
    marginal_logLik(fit_overall_full_withLength_gradient),
    marginal_logLik(fit_overall_full_withLength_comm)
  )
) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood,
    model = factor(model, levels=c("binary", "gradient", "comm"))
  )

p_main_likelihoods <- ggplot(main_likelihoods) +
  geom_point(
    aes(
      y="All Languages",
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

p_main_likelihoods

### ---- Analysis By Language ----
if (file.exists("models/fits_byLang_full_withLength.RData")) {
  load("models/fits_byLang_full_withLength.RData")
} else {
  
  ## Fit a model for each language
  ## NB: this in-loop tibble records the CONDITIONAL log-likelihood (scam's own
  ## logLik(), which plugs in the fitted random effects). It is kept only as a
  ## diagnostic; the figures below use the MARGINAL values rebuilt with
  ## marginal_stats_bylang() after this cache block.
  fits_byLang_full_withLength_stats <- tibble(
    language = character(),
    model = character(),
    log_likelihood_conditional = double()
  )
  
  fits_byLang_full_withLength <- list()
  for (language_iter in levels(crit_by_items$language)) {
    print(
      paste0("Fitting models for ", language_iter)
    )
    ## Get specific data
    crit_by_items_language_specific <- crit_by_items %>%
      filter(language == language_iter)
    
    if (language_iter == "en_exp") {
      ## Fit models
      binary <- scam(
        response_num ~ grammatical_numeric + s(n_words, bs="ps") +
          s(item, bs="re") + s(item, grammatical_numeric, bs="re") + s(item, n_words, bs="re"),
        data=crit_by_items_language_specific
      )
      
      gradient <- scam(
        response_num ~ s(dl_distance, bs="mpd") + s(n_words, bs="ps") + 
          s(item, bs="re") + s(item, dl_distance, bs="re") + s(item, n_words, bs="re"),
        data=crit_by_items_language_specific
      )
      
      comm <- scam(
        response_num ~ s(weighted_dl_distance, bs="mpd") + s(n_words, bs="mpd") +
          s(item, bs="re") + s(item, weighted_dl_distance, bs="re") + + s(item, n_words, bs="re"),
        data=crit_by_items_language_specific
      )
    } else {
      ## Fit models
      binary <- scam(
        response_num ~ grammatical_numeric + s(n_words, bs="ps") +
          s(item, bs="re") ,
        data=crit_by_items_language_specific
      )
      
      gradient <- scam(
        response_num ~ s(dl_distance, bs="mpd") + s(n_words, bs="ps") +
          s(item, bs="re"),
        data=crit_by_items_language_specific
      )
      
      comm <- scam(
        response_num ~ s(weighted_dl_distance, bs="mpd") + s(n_words, bs="mpd") +
          s(item, bs="re"),
        data=crit_by_items_language_specific
      )
    }
    
    
    
    fits_byLang_full_withLength[[language_iter]] <- list(
      binary = binary,
      gradient = gradient,
      comm = comm
    )
    
    fits_byLang_full_withLength_stats <- fits_byLang_full_withLength_stats %>%
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
  
  fits_byLang_full_withLength_stats <- fits_byLang_full_withLength_stats %>%
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
    fits_byLang_full_withLength,
    fits_byLang_full_withLength_stats,
    file="models/fits_byLang_full_withLength.RData"
  )
}

### ---- Figure based on individual-language likelihoods  
fits_byLang_full_withLength_stats <- marginal_stats_bylang(
  fits_byLang_full_withLength,
  c(binary = "Binary\nGrammaticality", gradient = "Gradient\nGrammaticality", comm = "Communication-based"))

fits_byLang_full_withLength_stats <- fits_byLang_full_withLength_stats %>%
  group_by(language) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood
  )

p_language_likelihoods <- ggplot(fits_byLang_full_withLength_stats) +
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

## Stack the combined "All Languages" panel (top) over the by-language panel
## (bottom), matching the main-text likelihood figure.
likelihood_differences_with_length_figure <- p_main_likelihoods / p_language_likelihoods +
  plot_layout(heights = c(1, 3))

ggsave("figures/likelihood_differences_with_length.png",
       plot = likelihood_differences_with_length_figure, width = 8, height = 8)

## ---- Model Predictions Figure ----

### ---- Baseline scatterplot (item means) ----
p_baseline <- ggplot(
  plot_data_by_item,
  aes(x = n_words, y = response, color = dl_distance_bucket)
) +
  geom_point(alpha = 0.3) +
  facet_wrap(~language, ncol = 1, strip.position = "left") +
  labs(x="Sentence Length in Words", y="Rating") +
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
    aes(x = n_words, y = response, color = dl_distance_bucket),
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

ggsave(p_gam, filename = "figures/main_figure_gam.png",
       width    = 6,
       height   = 16,
       units    = "in",
       dpi      = 300)

### ---- Model-based plots ----
#### ---- Grammatical model ----
#### ---- Predictions from grammatical model ----
grammatical_predictions <- setNames(
  lapply(names(fits_byLang_full_withLength), function(lang) {
    predict(
      fits_byLang_full_withLength[[lang]]$binary,
      newdata = get_prediction_grid(fits_byLang_full_withLength, lang, "binary", n_words = range(crit[crit$language == lang,]$n_words)),
      type = "response",
      se.fit = TRUE,
      exclude = "s(item)"
    )
  }
  ),
  names(fits_byLang_full_withLength)
)

## Combine list to tibble
grammatical_plot_data <- bind_rows(
  lapply(names(grammatical_predictions), function(lang) {
    preds <- grammatical_predictions[[lang]]
    grid  <- get_prediction_grid(fits_byLang_full_withLength, lang, "binary", n_words = range(crit[crit$language == lang,]$n_words))
    
    tibble(
      language = label_language(lang),
      n_words = grid$n_words,
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
      x = n_words,
      ymin = response - 2 * se,
      ymax = response + 2 * se,
      group = dl_distance_bucket
    ),
    alpha = 0.25
  ) +
  geom_line(
    data=grammatical_plot_data,
    aes(
      x = n_words,
      y = response,
      color = dl_distance_bucket
    ),
    linewidth = 1.5
  ) +
  facet_wrap(~language, ncol = 1) +
  labs(
    x = "Sentence Length in Words",
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

p_grammatical

ggsave(
  p_grammatical,
  filename = "figures/withLength_figure_grammatical.png",
  width    = 6,
  height   = 16,
  units    = "in",
  dpi      = 300
)

#### ---- Distance model ----
#### ---- Predictions from distance model ----
distance_predictions <- setNames(
  lapply(names(fits_byLang_full_withLength), function(lang) {
    predict(
      fits_byLang_full_withLength[[lang]]$gradient,
      newdata = get_prediction_grid(fits_byLang_full_withLength, lang, "gradient", n_words = range(crit[crit$language == lang,]$n_words)),
      type = "response",
      se.fit = TRUE,
      exclude = "s(item)"
    )
  }
  ),
  names(fits_byLang_full_withLength)
)

## Combine list to tibble
distance_plot_data <- bind_rows(
  lapply(names(distance_predictions), function(lang) {
    preds <- distance_predictions[[lang]]
    grid  <- get_prediction_grid(fits_byLang_full_withLength, lang, "gradient", n_words = range(crit[crit$language == lang,]$n_words))
    
    tibble(
      language = label_language(lang),
      n_words = grid$n_words,
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
      x = n_words,
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
      x = n_words,
      y = response,
      color = dl_distance_bucket
    ),
    linewidth = 1.5
  ) +
  facet_wrap(~language, ncol = 1) +
  labs(
    x = "Sentence Length in Words",
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
  filename = "figures/withLength_figure_distance.png",
  width    = 6,
  height   = 16,
  units    = "in",
  dpi      = 300
)

#### ---- Weighted distance model ----
#### ---- Predictions from weighted distance model ----
weighted_distance_predictions <- setNames(
  lapply(names(fits_byLang_full_withLength), function(lang) {
    predict(
      fits_byLang_full_withLength[[lang]]$comm,
      newdata = get_prediction_grid(fits_byLang_full_withLength, lang, "comm", n_words = range(crit[crit$language == lang,]$n_words)),
      type = "response",
      se.fit = TRUE,
      exclude = "s(item)"
    )
  }
  ),
  names(fits_byLang_full_withLength)
)

## Combine list to tibble
weighted_distance_plot_data <- bind_rows(
  lapply(names(weighted_distance_predictions), function(lang) {
    preds <- weighted_distance_predictions[[lang]]
    grid  <- get_prediction_grid(fits_byLang_full_withLength, lang, "comm", n_words = range(crit[crit$language == lang,]$n_words))
    
    tibble(
      language = label_language(lang),
      n_words = grid$n_words,
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
      x = n_words,
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
      x = n_words,
      y = response,
      color = dl_distance_bucket
    ),
    linewidth = 1.5
  ) +
  facet_wrap(~language, ncol = 1) +
  labs(
    x = "Sentence Length in Words",
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
  filename = "figures/withLength_figure_weighted_distance.png",
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
  filename = "figures/withLength_figure.png",
  plot     = main_figure,
  width    = 25,
  height   = 20,
  units    = "in",
  dpi      = 300
)

## Save as pdf
ggsave(
  filename = "figures/withLength_figure.pdf",
  plot     = main_figure,
  width    = 25,
  height   = 20,
  units    = "in"
)
