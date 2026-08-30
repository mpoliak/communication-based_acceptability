## =============================================================================
## SI Section 3: Quantifying information using surprisal instead of length
##
## Replaces utterance length with summed GPT surprisal throughout.
## Produces with_surprisal_figure.png and likelihood_differences_with_surprisal.png.
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

## ---- Main Analysis with Surprisal ----
## Operationalizes information content as summed token surprisal from monolingual
## Goldfish language models instead of utterance length in tokens.
## Produces with_surprisal_figure.png and likelihood_differences_with_surprisal.png.
### ---- Linear Model ----
lmer_fits <- lapply(unique(crit$language), function(current_language) {
  print(paste0("Processing language: ", current_language))
  lmer(
    response_num ~ dl_distance * sum_surprisal + 
      (dl_distance * sum_surprisal|participant_id) + 
      (dl_distance|item),
    data=crit %>% filter(
      language == current_language,
      dl_distance <= 3
    ) %>%
      mutate(
        sum_surprisal = scale(sum_surprisal)[,1],
        dl_distance = scale(dl_distance)[,1]
      ),
    control = lmerControl(
      optimizer ="bobyqa",
      optCtrl = list(maxfun = 1e7)
    )
  )
})

lapply(lmer_fits, function(model_fit) {
  model_fit@optinfo$conv
})

names(lmer_fits) <- unique(crit$language)

## Fit full RE for en_exp. The item id encodes the length version
## (en_experimental_62m); strip the trailing s/m/l to recover the 72 base items,
## each of which carries all 3 lengths x 4 grammaticality conditions. Length then
## varies within item, making the by-item surprisal slope identifiable.
lmer_fits[["en_exp"]] <- lmer(
  response_num ~ dl_distance * sum_surprisal +
    (dl_distance * sum_surprisal|participant_id) +
    (dl_distance * sum_surprisal|item),
  data=crit %>% filter(
    language == "en_exp",
    dl_distance <= 3
  ) %>%
    mutate(
      item = as.factor(gsub("[sml]$", "", item)),
      sum_surprisal = scale(sum_surprisal)[,1],
      dl_distance = scale(dl_distance)[,1]
    ),
  control = lmerControl(
    optimizer ="bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)

lapply(lmer_fits, summary)

## Get largest p-value from dl_distance
lapply(lmer_fits, function(model_fit) {
  summary(model_fit)$coefficients[2,5]
}) %>% unlist %>% max

## Get largest p-value from interaction
lapply(lmer_fits, function(model_fit) {
  summary(model_fit)$coefficients[4,5]
}) %>% unlist %>% max

### ---- Effects of Length in Well-formed Materials ----
well_formed <- crit %>% filter(dl_distance == 0)

## sum_surprisal is a between-item predictor, so the by-item random intercept is
## needed for its standard error to be correct. For these 8 languages surprisal is
## constant within an item, so only an intercept (not a slope) is estimable by item.
length_effects <- lapply(unique(well_formed$language), function(current_language) {
  print(paste0("Processing language: ", current_language))
  lmer(
    response_num ~ sum_surprisal + (sum_surprisal | participant_id) + (1 | item),
    data=well_formed %>% filter(
      language == current_language
    ) %>%
      mutate(
        sum_surprisal = scale(sum_surprisal)[,1]
      ),
    control = lmerControl(
      optimizer ="bobyqa",
      optCtrl = list(maxfun = 1e5)
    )
  )
})

names(length_effects) <- unique(well_formed$language)

## en_exp varied length within item, so a by-item surprisal slope is estimable
## once the item id is collapsed to the 72 base items.
length_effects[["en_exp"]] <- lmer(
  response_num ~ sum_surprisal +
    (sum_surprisal | participant_id) +
    (sum_surprisal | item),
  data=well_formed %>% filter(
    language == "en_exp"
  ) %>%
    mutate(
      item = as.factor(gsub("[sml]$", "", item)),
      sum_surprisal = scale(sum_surprisal)[,1]
    ),
  control = lmerControl(
    optimizer ="bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
)

lapply(length_effects, summary)

lapply(names(length_effects), function(current_language) {
  paste0(current_language,
         ifelse(
           summary(length_effects[[current_language]])$coefficients[2,1] > 0,
           " + ",
           " - "
         ),
         ifelse(
           summary(length_effects[[current_language]])$coefficients[2,5] < .05,
           "<.05",
           ">.05"
         )
         
  ) %>%
    print()
})


### ---- Surprisal Figure ----
unique_sentences <- crit_by_items %>%
  filter(condition == "original") %>%
  select(item, n_words, sum_surprisal) %>%
    distinct()

cor(unique_sentences$n_words, unique_sentences$sum_surprisal)

plot_data <- crit %>%
  mutate(
    dl_distance_bucket = cut(
      dl_distance, 
      breaks=c(-1, 0, 2, 5, Inf),
      labels=c("0", "1-2", "3-5", "6+"))
  )

table(plot_data$dl_distance_bucket, plot_data$dl_distance)

## Give full names to languages
plot_data <- plot_data %>%
  mutate(
    language = case_when(
      language == "da" ~ "Danish",
      language == "de" ~ "German",
      language == "en" ~ "English",
      language == "hi" ~ "Hindi",
      language == "fr" ~ "French",
      language == "ko" ~ "Korean",
      language == "ru" ~ "Russian",
      language == "zh" ~ "Mandarin",
      language == "en_exp" ~ "English\n(controlled)"
    ),
    language = factor(
      language,
      levels=c("Danish", "English", "French", "German", "Hindi", "Korean", 
               "Mandarin", "Russian", "English\n(controlled)")
    )
  )

plot_data_by_item <- plot_data %>%
  group_by(language, item, sum_surprisal, dl_distance_bucket) %>%
  summarize(
    response = mean(response_num)
  )

main_plot_surprisal <- ggplot(
  plot_data_by_item,
  aes(x=sum_surprisal, y=response, color=dl_distance_bucket)
) +
  geom_point(
    data=plot_data_by_item,
    alpha=.3
  ) +
  geom_smooth(
    se=T,
    method="gam",
    formula = y ~ s(x, bs = "tp")
  ) +
  facet_wrap(~language) +
  labs(
    x="Sum Surprisal",
    y="Rating",
    color="Damerau-Levenshtein Distance"
  ) +
  theme(legend.position="top") +
  scale_y_continuous(
    breaks = seq(1, 7, 2),
    labels = seq(1, 7, 2)
  ) +
  scale_color_manual(
    values=c("#EB7AA3", "#1E88E5", "#FFC107", "#004D40")
  )

main_plot_surprisal
ggsave("figures/all_languages_surprisal.png", width=12, height=8)  

### ---- Main Models ----
if (!file.exists("models/full_scam_models_surprisal.RData")) {
  mod_grammatical <- scam(
    response_num ~ grammatical_numeric + # + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, grammatical_numeric, bs="re"), # RE
    data=crit_by_items %>% 
      filter(language != "en_exp"),
    control = scam.control(trace = TRUE)
  )
  
  mod_distance <- scam(
    response_num ~ s(dl_distance, bs="mpd") + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, dl_distance, bs="re"), # RE
    data=crit_by_items %>% 
      filter(language != "en_exp"),
    control = scam.control(trace = TRUE)
  )
  
  mod_weighted_distance <- scam(
    response_num ~ s(weighted_dl_distance_surprisal, bs="mpd") + s(sum_surprisal, bs="mpd") + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, weighted_dl_distance_surprisal, bs="re") + 
      s(language, sum_surprisal, bs="re"), # RE
    data=crit_by_items %>% 
      filter(language != "en_exp"),
    control = scam.control(trace = TRUE)
  )
  
  save(
    mod_grammatical, 
    mod_distance, 
    mod_weighted_distance,
    file="models/full_scam_models_surprisal.RData"
  )
  
} else {
  load("models/full_scam_models_surprisal.RData")
}

## Evaluate Convergence
mod_grammatical$conv
mod_distance$conv
mod_weighted_distance$conv

## Extract the marginal log-likelihood of each model (see the note above: scam's
## own logLik()/AIC()/BIC() are conditional on the fitted random effects).
marginal_logLik(mod_grammatical)
marginal_logLik(mod_distance)
marginal_logLik(mod_weighted_distance)

## Save the figures
png("figures/mod_gradient.png", width=6, height=6, units="in", res=300)
plot(mod_distance, select = 1,
     shade = TRUE,            # optional: shaded CI instead of dashed lines
     shade.col = "lightblue", # color for shading
     main = "Gradient Grammaticality", # title
     xlab = "DL distance",    # x-axis label
     ylab = "Acceptability (centered)", # y-axis label
     rug = TRUE,              # keep the rug
     jit = TRUE               # jitter the rug ticks
)
dev.off()

png("figures/mod_communication_dlweighted_surprisal.png", width=6, height=6, units="in", res=300)
plot(mod_weighted_distance, select = 1,
     shade = TRUE,            # optional: shaded CI instead of dashed lines
     shade.col = "lightblue", # color for shading
     main = "Communication", # title
     xlab = "DL distance / sum surprisal",    # x-axis label
     ylab = "Acceptability (centered)", # y-axis label
     rug = TRUE,              # keep the rug
     jit = TRUE               # jitter the rug ticks
)
dev.off()

png("figures/mod_communication_sum_surprisal.png", width=6, height=6, units="in", res=300)
plot(mod_weighted_distance, select = 2,
     shade = TRUE,            # optional: shaded CI instead of dashed lines
     shade.col = "lightblue", # color for shading
     main = "Communication", # title
     xlab = "Sum surprisal",    # x-axis label
     ylab = "Acceptability (centered)", # y-axis label
     rug = TRUE,              # keep the rug
     jit = TRUE               # jitter the rug ticks
)
dev.off()

main_likelihoods <- tibble(
  model = c("grammatical", "distance", "weighted_distance"),
  log_likelihood = c(
    marginal_logLik(mod_grammatical),
    marginal_logLik(mod_distance),
    marginal_logLik(mod_weighted_distance)
  )
) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood,
    model = factor(model, levels=c("grammatical", "distance", "weighted_distance"))
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


## ---- Analysis By Language ----
### ---- Analysis By Language ----
if (file.exists("models/fits_byLang_full_with_surprisal.RData")) {
  load("models/fits_byLang_full_with_surprisal.RData")
} else {
  
  ## Fit a model for each language
  ## Conditional log-likelihoods, diagnostic only -- see the note in Section 1.
  ## This rebuild matters most here: under the CONDITIONAL likelihood the
  ## gradient model beats the communication model in Danish and German, which is
  ## reversed once the item random effects are integrated out.
  fits_byLang_full_with_surprisal_stats <- tibble(
    language = character(),
    model = character(),
    log_likelihood_conditional = double()
  )
  
  fits_byLang_full_with_surprisal <- list()
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
          s(item, bs="re") + s(item, weighted_dl_distance_surprisal, bs="re") + + s(item, sum_surprisal, bs="re"),
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
    
    fits_byLang_full_with_surprisal[[language_iter]] <- list(
      binary = binary,
      gradient = gradient,
      comm = comm
    )
    
    fits_byLang_full_with_surprisal_stats <- fits_byLang_full_with_surprisal_stats %>%
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
  
  fits_byLang_full_with_surprisal_stats <- fits_byLang_full_with_surprisal_stats %>%
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
    fits_byLang_full_with_surprisal,
    fits_byLang_full_with_surprisal_stats,
    file="models/fits_byLang_full_with_surprisal.RData"
  )
}

### ---- Figure based on individual-language likelihoods  
fits_byLang_full_with_surprisal_stats <- marginal_stats_bylang(
  fits_byLang_full_with_surprisal,
  c(binary = "Binary\nGrammaticality", gradient = "Gradient\nGrammaticality", comm = "Communication-based"))

fits_byLang_full_with_surprisal_stats <- fits_byLang_full_with_surprisal_stats %>%
  group_by(language) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood
  )

p_language_likelihoods <- ggplot(fits_byLang_full_with_surprisal_stats) +
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
likelihood_differences_with_surprisal_figure <- p_main_likelihoods / p_language_likelihoods +
  plot_layout(heights = c(1, 3))

ggsave("figures/likelihood_differences_with_surprisal.png",
       plot = likelihood_differences_with_surprisal_figure, width = 8, height = 8)

## ---- Model Predictions Figure ----
#### ---- Helper Functions 
YLIM <- c(0, 8)
DL_COLORS <- c("0"="#EB7AA3", "1-2"="#1E88E5", "3"="#FFC107", "4+"="#BBBBBB")

get_prediction_grid <- function(
    fits,
    lang,
    model_type,
    sum_surprisal=c(4, 40), 
    dl_distance=c(0, 1, 3)
) {
  expand.grid(
    sum_surprisal = seq(sum_surprisal[1] - .25, sum_surprisal[2] + .25, 1),
    dl_distance = dl_distance
  ) %>%
    mutate(
      weighted_dl_distance_surprisal = dl_distance / sum_surprisal,
      grammatical_numeric = ifelse(dl_distance == 0, .5, -.5),
      item=factor(fits[[lang]][[model_type]]$model$item[1], levels=levels(fits[[lang]][[model_type]]$model$item))
    )
}

### ---- Baseline data for plotting (item means) ----
plot_data_by_item <- crit %>% 
  mutate(
    weighted_dl_distance = crit$dl_distance / crit$sum_surprisal,
    grammatical = factor(
      ifelse(dl_distance == 0, "grammatical", "ungrammatical"),
      levels=c("ungrammatical", "grammatical")),
    item = as.factor(item),
    language = as.factor(language),
    response = as.factor(response),
    response_num = as.numeric(as.character(response)),
    participant_id = as.factor(participant_id),
    item = as.factor(item),
    language = as.factor(language)
  ) %>%
  group_by(language, item, sum_surprisal, dl_distance, weighted_dl_distance_surprisal, 
           grammatical, condition) %>%
  summarize(
    response_num = mean(response_num)
  ) %>%
  mutate(
    grammatical_numeric = ifelse(grammatical == "grammatical", .5, -.5)
  ) %>%
  dplyr::group_by(language, item, sum_surprisal, dl_distance) %>%
  dplyr::summarise(response = mean(response_num), .groups = "drop") %>%
  dplyr::mutate(
    dl_distance_bucket = bucket_dl(dl_distance),
    language           = label_language(language)
  )

### ---- Baseline scatterplot (item means) ----
p_baseline <- ggplot(
  plot_data_by_item,
  aes(x = sum_surprisal, y = response, color = dl_distance_bucket)
) +
  geom_point(alpha = 0.3) +
  facet_wrap(~language, ncol = 1, strip.position = "left") +
  labs(x="Sentence surprisal (sum)", y="Rating") +
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

ggsave(p_gam, filename = "figures/main_figure_gam_surprisal.png",
       width    = 6,
       height   = 16,
       units    = "in",
       dpi      = 300)

### ---- Model-based plots ----
#### ---- Grammatical model ----
#### ---- Predictions from grammatical model ----
grammatical_predictions <- setNames(
  lapply(names(fits_byLang_full_with_surprisal), function(lang) {
    predict(
      fits_byLang_full_with_surprisal[[lang]]$binary,
      newdata = get_prediction_grid(fits_byLang_full_with_surprisal, lang, "binary", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal)),
      type = "response",
      se.fit = TRUE,
      exclude = "s(item)"
    )
  }
  ),
  names(fits_byLang_full_with_surprisal)
)

## Combine list to tibble
grammatical_plot_data <- bind_rows(
  lapply(names(grammatical_predictions), function(lang) {
    preds <- grammatical_predictions[[lang]]
    grid  <- get_prediction_grid(fits_byLang_full_with_surprisal, lang, "binary", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal))
    
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
    x = "Sentence surprisal (sum)",
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
  filename = "figures/with_surprisal_figure_grammatical.png",
  width    = 6,
  height   = 16,
  units    = "in",
  dpi      = 300
)

#### ---- Distance model ----
#### ---- Predictions from distance model ----
distance_predictions <- setNames(
  lapply(names(fits_byLang_full_with_surprisal), function(lang) {
    predict(
      fits_byLang_full_with_surprisal[[lang]]$gradient,
      newdata = get_prediction_grid(fits_byLang_full_with_surprisal, lang, "gradient", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal)),
      type = "response",
      se.fit = TRUE,
      exclude = "s(item)"
    )
  }
  ),
  names(fits_byLang_full_with_surprisal)
)

## Combine list to tibble
distance_plot_data <- bind_rows(
  lapply(names(distance_predictions), function(lang) {
    preds <- distance_predictions[[lang]]
    grid  <- get_prediction_grid(fits_byLang_full_with_surprisal, lang, "gradient", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal))
    
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
    x = "Sentence surprisal (sum)",
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
  filename = "figures/with_surprisal_figure_distance.png",
  width    = 6,
  height   = 16,
  units    = "in",
  dpi      = 300
)

#### ---- Weighted distance model ----
#### ---- Predictions from weighted distance model ----
weighted_distance_predictions <- setNames(
  lapply(names(fits_byLang_full_with_surprisal), function(lang) {
    predict(
      fits_byLang_full_with_surprisal[[lang]]$comm,
      newdata = get_prediction_grid(fits_byLang_full_with_surprisal, lang, "comm", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal)),
      type = "response",
      se.fit = TRUE,
      exclude = "s(item)"
    )
  }
  ),
  names(fits_byLang_full_with_surprisal)
)

## Combine list to tibble
weighted_distance_plot_data <- bind_rows(
  lapply(names(weighted_distance_predictions), function(lang) {
    preds <- weighted_distance_predictions[[lang]]
    grid  <- get_prediction_grid(fits_byLang_full_with_surprisal, lang, "comm", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal))
    
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
    x = "Sentence surprisal (sum)",
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
  filename = "figures/with_surprisal_figure_weighted_distance.png",
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
  filename = "figures/with_surprisal_figure.png",
  plot     = main_figure,
  width    = 25,
  height   = 20,
  units    = "in",
  dpi      = 300
)

## Save as pdf
ggsave(
  filename = "figures/with_surprisal_figure.pdf",
  plot     = main_figure,
  width    = 25,
  height   = 20,
  units    = "in"
)




### ---- Main Figure with Model Predictions (Surprisal) ----

# Build a prediction grid: 36 surprisal values x 4 DL distances, per language
model_predictions_surprisal <- crit_by_items %>% 
  ungroup() %>%
  select(language) %>%
  distinct() %>%
  dplyr::slice(rep(1:n(), each = 36*4)) %>%
  mutate(
    sum_surprisal = rep(
      seq(
        from = quantile(crit_by_items$sum_surprisal, 0, na.rm = TRUE),
        to   = quantile(crit_by_items$sum_surprisal, 1, na.rm = TRUE),
        length.out = 36
      ),
      times = n_distinct(crit_by_items$language) * 4
    ),
    dl_distance = rep(c(0, 1, 3, 20), times = 36 * n_distinct(crit_by_items$language)),
    # Need these because the model includes them:
    n_words = mean(crit_by_items$n_words, na.rm = TRUE),
    weighted_dl_distance_surprisal = dl_distance / sum_surprisal,
    ## The combined model carries s(item, bs="re"), so predict() needs an item
    ## column to build the model frame. Supply any valid level and exclude the
    ## term below, giving predictions at the population (item-average) level.
    item = factor(levels(mod_weighted_distance$model$item)[1],
                  levels = levels(mod_weighted_distance$model$item))
  ) %>%
  filter(language != "en_exp")

pred_surp <- predict(
  mod_weighted_distance,
  newdata = model_predictions_surprisal,
  se.fit = TRUE,
  exclude = "s(item)"
)

model_predictions_surprisal <- model_predictions_surprisal %>%
  mutate(
    response = pred_surp$fit,
    se = pred_surp$se.fit,
    lower = response - 1.96 * se,
    upper = response + 1.96 * se
  )

# English controlled model (fit separately)
en_exp_predictions_surprisal <- tibble(
  sum_surprisal = rep(
    seq(
      from = quantile(crit_by_items$sum_surprisal, 0, na.rm = TRUE),
      to   = quantile(crit_by_items$sum_surprisal, 1, na.rm = TRUE),
      length.out = 36
    ),
    times = 4
  ),
  dl_distance = rep(c(0, 1, 3, 20), times = 36),
  n_words = mean(crit_by_items$n_words, na.rm = TRUE),
  weighted_dl_distance_surprisal = dl_distance / sum_surprisal,
  language = as.factor("en_exp"),
  item = factor(
    fits_byLang_full_with_surprisal[["en_exp"]][["comm"]]$model$item[1],
    levels = levels(fits_byLang_full_with_surprisal[["en_exp"]][["comm"]]$model$item)
  )
)

pred_en_exp_surp <- predict(
  fits_byLang_full_with_surprisal[["en_exp"]][["comm"]],
  newdata = en_exp_predictions_surprisal,
  se.fit = TRUE,
  exclude = "s(item)"
)

en_exp_predictions_surprisal <- en_exp_predictions_surprisal %>%
  mutate(
    response = pred_en_exp_surp$fit,
    se = pred_en_exp_surp$se.fit,
    lower = response - 1.96 * se,
    upper = response + 1.96 * se
  )

model_predictions_surprisal <- bind_rows(
  model_predictions_surprisal,
  en_exp_predictions_surprisal
) %>%
  mutate(
    language = case_when(
      language == "da" ~ "Danish",
      language == "de" ~ "German",
      language == "en" ~ "English",
      language == "hi" ~ "Hindi",
      language == "fr" ~ "French",
      language == "ko" ~ "Korean",
      language == "ru" ~ "Russian",
      language == "zh" ~ "Mandarin",
      language == "en_exp" ~ "English\n(controlled)"
    ),
    language = factor(
      language,
      levels=c("Danish", "English", "French", "German", "Hindi", "Korean", 
               "Mandarin", "Russian", "English\n(controlled)")
    ),
    dl_distance_bucket = bucket_dl(dl_distance)
  )

ggplot(
  model_predictions_surprisal,
  aes(x = sum_surprisal, y = response, color = dl_distance_bucket)
) +
  geom_point(data = plot_data_by_item, alpha = .3) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper, group = dl_distance_bucket),
    alpha = .3,
    color = NA
  ) +
  geom_line(linewidth = 1.5) +
  facet_wrap(~language) +
  labs(
    x = "Sum Surprisal",
    y = "Rating",
    color = "Damerau-Levenshtein Distance"
  ) +
  theme(legend.position="top") +
  scale_y_continuous(
    breaks = seq(1, 7, 2),
    labels = seq(1, 7, 2)
  ) +
  scale_color_manual(values=c("#EB7AA3", "#1E88E5", "#FFC107", "#004D40"))

ggsave("figures/all_languages_with_model_predictions_surprisal.png", width=12, height=8)
