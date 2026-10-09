## =============================================================================
## SI Section 3: Quantifying information using length in tokens instead of
## surprisal
##
## The main text measures a sentence's information as its summed surprisal. This
## section repeats the main analyses with information measured as sentence
## length in words instead, which is how the submitted paper measured it:
##   - mixed-effects models of DL distance x length, per experiment
##   - the effect of length on well-formed sentences, tested across experiments
##   - the five-model comparison (Binary, Gradient, Communication, SLOR, MORCELA),
##     per experiment and pooled over all nine experiments, with the
##     Communication model's information term on length
##   - the curve figure: ratings and model predictions over length
## Produces likelihood_differences_length.png, curves_length.png/.pdf and
## length_mod_*.png, and writes its statistics to outputs/s3_length_*.csv.
## =============================================================================

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

## ---- Linear Mixed-Effects Models ----
## The submitted paper's main-text regression: ratings on DL distance, length and
## their interaction, fit separately for each experiment with both predictors
## z-scored within the experiment. en_exp is fit on its own below, with a fuller
## random-effects structure, so it is left out of this loop.
### ---- Linear Model ----
lmer_languages <- setdiff(as.character(unique(crit$language)), "en_exp")

lmer_fits <- lapply(lmer_languages, function(current_language) {
  print(paste0("Processing language: ", current_language))
  lmer(
  response_num ~ dl_distance * n_words +
    (dl_distance * n_words|participant_id) +
    (dl_distance|item),
  data=crit %>% filter(
    language == current_language,
    dl_distance <= 3
    ) %>%
    mutate(
      n_words = scale(n_words)[,1],
      dl_distance = scale(dl_distance)[,1]
    ) %>%
    droplevels(),
  control = lmerControl(
    optimizer ="bobyqa",
    optCtrl = list(maxfun = 1e5)
  )
  )
})

names(lmer_fits) <- lmer_languages

## Fit full RE for en_exp.
## In this experiment each of the 72 base items was written in 3 lengths (short,
## medium, long), and each of those crossed with the 4 grammaticality conditions.
## The item id in the data encodes the length version (e.g. en_experimental_62m),
## so we strip the trailing s/m/l to recover the 72 base items. Length then varies
## WITHIN item, which is what makes the by-item n_words slope identifiable.
lmer_fits[["en_exp"]] <- lmer(
  response_num ~ dl_distance * n_words +
    (dl_distance * n_words|participant_id) +
    (dl_distance * n_words|item),
  data=crit %>% filter(
    language == "en_exp",
    dl_distance <= 3
  ) %>%
    mutate(
      item = as.factor(gsub("[sml]$", "", item)),
      n_words = scale(n_words)[,1],
      dl_distance = scale(dl_distance)[,1]
    ) %>%
    droplevels(),
  ## nlminb rather than bobyqa, as in analysis/main_analysis.R: with this
  ## random-effects structure bobyqa stops ~28 log-likelihood units short of the
  ## optimum and fails the convergence check. The fit is still singular.
  control = lmerControl(
    optimizer ="nlminbwrap",
    optCtrl = list(eval.max = 1e5, iter.max = 1e5)
  )
)

lapply(lmer_fits, summary)

## Get largest p-value
## DL Distance
lapply(lmer_fits, function(model_fit) {
  summary(model_fit)$coefficients[2,5]
}) %>% unlist %>% max

## Length
lapply(lmer_fits, function(model_fit) {
  summary(model_fit)$coefficients[3,5]
}) %>% unlist %>% max

## Interaction
lapply(lmer_fits, function(model_fit) {
  summary(model_fit)$coefficients[4,5]
}) %>% unlist %>% max

## Save the fixed effects of every model, with lme4's convergence messages (empty
## when the fit converged cleanly)
lmm_table <- bind_rows(lapply(names(lmer_fits), function(current_language) {
  model_fit <- lmer_fits[[current_language]]
  coefs <- summary(model_fit)$coefficients
  tibble(
    language = current_language,
    term = rownames(coefs),
    estimate = coefs[, "Estimate"],
    std_error = coefs[, "Std. Error"],
    df = coefs[, "df"],
    t_value = coefs[, "t value"],
    p_value = coefs[, "Pr(>|t|)"],
    singular = isSingular(model_fit),
    convergence_messages = paste(model_fit@optinfo$conv$lme4$messages, collapse = "; ")
  )
}))

write_csv(lmm_table, "outputs/s3_length_lmm.csv")

### ---- Effects of Length in Well-formed Materials ----
well_formed <- crit %>% filter(dl_distance == 0)

## n_words is a between-item predictor, so the by-item random intercept is needed
## for its standard error to be correct. For these 8 languages length is constant
## within an item, so only an intercept (not a slope) can be estimated by item.
## Length is left in raw words here, as in the paper, so the slopes are in rating
## points per word. en_exp is fit on its own below.
length_languages <- setdiff(as.character(unique(well_formed$language)), "en_exp")

length_effects <- lapply(length_languages, function(current_language) {
  print(paste0("Processing language: ", current_language))
  lmer(
    response_num ~ n_words + (n_words | participant_id) + (1 | item),
    data=well_formed %>% filter(
      language == current_language
    ) %>%
      droplevels(),
    control = lmerControl(
      optimizer ="bobyqa",
      optCtrl = list(maxfun = 1e5)
    )
  )
})

names(length_effects) <- length_languages

## en_exp manipulated length within item, so here length varies within the 72 base
## items and a by-item length slope is estimable (see the note above on item ids).
length_effects[["en_exp"]] <- lmer(
  response_num ~ n_words +
    (n_words | participant_id) +
    (n_words | item),
  data=well_formed %>% filter(
    language == "en_exp"
  ) %>%
    mutate(
      item = as.factor(gsub("[sml]$", "", item))
    ) %>%
    droplevels(),
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

## Save the length slope of every experiment
length_effects_table <- bind_rows(lapply(names(length_effects), function(current_language) {
  model_fit <- length_effects[[current_language]]
  coefs <- summary(model_fit)$coefficients
  tibble(
    language = current_language,
    estimate = coefs["n_words", "Estimate"],
    std_error = coefs["n_words", "Std. Error"],
    df = coefs["n_words", "df"],
    t_value = coefs["n_words", "t value"],
    p_value = coefs["n_words", "Pr(>|t|)"],
    singular = isSingular(model_fit),
    convergence_messages = paste(model_fit@optinfo$conv$lme4$messages, collapse = "; ")
  )
}))

write_csv(length_effects_table, "outputs/s3_length_error_free_slopes.csv")

## Test whether the length slopes differ from zero across the experiments
length_betas <- sapply(length_effects,function(m) fixef(m)[2])
length_t_test <- t.test(length_betas)
length_t_test

write_csv(
  tibble(
    n_experiments = length(length_betas),
    mean_slope = unname(length_t_test$estimate),
    t_value = unname(length_t_test$statistic),
    df = unname(length_t_test$parameter),
    p_value = length_t_test$p.value,
    conf_low = length_t_test$conf.int[1],
    conf_high = length_t_test$conf.int[2]
  ),
  "outputs/s3_length_error_free_ttest.csv"
)

## ---- Five-Model Comparison ----
## Binary, Gradient and Communication are the three models from the submitted
## paper; SLOR and MORCELA are the two language-model baselines from the main
## text. Every model is fit exactly as in the earlier SLOR/MORCELA comparison
## (removed in the 2026-10-09 restructure), whose Communication model also took
## length as its information term, so the log-likelihoods here reproduce the ones
## it reported. Only the Communication
## model differs from the main text's: here it is a function of DL distance per
## word and of length, there of DL distance per bit and of surprisal.
##
## The item means come from setup.R (distances above 3 dropped, en_exp item ids
## already stripped of their length suffix, so the three length versions of each
## base item share one item random effect). Rebuilding the two factors gives them
## sorted levels, as in that earlier comparison.
length_model_data <- crit_by_items %>%
  ungroup() %>%
  mutate(
    item = as.factor(as.character(item)),
    language = as.factor(as.character(language))
  )

## One experiment's data. MORCELA's two weights are learned separately for each
## experiment, since there is no reason a single beta and gamma should suit all
## of them.
language_data_length <- function(language_iter) {
  language_specific <- length_model_data %>%
    filter(language == language_iter) %>% droplevels()

  language_specific$morcela <- morcela_score(
    language_specific, morcela_weights(language_specific)
  )

  language_specific
}

### ---- Analysis By Language ----
if (file.exists("models/s3_length_by_language_5models.RData")) {
  load("models/s3_length_by_language_5models.RData")
} else {

  ## Fit a model for each language. SLOR and MORCELA get a monotonically
  ## INCREASING smooth: a higher score means the sentence is more probable than
  ## its words alone would predict, and so should be judged more acceptable.
  models_by_language_length <- list()
  for (language_iter in levels(length_model_data$language)) {
    print(
      paste0("Fitting models for ", language_iter)
    )

    ## Get specific data
    crit_by_items_language_specific <- language_data_length(language_iter)

    if (language_iter == "en_exp") {
      ## Fit models
      grammatical <- scam(
        response_num ~ grammatical_numeric +
          s(item, bs="re") + s(item, grammatical_numeric, bs="re"),
        data=crit_by_items_language_specific
      )

      distance <- scam(
        response_num ~ s(dl_distance, bs="mpd") +
          s(item, bs="re") + s(item, dl_distance, bs="re"),
        data=crit_by_items_language_specific
      )

      weighted_distance <- scam(
        response_num ~ s(weighted_dl_distance, bs="mpd") + s(n_words, bs="mpd") +
          s(item, bs="re") + s(item, weighted_dl_distance, bs="re") + s(item, n_words, bs="re"),
        data=crit_by_items_language_specific
      )

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
      grammatical <- scam(
        response_num ~ grammatical_numeric +
          s(item, bs="re") ,
        data=crit_by_items_language_specific
      )

      distance <- scam(
        response_num ~ s(dl_distance, bs="mpd") +
          s(item, bs="re"),
        data=crit_by_items_language_specific
      )

      weighted_distance <- scam(
        response_num ~ s(weighted_dl_distance, bs="mpd") + s(n_words, bs="mpd") +
          s(item, bs="re"),
        data=crit_by_items_language_specific
      )

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

    models_by_language_length[[language_iter]] <- list(
      grammatical = grammatical,
      distance = distance,
      weighted_distance = weighted_distance,
      slor = slor,
      morcela = morcela
    )

  }

  save(
    models_by_language_length,
    file="models/s3_length_by_language_5models.RData"
  )
}

## Make sure the cache holds what this script fits: all nine experiments, five
## models each, and a Communication model whose information term is length.
smooth_labels <- function(m) vapply(m$smooth, function(s) s$label, character(1))

stopifnot(
  setequal(names(models_by_language_length), levels(length_model_data$language)),
  all(vapply(models_by_language_length, function(M) {
    identical(names(M), c("grammatical", "distance", "weighted_distance", "slor", "morcela")) &&
      all(c("s(weighted_dl_distance)", "s(n_words)") %in% smooth_labels(M$weighted_distance))
  }, logical(1)))
)

## Evaluate Convergence
lapply(models_by_language_length, function(M) sapply(M, function(m) m$conv))

### ---- Analysis Across All Experiments ----
## Each model is also fit once to all nine experiments at the same time, with a
## random effect for language on top of the one for item. The language random
## effects are really experiment random effects, with English and controlled
## English as two separate levels. MORCELA's two weights are estimated once
## here, on the pooled data: this is one model of all nine experiments, so it
## gets one beta and one gamma.
crit_by_items_all_experiments <- length_model_data

crit_by_items_all_experiments$morcela <- morcela_score(
  crit_by_items_all_experiments, morcela_weights(crit_by_items_all_experiments)
)

if (file.exists("models/s3_length_all_experiments_5models.RData")) {
  load("models/s3_length_all_experiments_5models.RData")
} else {

  print("Fitting models for all experiments together")

  mod_grammatical <- scam(
    response_num ~ grammatical_numeric +  # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, grammatical_numeric, bs="re"), # RE
    data=crit_by_items_all_experiments,
    control = scam.control(trace = TRUE)
  )

  mod_distance <- scam(
    response_num ~ s(dl_distance, bs="mpd") + # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, dl_distance, bs="re"), # RE
    data=crit_by_items_all_experiments,
    control = scam.control(trace = TRUE)
  )

  mod_weighted_distance <- scam(
    response_num ~ s(weighted_dl_distance, bs="mpd") + s(n_words, bs="mpd") + # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, weighted_dl_distance, bs="re") +
      s(language, n_words, bs="re"), # RE
    data=crit_by_items_all_experiments,
    control = scam.control(trace = TRUE)
  )

  mod_slor <- scam(
    response_num ~ s(slor, bs="mpi") + # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, slor, bs="re"), # RE
    data=crit_by_items_all_experiments,
    control = scam.control(trace = TRUE)
  )

  mod_morcela <- scam(
    response_num ~ s(morcela, bs="mpi") + # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, morcela, bs="re"), # RE
    data=crit_by_items_all_experiments,
    control = scam.control(trace = TRUE)
  )

  save(
    mod_grammatical,
    mod_distance,
    mod_weighted_distance,
    mod_slor,
    mod_morcela,
    file="models/s3_length_all_experiments_5models.RData"
  )

}

## Make sure the cache holds the pooled length-based fits of all nine experiments
stopifnot(
  nrow(mod_weighted_distance$model) == nrow(crit_by_items_all_experiments),
  nlevels(mod_weighted_distance$model$language) == 9,
  all(c("s(weighted_dl_distance)", "s(n_words)") %in% smooth_labels(mod_weighted_distance))
)

## Evaluate Convergence
mod_grammatical$conv
mod_distance$conv
mod_weighted_distance$conv
mod_slor$conv
mod_morcela$conv

## Save the figures
## The rug ticks are jittered at random; a fixed seed makes the figure reproducible
set.seed(2026)
png("figures/length_mod_gradient_dl.png", width=6, height=6, units="in", res=300)
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

## The rug ticks are jittered at random; a fixed seed makes the figure reproducible
set.seed(2026)
png("figures/length_mod_communication_dlweighted.png", width=6, height=6, units="in", res=300)
plot(mod_weighted_distance, select = 1,
     shade = TRUE,            # optional: shaded CI instead of dashed lines
     shade.col = "lightblue", # color for shading
     main = "Communication", # title
     xlab = "DL distance / length",    # x-axis label
     ylab = "Acceptability (centered)", # y-axis label
     rug = TRUE,              # keep the rug
     jit = TRUE               # jitter the rug ticks
)
dev.off()

## The rug ticks are jittered at random; a fixed seed makes the figure reproducible
set.seed(2026)
png("figures/length_mod_communication_nwords.png", width=6, height=6, units="in", res=300)
plot(mod_weighted_distance, select = 2,
         shade = TRUE,            # optional: shaded CI instead of dashed lines
         shade.col = "lightblue", # color for shading
         main = "Communication", # title
         xlab = "Length",    # x-axis label
         ylab = "Acceptability (centered)", # y-axis label
         rug = TRUE,              # keep the rug
         jit = TRUE               # jitter the rug ticks
)
dev.off()

### ---- Marginal log-likelihoods ----
## Build the by-language log-likelihood table using the MARGINAL log-likelihood
## (integrating the item random effects out) rather than the conditional value
## logLik() returns. Done here, after the cache block, so it applies whether the
## models were just fit or loaded from disk.
stats_models_by_language <- do.call(rbind, lapply(names(models_by_language_length), function(lg) {
  M <- models_by_language_length[[lg]]
  tibble(language = lg, model = names(M),
         log_likelihood = vapply(M, marginal_logLik, numeric(1)),
         converged = vapply(M, function(m) isTRUE(m$conv), logical(1)))
})) %>%
  mutate(
    language_code = language,
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
      levels=c("English\n(controlled)", "Russian", "Mandarin",
               "Korean", "Hindi", "German", "French", "English", "Danish")
    ),
    model = case_when(
      model == "grammatical" ~ "Binary",
      model == "distance" ~ "Gradient",
      model == "weighted_distance" ~ "Communication",
      model == "slor" ~ "SLOR",
      model == "morcela" ~ "MORCELA"
    ) %>%
      factor(levels=MODEL_LEVELS)
  ) %>%
  group_by(language) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood
  ) %>%
  ungroup()

write_csv(stats_models_by_language, "outputs/s3_length_likelihoods_by_language.csv")

main_likelihoods <- tibble(
  model = c("Binary", "Gradient", "Communication", "SLOR", "MORCELA"),
  log_likelihood = c(
    marginal_logLik(mod_grammatical),
    marginal_logLik(mod_distance),
    marginal_logLik(mod_weighted_distance),
    marginal_logLik(mod_slor),
    marginal_logLik(mod_morcela)
  ),
  converged = c(
    isTRUE(mod_grammatical$conv),
    isTRUE(mod_distance$conv),
    isTRUE(mod_weighted_distance$conv),
    isTRUE(mod_slor$conv),
    isTRUE(mod_morcela$conv)
  )
) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood,
    model = factor(model, levels=MODEL_LEVELS)
  )

write_csv(main_likelihoods, "outputs/s3_length_likelihoods_all_experiments.csv")

### ---- Likelihood differences figure ----
## Drawn exactly like the main text's five-model figure, with the pooled "All
## Experiments" row on top of the by-experiment rows.
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

p_language_likelihoods <- ggplot(stats_models_by_language) +
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
  scale_y_discrete(name="")

## Stack the combined "All Experiments" panel (top) over the by-language panel
## (bottom).
likelihood_differences_figure <- p_main_likelihoods / p_language_likelihoods +
  plot_layout(heights = c(1, 3))

ggsave("figures/likelihood_differences_length.png",
       plot = likelihood_differences_figure, width = 8, height = 8)

## ---- Curve Figure ----
## The submitted paper's main figure: item-mean ratings over length, with the
## Binary, Gradient and Communication models' predictions drawn as curves for
## DL distances 0, 1 and 3. The predictions come from the by-experiment fits
## above. YLIM, the DL colours and buckets, label_language() and
## scale_dl_color() come from setup.R.
#### ---- Helper Functions
## The level of `item` is only a placeholder so predict() sees every term in the
## model; all of the item random effects are excluded from the predictions
## themselves (random_effect_labels(), in setup.R), so the curves are for an
## average item. For en_exp that includes its by-item slopes.
get_prediction_grid_length <- function(
    fits,
    lang,
    model_type,
    n_words=c(4, 40),
    dl_distance=c(0, 1, 3)
) {
  expand.grid(
    n_words = seq(n_words[1] - .25, n_words[2] + .25, 1),
    dl_distance = dl_distance
  ) %>%
    mutate(
      weighted_dl_distance = dl_distance / n_words,
      grammatical_numeric = ifelse(dl_distance == 0, .5, -.5),
      item=factor(fits[[lang]][[model_type]]$model$item[1], levels=levels(fits[[lang]][[model_type]]$model$item))
    )
}

### ---- Baseline data for plotting (item means) ----
## Same aggregation as the model data, but deliberately WITHOUT the
## filter(dl_distance <= 3): the scatter shows the shuffled items as well, in the
## grey "4+" bucket, as context. The models (and hence the fitted curves drawn
## over these points) are still fit on distances up to 3 only.
plot_data_by_item_length <- crit %>%
  mutate(
    weighted_dl_distance = crit$dl_distance / crit$n_words,
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
  group_by(language, item, n_words, dl_distance, weighted_dl_distance,
           grammatical, condition) %>%
  summarize(
    response_num = mean(response_num)
  ) %>%
  mutate(
    grammatical_numeric = ifelse(grammatical == "grammatical", .5, -.5)
  ) %>%
  dplyr::group_by(language, item, n_words, dl_distance) %>%
  dplyr::summarise(response = mean(response_num), .groups = "drop") %>%
  dplyr::mutate(
    dl_distance_bucket = bucket_dl(dl_distance),
    language           = label_language(language)
  )

### ---- Baseline scatterplot (item means) ----
p_baseline <- ggplot(
  plot_data_by_item_length,
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

ggsave(p_gam, filename = "figures/curves_length_gam.png",
       width    = 6,
       height   = 16,
       units    = "in",
       dpi      = 300)

### ---- Model-based plots ----
#### ---- Grammatical model ----
#### ---- Predictions from grammatical model ----
grammatical_predictions <- setNames(
    lapply(names(models_by_language_length), function(lang) {
      predict(
        models_by_language_length[[lang]]$grammatical,
        newdata = get_prediction_grid_length(models_by_language_length, lang, "grammatical", n_words = range(crit[crit$language == lang,]$n_words)),
        type = "response",
        se.fit = TRUE,
        exclude = random_effect_labels(models_by_language_length[[lang]]$grammatical)
      )
    }
  ),
  names(models_by_language_length)
)

## Combine list to tibble
grammatical_plot_data <- bind_rows(
  lapply(names(grammatical_predictions), function(lang) {
    preds <- grammatical_predictions[[lang]]
    grid  <- get_prediction_grid_length(models_by_language_length, lang, "grammatical", n_words = range(crit[crit$language == lang,]$n_words))

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
  filename = "figures/curves_length_grammatical.png",
  width    = 6,
  height   = 16,
  units    = "in",
  dpi      = 300
)

#### ---- Distance model ----
#### ---- Predictions from distance model ----
distance_predictions <- setNames(
  lapply(names(models_by_language_length), function(lang) {
    predict(
      models_by_language_length[[lang]]$distance,
      newdata = get_prediction_grid_length(models_by_language_length, lang, "distance", n_words = range(crit[crit$language == lang,]$n_words)),
      type = "response",
      se.fit = TRUE,
      exclude = random_effect_labels(models_by_language_length[[lang]]$distance)
    )
  }
  ),
  names(models_by_language_length)
)

## Combine list to tibble
distance_plot_data <- bind_rows(
  lapply(names(distance_predictions), function(lang) {
    preds <- distance_predictions[[lang]]
    grid  <- get_prediction_grid_length(models_by_language_length, lang, "distance", n_words = range(crit[crit$language == lang,]$n_words))

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
  filename = "figures/curves_length_distance.png",
  width    = 6,
  height   = 16,
  units    = "in",
  dpi      = 300
)

#### ---- Weighted distance model ----
#### ---- Predictions from weighted distance model ----
weighted_distance_predictions <- setNames(
  lapply(names(models_by_language_length), function(lang) {
    predict(
      models_by_language_length[[lang]]$weighted_distance,
      newdata = get_prediction_grid_length(models_by_language_length, lang, "weighted_distance", n_words = range(crit[crit$language == lang,]$n_words)),
      type = "response",
      se.fit = TRUE,
      exclude = random_effect_labels(models_by_language_length[[lang]]$weighted_distance)
    )
  }
  ),
  names(models_by_language_length)
)

## Combine list to tibble
weighted_distance_plot_data <- bind_rows(
  lapply(names(weighted_distance_predictions), function(lang) {
    preds <- weighted_distance_predictions[[lang]]
    grid  <- get_prediction_grid_length(models_by_language_length, lang, "weighted_distance", n_words = range(crit[crit$language == lang,]$n_words))

    tibble(
      language = label_language(lang),
      n_words = grid$n_words,
      dl_distance_bucket = bucket_dl(grid$dl_distance),
      response = preds$fit,
      se = preds$se.fit
    )
  })
)

##### ---- Plot weighted distance (communication-based) model ----
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
  filename = "figures/curves_length_weighted_distance.png",
  width    = 6,
  height   = 16,
  units    = "in",
  dpi      = 300
)

### ---- Combine plots into the curve figure ----
curves_length_figure = (p_gam |  p_weighted_distance | p_distance |p_grammatical) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
        )

curves_length_figure

ggsave(
  filename = "figures/curves_length.png",
  plot     = curves_length_figure,
  width    = 25,
  height   = 20,
  units    = "in",
  dpi      = 300
)

## Save as pdf
ggsave(
  filename = "figures/curves_length.pdf",
  plot     = curves_length_figure,
  width    = 25,
  height   = 20,
  units    = "in"
)
