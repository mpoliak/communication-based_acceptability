## =============================================================================
## SI Section 3: Quantifying information using length in tokens instead of
## surprisal
##
## The main text measures a sentence's information as its summed surprisal. This
## section repeats the main analyses with information measured as sentence
## length in tokens instead:
##   - mixed-effects models of DL distance x length, per experiment
##   - the effect of length on well-formed sentences, tested across experiments
##   - the five-model comparison (Binary, Gradient, Communication, SLOR, MORCELA),
##     per experiment and pooled over the eight language experiments (controlled
##     English left out of the pool), with the Communication model's information
##     term on length
##   - the curve figure: ratings and model predictions over length
## Produces Figure S5 (figures/curves_length.png) and Figure S6
## (figures/likelihood_differences_length.png), and writes its statistics to
## outputs/s3_length_*.csv. Reads ../analysis/outputs/likelihoods_*.csv, so run
## ../analysis/main_analysis.R first.
## =============================================================================

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

## ---- Linear Mixed-Effects Models ----
## The main text's regression with length in place of surprisal: ratings on DL
## distance, length and their interaction, fit separately for each experiment
## with both predictors z-scored within the experiment. en_exp is fit on its own
## below, with a fuller random-effects structure, so it is left out of this loop.
## As in ../analysis/main_analysis.R, each experiment uses the optimizer that
## allFit() shows at the best log-likelihood: nlminb, except bobyqa for English.
## There nlminb stops 1.7 REML log-likelihood units short of the optimum that
## bobyqa and both nloptwrap optimizers reach (within 1.5e-6), with a negative
## eigenvalue of the Hessian. In Mandarin, nlminb reports "singular convergence
## (7)" at the same optimum that bobyqa finds.
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
  control = if (current_language == "en") {
    lmerControl(
      optimizer ="bobyqa",
      optCtrl = list(maxfun = 1e5)
    )
  } else {
    lmerControl(
      optimizer ="nlminbwrap",
      optCtrl = list(eval.max = 1e5, iter.max = 1e5)
    )
  }
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
  ## With this random-effects structure bobyqa stops ~28 log-likelihood units
  ## short of the optimum and fails the convergence check. The fit is still
  ## singular, and nlminb reports convergence code 1, "singular convergence
  ## (7)", as in Mandarin.
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
## Length is left in raw tokens here (the main text z-scores surprisal), so the
## slopes are in rating points per token. en_exp is fit on its own below.
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
## Binary, Gradient and Communication are the rule-based and communication-based
## models of the main text; SLOR and MORCELA are its two language-model
## baselines. Only the Communication model differs from the main text's: here
## it is a function of DL distance per token and of length, there of DL
## distance per bit and of surprisal.
##
## The item means are setup.R's crit_by_items (distances above 3 dropped, en_exp
## item ids already stripped of their length suffix, so the three length
## versions of each base item share one item random effect).

### ---- Analysis By Language ----
if (file.exists("models/s3_length_by_language_5models.RData")) {
  load("models/s3_length_by_language_5models.RData")
} else {

  ## Fit a model for each language. SLOR and MORCELA get a monotonically
  ## INCREASING smooth: a higher score means the sentence is more probable than
  ## its words alone would predict, and so should be judged more acceptable.
  models_by_language_length <- list()
  for (language_iter in levels(crit_by_items$language)) {
    print(
      paste0("Fitting models for ", language_iter)
    )

    ## Get specific data. MORCELA's two weights are learned separately for each
    ## experiment, since there is no reason a single beta and gamma should suit
    ## all of them.
    crit_by_items_language_specific <- crit_by_items %>%
      filter(language == language_iter) %>%
      droplevels()

    crit_by_items_language_specific$morcela <- morcela_score(
      crit_by_items_language_specific, morcela_weights(crit_by_items_language_specific)
    )

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
stopifnot(
  setequal(names(models_by_language_length), levels(crit_by_items$language)),
  all(vapply(models_by_language_length, function(M) {
    identical(names(M), c("grammatical", "distance", "weighted_distance", "slor", "morcela")) &&
      all(c("s(weighted_dl_distance)", "s(n_words)") %in%
            vapply(M$weighted_distance$smooth, function(s) s$label, character(1)))
  }, logical(1)))
)

## Evaluate Convergence
print(lapply(models_by_language_length, function(M) sapply(M, function(m) m$conv)))

### ---- Analysis Across All Languages ----
## Each model is also fit once to the eight language experiments at the same
## time, with a random intercept for language (8 levels) and a by-language random
## slope for each predictor, on top of the item random intercept.
## Controlled English is left out of the pool and fit only on its own, above.
## MORCELA's two weights are estimated once here, on the pooled data: this is one
## model of the eight language experiments, so it gets one beta and one gamma.
crit_by_items_all_experiments <- crit_by_items %>%
  filter(language != "en_exp") %>%
  droplevels()

crit_by_items_all_experiments$morcela <- morcela_score(
  crit_by_items_all_experiments, morcela_weights(crit_by_items_all_experiments)
)

if (file.exists("models/s3_length_all_experiments_5models.RData")) {
  load("models/s3_length_all_experiments_5models.RData")
} else {

  print("Fitting models for all languages together")

  mod_grammatical <- scam(
    response_num ~ grammatical_numeric +  # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, grammatical_numeric, bs="re"), # RE
    data=crit_by_items_all_experiments
  )

  mod_distance <- scam(
    response_num ~ s(dl_distance, bs="mpd") + # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, dl_distance, bs="re"), # RE
    data=crit_by_items_all_experiments
  )

  mod_weighted_distance <- scam(
    response_num ~ s(weighted_dl_distance, bs="mpd") + s(n_words, bs="mpd") + # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, weighted_dl_distance, bs="re") +
      s(language, n_words, bs="re"), # RE
    data=crit_by_items_all_experiments
  )

  mod_slor <- scam(
    response_num ~ s(slor, bs="mpi") + # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, slor, bs="re"), # RE
    data=crit_by_items_all_experiments
  )

  mod_morcela <- scam(
    response_num ~ s(morcela, bs="mpi") + # Fixed effects
      s(item, bs="re") +
      s(language, bs="re") + s(language, morcela, bs="re"), # RE
    data=crit_by_items_all_experiments
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

## Make sure the cache holds the pooled length-based fits of the eight language
## experiments, without controlled English
stopifnot(
  nrow(mod_weighted_distance$model) == nrow(crit_by_items_all_experiments),
  nlevels(mod_weighted_distance$model$language) == 8,
  !"en_exp" %in% levels(mod_weighted_distance$model$language),
  all(c("s(weighted_dl_distance)", "s(n_words)") %in%
        vapply(mod_weighted_distance$smooth, function(s) s$label, character(1)))
)

## Evaluate Convergence
print(mod_grammatical$conv)
print(mod_distance$conv)
print(mod_weighted_distance$conv)
print(mod_slor$conv)
print(mod_morcela$conv)

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
    ## In reverse, as the figure's y axis runs from the bottom up
    language = fct_rev(label_language(language)),
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

## The figure label "English\n(controlled)" goes on one line in the CSV
write_csv(
  stats_models_by_language %>% mutate(language = gsub("\n", " ", language)),
  "outputs/s3_length_likelihoods_by_language.csv"
)

pooled_likelihoods <- tibble(
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

write_csv(pooled_likelihoods, "outputs/s3_length_likelihoods_all_experiments.csv")

## Binary, Gradient, SLOR and MORCELA have no information term, so they are the
## main text's models, and only Communication changes. The change: these fits
## (length) minus the main text's (surprisal; ../analysis/outputs, written by
## ../analysis/main_analysis.R), for the eight language experiments together
## (both pools leave out controlled English) and then in each experiment
pooled_likelihoods %>%
  mutate(model = as.character(model)) %>%
  left_join(
    read_csv("../analysis/outputs/likelihoods_all_experiments.csv", show_col_types = FALSE) %>%
      select(model, main_text_log_likelihood = log_likelihood),
    by = "model"
  ) %>%
  mutate(length_minus_surprisal = log_likelihood - main_text_log_likelihood) %>%
  select(model, log_likelihood, main_text_log_likelihood, length_minus_surprisal) %>%
  print()

stats_models_by_language %>%
  mutate(model = as.character(model)) %>%
  left_join(
    read_csv("../analysis/outputs/likelihoods_by_language.csv", show_col_types = FALSE) %>%
      select(language_code = experiment, model, main_text_log_likelihood = log_likelihood),
    by = c("language_code", "model")
  ) %>%
  mutate(length_minus_surprisal = log_likelihood - main_text_log_likelihood) %>%
  select(language_code, model, length_minus_surprisal) %>%
  pivot_wider(names_from = model, values_from = length_minus_surprisal) %>%
  print(width = Inf)

### ---- Likelihood differences figure ----
## Drawn exactly like the main text's five-model figure, with the pooled "All
## Languages" row on top of the by-experiment rows.
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

## Stack the combined "All Languages" panel (top) over the by-language panel
## (bottom).
likelihood_differences_figure <- p_pooled_likelihoods / p_language_likelihoods +
  plot_layout(heights = c(1, 3))

ggsave("figures/likelihood_differences_length.png",
       plot = likelihood_differences_figure, width = 8, height = 8)

## ---- Curve Figure (Figure S5) ----
## The main text's Figure 1 with length in place of surprisal: item-mean ratings
## over length, with the Binary, Gradient and Communication models' predictions
## (from the by-experiment fits above) drawn as curves for DL distances 0, 1 and
## 3. Drawn like Figure 1 (../analysis/main_analysis.R), at its print size
## (theme_print_8pt, setup.R), with the models in the order of the other
## figures. Predictions are population-level: random_effect_labels() drops every
## item effect, including the by-item slopes of the controlled English models.

### ---- Baseline scatterplot (item means) ----
## plot_data_by_item (setup.R) deliberately keeps the items at DL >= 4: the
## scatter and the human-judgment smooths show them in the grey "4+" bucket, as
## context. The models (and hence their fitted curves) are still fit on
## distances up to 3 only. Length, like surprisal, is a property of the item.
p_baseline <- ggplot(
  plot_data_by_item,
  aes(x = n_words, y = response, color = dl_distance_bucket)
) +
  geom_point(alpha = 0.3, size = 0.6, stroke = 0) +
  facet_wrap(~language, ncol = 1, strip.position = "left") +
  labs(x = "Sentence length (tokens)", y = NULL) +
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
    aes(x = n_words, y = response, color = dl_distance_bucket),
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
## Each model's predictions over the experiment's range of length, at DL
## distances 0, 1 and 3, one experiment at a time
binary_plot_data <- bind_rows(lapply(names(models_by_language_length), function(lang) {
  fit <- models_by_language_length[[lang]]$grammatical
  grid <- get_prediction_grid(
    models_by_language_length, lang, "grammatical",
    n_words = range(crit$n_words[crit$language == lang])
  )
  preds <- predict(fit, newdata = grid, type = "response", se.fit = TRUE,
                   exclude = random_effect_labels(fit))
  tibble(
    language = label_language(lang),
    n_words = grid$n_words,
    dl_distance_bucket = bucket_dl(grid$dl_distance),
    response = preds$fit,
    se = preds$se.fit
  )
}))

gradient_plot_data <- bind_rows(lapply(names(models_by_language_length), function(lang) {
  fit <- models_by_language_length[[lang]]$distance
  grid <- get_prediction_grid(
    models_by_language_length, lang, "distance",
    n_words = range(crit$n_words[crit$language == lang])
  )
  preds <- predict(fit, newdata = grid, type = "response", se.fit = TRUE,
                   exclude = random_effect_labels(fit))
  tibble(
    language = label_language(lang),
    n_words = grid$n_words,
    dl_distance_bucket = bucket_dl(grid$dl_distance),
    response = preds$fit,
    se = preds$se.fit
  )
}))

communication_plot_data <- bind_rows(lapply(names(models_by_language_length), function(lang) {
  fit <- models_by_language_length[[lang]]$weighted_distance
  grid <- get_prediction_grid(
    models_by_language_length, lang, "weighted_distance",
    n_words = range(crit$n_words[crit$language == lang])
  )
  preds <- predict(fit, newdata = grid, type = "response", se.fit = TRUE,
                   exclude = random_effect_labels(fit))
  tibble(
    language = label_language(lang),
    n_words = grid$n_words,
    dl_distance_bucket = bucket_dl(grid$dl_distance),
    response = preds$fit,
    se = preds$se.fit
  )
}))

### ---- Plot binary grammaticality model ----
p_binary <- p_baseline +
  ggnewscale::new_scale_color() +   # <- reset color scale AFTER baseline
  geom_ribbon(
    data = binary_plot_data,
    aes(
      x = n_words,
      ymin = response - 2 * se,
      ymax = response + 2 * se,
      group = dl_distance_bucket
    ),
    alpha = 0.25
  ) +
  geom_line(
    data = binary_plot_data,
    aes(
      x = n_words,
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
      x = n_words,
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
      x = n_words,
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
      x = n_words,
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
      x = n_words,
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
  ggtitle(MODEL_LABELS[3]) +
  scale_y_continuous(limits = YLIM, breaks = NULL) +
  ggplot2::guides(color = "none")

### ---- Combine the panels into Figure S5 ----
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
  filename = "figures/curves_length.png",
  plot     = prediction_figure,
  width    = 7,
  height   = 7.5,
  units    = "in",
  dpi      = 300
)
