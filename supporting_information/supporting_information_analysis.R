## =============================================================================
## Supporting Information analyses for
##   "Evidence for a communication-based theory of linguistic acceptability judgments"
##
## Organized to mirror the SI document:
##   Section 1 - Controlling for Number of Variables (length as a predictor)
##   Section 2 - Analyses including the shuffled condition
##   Section 3 - Quantifying Information Using Surprisal instead of Length
##   Section 4 - Error awareness (English): how well participants detect the
##               errors we introduce, and whether the acceptability results
##               hold up once detection is taken into account
## followed by 'Diagnostic and Additional Code' (materials & participant
## descriptives, parametric linear models, the maximal-exclusion robustness
## check, and diagnostic comparison figures) that is not part of the SI document.
##
## Run in a UTF-8 locale (default in RStudio; for Rscript use e.g.
## LC_ALL=en_US.UTF-8). Fitted models are cached under models/; delete that
## folder's contents to refit from scratch.
## =============================================================================

library(tidyverse)
library(lme4)
library(lmerTest)
library(mgcv)
library(scam)
library(patchwork)

## Set plotting preferences
ggplot2::theme_set(ggplot2::theme_bw(base_size=18))
ggplot2::theme_update(panel.grid = ggplot2::element_blank(), 
                      strip.background = ggplot2::element_blank(),
                      legend.key = ggplot2::element_blank(),
                      panel.border = ggplot2::element_blank(),
                      axis.line = ggplot2::element_line(),
                      strip.text = ggplot2::element_text(face = "bold"),
                      plot.title = element_text(hjust = 0.5))

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

## Create the output folders if they are absent. models/ holds git-ignored model
## caches, so it does not exist in a fresh clone and save() below would fail.
dir.create("models", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)

## Helper functions
label_language <- function(language_code) {
  factor(dplyr::case_when(
    language_code == "da" ~ "Danish",
    language_code == "de" ~ "German",
    language_code == "en" ~ "English",
    language_code == "hi" ~ "Hindi",
    language_code == "fr" ~ "French",
    language_code == "ko" ~ "Korean",
    language_code == "ru" ~ "Russian",
    language_code == "zh" ~ "Mandarin",
    language_code == "en_exp" ~ "English\n(controlled)"
  ),
  ## Use the reverse of this order: levels=c("English\n(controlled)", "Russian", "Mandarin", "Korean", "Hindi", "German", "French", "English", "Danish")
  levels=c("Danish", "English", "French", "German", "Hindi", "Korean", 
           "Mandarin", "Russian", "English\n(controlled)")
  )
}


## ---- Marginal log-likelihood for scam/gam models with random effects ----
## scam/gam's logLik() is CONDITIONAL on the fitted random effects b; this returns
## the MARGINAL log-likelihood (integrating the b's out), i.e. the quantity that
## lmer() reports, so that model comparison does not credit a model for its fitted
## item- or language-specific random effects. Validated against lmer() for every
## random-effects structure used here (item/language intercepts and slopes) in
## ../marginal_loglik_test.R (outside this repo).
marginal_logLik <- function(m) {
  scale  <- m$sig2
  y      <- if (!is.null(m$y)) m$y else model.response(m$model)
  n      <- length(y)
  re_idx <- which(vapply(m$smooth, function(s) inherits(s, "random.effect"), logical(1)))
  if (!length(re_idx)) return(as.numeric(logLik(m)))       # no REs: nothing to marginalize
  re_lab <- vapply(m$smooth[re_idx], function(s) s$label, character(1))
  mu <- as.numeric(predict(m, exclude = re_lab))           # fixed part with random effects off
  Xp <- predict(m, type = "lpmatrix")                      # only the random-effect columns used below
  V  <- diag(scale, n)                                     # residual variance
  for (k in re_idx) {
    sm   <- m$smooth[[k]]
    Zk   <- Xp[, sm$first.para:sm$last.para, drop = FALSE]
    vark <- scale / (m$sp[[sm$label]] * diag(sm$S[[1]]))   # RE coefficient variance (accounts for
    V    <- V + Zk %*% (vark * t(Zk))                      #   mgcv's penalty rescaling of slopes)
  }
  r  <- y - mu
  ch <- chol(V)
  -0.5 * (n * log(2 * pi) + 2 * sum(log(diag(ch))) + sum(backsolve(ch, r, transpose = TRUE)^2))
}

## ---- By-language stats tibble with MARGINAL log-likelihoods ----
## Rebuilds a per-language stats tibble (for the likelihood-difference figures)
## from a fitted model list, using marginal_logLik() rather than the conditional
## logLik(). Called after each cache block so the figures use marginal values
## whether the models were just fit or loaded from disk.
marginal_stats_bylang <- function(model_list, model_map) {
  do.call(rbind, lapply(names(model_list), function(lg) {
    M <- model_list[[lg]]
    tibble(language = lg, model = names(M),
           log_likelihood = vapply(M, marginal_logLik, numeric(1)))
  })) %>%
    mutate(language = label_language(language),
           model = factor(model_map[model], levels = unname(model_map)))
}

## ---- Read and Prepare Data ----
## Read full data
dat <- read_csv("../data/full_data.csv")
surprisals <- read_csv("goldfish_surprisals/surprisal_goldfish.csv")

## Examine comprehension questions
comprehension <- dat %>% filter(grepl("comprehension", experiment))

dat <- dat %>% 
  left_join(
    surprisals %>% select(-sentence) %>% mutate(language = ifelse(language == "en_experimental", "en_exp", language)),
    by=c("language", "item")
  )

## Exclude participants with a low comprehension score
exclusions <- comprehension %>%
  group_by(participant_id) %>%
  mutate(
    correct_response = literal_response == response
    ) %>%
  summarize(
    accuracy = mean(correct_response)
  ) %>% 
  filter(
    accuracy < .8
  ) %>%
  pull(participant_id)

dat <- dat %>% 
  filter(!participant_id %in% exclusions)

## Get critical items
crit <- dat %>% 
  filter(experiment == "exp_b") %>%
  mutate(
    participant_id = as.factor(participant_id),
    item = as.factor(item),
    language = as.factor(language)
  )

## Remove nonsense conditions
crit <- crit %>% 
  filter(!grepl("non_sense", condition))

## Count number of participants per language
n_participants <- crit %>% 
  group_by(language) %>%
  summarize(
    n = n_distinct(participant_id)
  )

n_participants

table(crit$participant_id, crit$language)

## Verify the number of sentences per each language
n_sentences <- crit %>%
  group_by(language) %>%
  summarize(
    n = n_distinct(original_sentence)
  )

n_sentences

## Convert responses to numbers
crit$response_num <- as.numeric(substr(crit$response,1,1))

## Prepare the data frame
crit <- crit %>% 
  mutate(
    weighted_dl_distance = crit$dl_distance / crit$n_words,
    weighted_dl_distance_surprisal = crit$dl_distance / crit$sum_surprisal,
    grammatical = factor(
      ifelse(dl_distance == 0, "grammatical", "ungrammatical"),
      levels=c("ungrammatical", "grammatical")),
    item = as.factor(item),
    language = as.factor(language),
    response = as.factor(response_num)
  ) %>%
  select(
    participant_id, language, item, response, n_words, dl_distance, 
    weighted_dl_distance, grammatical, response_num, first_language,
    english_first_language, french_first_language, sum_surprisal, condition,
    weighted_dl_distance_surprisal
  )

## #############################################################################
## SI SECTION 1: Controlling for Number of Variables
## #############################################################################

## ---- Main Analysis with Length in all models ----
## Adds utterance length (P-spline) as a predictor to the binary- and gradient-
## grammaticality models, matching the communication model's predictor count.
## Produces withLength_figure.png and likelihood_differences_with_length.png.
crit_by_items_full <- crit %>% 
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
           grammatical, sum_surprisal, condition, weighted_dl_distance_surprisal) %>%
  summarize(
    response_num = mean(response_num)
  ) %>%
  mutate(
    grammatical_numeric = ifelse(grammatical == "grammatical", .5, -.5),
    item = as.factor(gsub("[sml]$", "", item))
  )

crit_by_items <- crit_by_items_full %>%
  filter(
    dl_distance <= 3
  )

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
#### ---- Helper Functions 
YLIM <- c(0, 8)

DL_COLORS <- c("0"="#EB7AA3", "1-2"="#1E88E5", "3"="#FFC107", "4+"="#BBBBBB")
DL_BUCKET_LEVELS <- c("0", "1-2", "3", "4+")

bucket_dl <- function(dl_distance) {
  dplyr::case_when(
    dl_distance == 0 ~ "0",
    dl_distance %in% 1:2 ~ "1-2",
    dl_distance == 3 ~ "3",
    dl_distance >= 4 ~ "4+"
  ) %>%
    factor(levels = DL_BUCKET_LEVELS)
}

scale_dl_color <- function() {
  ggplot2::scale_color_manual(
    name   = "Damerau-Levenshtein Distance",
    values = DL_COLORS,
    breaks = DL_BUCKET_LEVELS,
    limits = DL_BUCKET_LEVELS,
    drop   = FALSE,
    guide  = ggplot2::guide_legend(
      override.aes = list(alpha = 1, linewidth = 0, size = 3)
    )
  )
}


get_prediction_grid <- function(
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
      item=factor(fits_byLang_full_withLength[[lang]][[model_type]]$model$item[1], levels=levels(fits_byLang_full_withLength[[lang]][[model_type]]$model$item))
    )
}


### ---- Baseline data for plotting (item means) ----
plot_data_by_item <- crit %>% 
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
      newdata = get_prediction_grid(lang, "binary", n_words = range(crit[crit$language == lang,]$n_words)),
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
    grid  <- get_prediction_grid(lang, "binary", n_words = range(crit[crit$language == lang,]$n_words))
    
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
      newdata = get_prediction_grid(lang, "gradient", n_words = range(crit[crit$language == lang,]$n_words)),
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
    grid  <- get_prediction_grid(lang, "gradient", n_words = range(crit[crit$language == lang,]$n_words))
    
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
      newdata = get_prediction_grid(lang, "comm", n_words = range(crit[crit$language == lang,]$n_words)),
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
    grid  <- get_prediction_grid(lang, "comm", n_words = range(crit[crit$language == lang,]$n_words))
    
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

## #############################################################################
## SI SECTION 2: Analyses including the shuffled condition
## #############################################################################

## ---- Model Comparison with the Shuffled Condition ----
## Repeats the main model comparison on the full data, including the shuffled
## (DL distance >= 4) condition that the main text excludes.
## Produces with_shuffle_figure.png and likelihood_differences_with_shuffled.png.
if (!file.exists("models/main_model_fits_with_shuffled.RData")) {
  model_grammatical_incl_shuffled <- scam(
    response_num ~ grammatical_numeric + # + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, grammatical_numeric, bs="re"), # RE
    data=crit_by_items_full %>% 
      filter(language != "en_exp")
  )
  
  model_distance_incl_shuffled <- scam(
    response_num ~ s(dl_distance, bs="mpd") + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, dl_distance, bs="re"), # RE
    data=crit_by_items_full %>% 
      filter(language != "en_exp")
  )
  
  model_weighted_distance_incl_shuffled <- scam(
    response_num ~ s(weighted_dl_distance, bs="mpd") + s(n_words, bs="mpd") + # Fixed effects
      s(item, bs="re") + s(language, bs="re") + s(language, weighted_dl_distance, bs="re") + 
      s(language, n_words, bs="re"), # RE
    data=crit_by_items_full %>% 
      filter(language != "en_exp")
  )
  
  
  save(
    model_grammatical_incl_shuffled, 
    model_distance_incl_shuffled, 
    model_weighted_distance_incl_shuffled,
    file="models/main_model_fits_with_shuffled.RData"
  )
  
} else {
  load("models/main_model_fits_with_shuffled.RData")
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

p_incl_shuffled_main_likelihoods <- ggplot(incl_shuffled_main_likelihoods) +
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

### ---- Analysis By Language ----
if (file.exists("models/fits_byLang_full_with_shuffle.RData")) {
  load("models/fits_byLang_full_with_shuffle.RData")
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
        response_num ~ s(weighted_dl_distance, bs="mpd") + s(n_words, bs="mpd") +
          s(item, bs="re") + s(item, weighted_dl_distance, bs="re") + + s(item, n_words, bs="re"),
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
        response_num ~ s(weighted_dl_distance, bs="mpd") + s(n_words, bs="mpd") +
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
    file="models/fits_byLang_full_with_shuffle.RData"
  )
}

fits_byLang_full_with_shuffle_stats <- marginal_stats_bylang(
  fits_byLang_full_with_shuffle,
  c(binary = "Binary\nGrammaticality", gradient = "Gradient\nGrammaticality", comm = "Communication-based"))

fits_byLang_full_with_shuffle_stats <- fits_byLang_full_with_shuffle_stats %>%
  group_by(language) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood
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

## Stack the combined "All Languages" panel (top) over the by-language panel
## (bottom), matching the main-text likelihood figure.
likelihood_differences_with_shuffled_figure <- p_incl_shuffled_main_likelihoods / p_language_likelihoods +
  plot_layout(heights = c(1, 3))

ggsave("figures/likelihood_differences_with_shuffled.png",
       plot = likelihood_differences_with_shuffled_figure, width = 8, height = 8)


### ---- Baseline scatterplot (item means) ----
DL_COLORS <- c("0"="#EB7AA3", "1-2"="#1E88E5", "3"="#FFC107", "4+"="#004D40")

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


### ---- Model-based plots ----
#### ---- Grammatical model ----
#### ---- Predictions from grammatical model ----
grammatical_predictions <- setNames(
  lapply(names(fits_byLang_full_with_shuffle), function(lang) {
    predict(
      fits_byLang_full_with_shuffle[[lang]]$binary,
      newdata = get_prediction_grid(lang, "binary", n_words = range(crit[crit$language == lang,]$n_words),
                                    dl_distance=c(0, 1, 3, 10)),
      type = "response",
      se.fit = TRUE,
      exclude = "s(item)"
    )
  }
  ),
  names(fits_byLang_full_with_shuffle)
)

## Combine list to tibble
grammatical_plot_data <- bind_rows(
  lapply(names(grammatical_predictions), function(lang) {
    preds <- grammatical_predictions[[lang]]
    grid  <- get_prediction_grid(lang, "binary", n_words = range(crit[crit$language == lang,]$n_words),
                                 dl_distance=c(0, 1, 3, 10))
    
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
      newdata = get_prediction_grid(lang, "gradient", n_words = range(crit[crit$language == lang,]$n_words),
                                    dl_distance=c(0, 1, 3, 10)),
      type = "response",
      se.fit = TRUE,
      exclude = "s(item)"
    )
  }
  ),
  names(fits_byLang_full_with_shuffle)
)

## Combine list to tibble
distance_plot_data <- bind_rows(
  lapply(names(distance_predictions), function(lang) {
    preds <- distance_predictions[[lang]]
    grid  <- get_prediction_grid(lang, "gradient", n_words = range(crit[crit$language == lang,]$n_words),
                                 dl_distance=c(0, 1, 3, 10))
    
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
      newdata = get_prediction_grid(lang, "comm", n_words = range(crit[crit$language == lang,]$n_words),
                                    dl_distance=c(0, 1, 3, 10)),
      type = "response",
      se.fit = TRUE,
      exclude = "s(item)"
    )
  }
  ),
  names(fits_byLang_full_with_shuffle)
)

## Combine list to tibble
weighted_distance_plot_data <- bind_rows(
  lapply(names(weighted_distance_predictions), function(lang) {
    preds <- weighted_distance_predictions[[lang]]
    grid  <- get_prediction_grid(lang, "comm", n_words = range(crit[crit$language == lang,]$n_words),
                                 dl_distance=c(0, 1, 3, 10))
    
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


## #############################################################################
## SI SECTION 3: Quantifying Information Using Surprisal instead of Length
## #############################################################################

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
    method="lm"
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
      item=factor(fits_byLang_full_with_surprisal[[lang]][[model_type]]$model$item[1], levels=levels(fits_byLang_full_with_surprisal[[lang]][[model_type]]$model$item))
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
      newdata = get_prediction_grid(lang, "binary", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal)),
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
    grid  <- get_prediction_grid(lang, "binary", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal))
    
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
      newdata = get_prediction_grid(lang, "gradient", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal)),
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
    grid  <- get_prediction_grid(lang, "gradient", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal))
    
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
      newdata = get_prediction_grid(lang, "comm", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal)),
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
    grid  <- get_prediction_grid(lang, "comm", sum_surprisal = range(crit[crit$language == lang,]$sum_surprisal))
    
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


## #############################################################################
## SI SECTION 4: Error awareness (English)
## #############################################################################
## A separate English experiment with 80 participants, on the same materials as
## the original English experiment. On every trial participants both rated the
## sentence's acceptability (exp_b) and answered "Are there any errors in this
## sentence?" (exp_b_comprehension). `literal_response` holds the correct Yes/No
## answer -- "No" for the original condition, "Yes" once any adjacent-token
## exchanges were introduced -- so a trial counts as correctly judged exactly
## when response == literal_response.
##
## The section asks, in order:
##   4.1  How well do participants detect the errors we introduced?
##   4.2  Are acceptability ratings independent of whether the error was noticed?
##   4.3  Do ratings on correctly-judged trials look like the original experiment?
##   4.4  Do the paper's three models rank the same way on those trials?
##   4.5  Does that ranking survive adding back wrongly-flagged original sentences?
## 4.4/4.5 also fit the same three models to the unfiltered trials and to the
## original English experiment, because restricting to correctly-judged trials
## selects on a variable correlated with the rating and the rankings cannot be
## read without those baselines.
##
## Every model is restricted to dl_distance <= 3, matching the main analyses (in
## practice this keeps original/swap_1/swap_3 and drops all but a handful of
## shuffled trials). Produces error_awareness_detection_by_length.png,
## error_awareness_acceptability_by_detection.png,
## error_awareness_acceptability_by_version.png and
## error_awareness_likelihood_differences.png.

## ---- Shared data preparation ----
## DL distance and sentence length are not recorded in the error-awareness data
## files: the experiment reuses the original English materials verbatim, so both
## are joined in from `crit` (which already applies the >= 80% comprehension
## exclusion) rather than recomputed. item ids are "en_17" there and "17" here.
ea_key <- crit %>%
  filter(language == "en") %>%
  distinct(item, condition, dl_distance, n_words) %>%
  mutate(item = as.integer(str_extract(as.character(item), "[0-9]+"))) %>%
  rename(cond = condition)

ea_raw <- read_csv("../data/en_error_awareness/data_en_error_awareness.csv")

## One row per participant x item: the acceptability rating joined to that same
## trial's Yes/No error judgment. This single table drives the whole section.
ea_trials <- ea_raw %>%
  filter(experiment == "exp_b") %>%
  transmute(participant_id, item = as.integer(item), cond,
            acceptability = as.numeric(substr(response, 1, 1))) %>%
  inner_join(
    ea_raw %>%
      filter(experiment == "exp_b_comprehension") %>%
      transmute(participant_id, item = as.integer(item), cond,
                said_error = response == "Yes",
                detected = response == literal_response),
    by = c("participant_id", "item", "cond")) %>%
  left_join(ea_key, by = c("item", "cond")) %>%
  mutate(
    condition_label = factor(cond,
      levels = c("original", "swap_1", "swap_3", "random"),
      labels = c("Original", "1 swap", "3 swaps", "Shuffled")),
    detection_label = factor(
      ifelse(detected, "Error judged correctly", "Error judged incorrectly"),
      levels = c("Error judged correctly", "Error judged incorrectly")),
    weighted_dl_distance = dl_distance / n_words,
    grammatical = factor(ifelse(dl_distance == 0, "grammatical", "ungrammatical"),
                         levels = c("ungrammatical", "grammatical")))

## Colours matching the main-text DL-distance buckets (0 / 1-2 / 3 / 4+), mapped
## to the four error conditions so this section matches the rest of the study.
error_colors <- c("Original" = "#EB7AA3", "1 swap" = "#1E88E5",
                  "3 swaps" = "#FFC107", "Shuffled" = "#BBBBBB")

## Every figure in this section plots item means, so collapse them in one place.
item_means <- function(d, value, ...) {
  d %>%
    group_by(...) %>%
    summarize(value = mean({{ value }}), n_trials = n(), .groups = "drop")
}

## Predictors are z-scored inside each fitting subset: it helps the maximal
## random-effect structures converge, and a linear rescaling of a predictor
## leaves the log-likelihood unchanged (the random-effect covariances are
## unstructured), so the Section 4.4/4.5 comparisons are unaffected.
zscore <- function(x) as.numeric(scale(x))

lmer_ctrl  <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
glmer_ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))

## Genuine convergence failure, as opposed to a boundary (singular) fit. A
## singular fit is still the maximum-likelihood solution -- it just puts a
## variance component or correlation on the boundary -- so it is not a reason to
## simplify the random-effect structure and is ignored here.
converged <- function(fit) {
  msgs <- fit@optinfo$conv$lme4$messages
  length(msgs[!grepl("singular", msgs, ignore.case = TRUE)]) == 0
}

## Keep the maximal random-effect structure the design justifies, and step down
## `formulas` only if the optimizer actually fails to converge on it. Reports
## which structure was used so any simplification is visible rather than
## buried in a comment.
fit_maximal <- function(formulas, fitter, label) {
  for (i in seq_along(formulas)) {
    fit <- suppressWarnings(fitter(formulas[[i]]))
    ok <- converged(fit)
    if (ok || i == length(formulas)) {
      note <- if (!ok) " -- still not converged, reported as is"
              else if (i == 1) " (maximal)" else ""
      cat(sprintf("[%s] random-effect structure %d of %d%s\n", label, i,
                  length(formulas), note))
      return(fit)
    }
  }
}

## The models restrict to dl_distance <= 3, as the main analyses do.
ea_model_trials <- ea_trials %>%
  filter(dl_distance <= 3) %>%
  mutate(participant_id = factor(participant_id), item = factor(item)) %>%
  droplevels()


## ---- 4.1 Error detection: descriptives, figure, and model ----
cat("\n===== 4.1 Error detection =====\n")

## Participants, items and trials
ea_trials %>%
  summarize(participants = n_distinct(participant_id), items = n_distinct(item),
            trials = n()) %>%
  print()

ea_trials %>% count(condition_label) %>% print()

## Detection accuracy by condition. For the original condition "accuracy" is one
## minus the false-alarm rate (correctly saying there is no error); for the three
## error conditions it is the hit rate.
ea_trials %>%
  group_by(condition_label) %>%
  summarize(n = n(), accuracy = mean(detected), said_error = mean(said_error),
            mean_rating = mean(acceptability), .groups = "drop") %>%
  print(n = Inf)

ea_trials %>%
  summarize(overall_accuracy = mean(detected),
            false_alarm_rate = mean(said_error[cond == "original"]),
            hit_rate = mean(said_error[cond != "original"])) %>%
  print()

## Spread of by-item detection accuracy (the dots in the figure below)
ea_trials %>%
  item_means(detected, condition_label, item) %>%
  group_by(condition_label) %>%
  summarize(min = min(value), median = median(value), max = max(value),
            items_at_ceiling = sum(value == 1), .groups = "drop") %>%
  print(n = Inf)

## Figure: by-item detection accuracy against sentence length, one curve per
## condition. Dots are item means; the curve is a binomial GAM over the
## trial-level 0/1 data. The basis is deliberately small and shrinkage-penalised
## (k = 4, bs = "cs") so the curves describe the trend rather than tracking
## individual items -- the shuffled condition sits at 99% accuracy and an
## unconstrained smooth there is all wiggle and no signal.
p_detection <- ggplot(mapping = aes(x = n_words, color = condition_label,
                                    fill = condition_label)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey70") +
  geom_point(data = item_means(ea_trials, detected, item, n_words,
                               condition_label),
             aes(y = value), alpha = 0.35) +
  geom_smooth(data = ea_trials %>% mutate(correct = as.integer(detected)),
              aes(y = correct), method = "gam",
              method.args = list(family = "binomial"),
              formula = y ~ s(x, bs = "cs", k = 4)) +
  scale_color_manual(name = "Condition", values = error_colors,
                     aesthetics = c("colour", "fill")) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Error-detection accuracy", labels = scales::percent,
                     breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1))
p_detection
ggsave("figures/error_awareness_detection_by_length.png", p_detection,
       width = 8, height = 5.5)

## Model: does the chance of noticing an error depend on how much was displaced
## and on how long the sentence is? Same fixed and random structure as the
## main-analysis acceptability model, in binomial form.
detection_model <- fit_maximal(
  list(
    detected ~ dl_distance_z * n_words_z +
      (dl_distance_z * n_words_z | participant_id) + (dl_distance_z | item),
    detected ~ dl_distance_z * n_words_z +
      (dl_distance_z + n_words_z | participant_id) + (dl_distance_z | item),
    detected ~ dl_distance_z * n_words_z +
      (dl_distance_z + n_words_z || participant_id) + (dl_distance_z || item),
    detected ~ dl_distance_z * n_words_z +
      (dl_distance_z + n_words_z || participant_id) + (1 | item)),
  function(f) glmer(f, family = binomial, control = glmer_ctrl,
                    data = ea_model_trials %>%
                      mutate(dl_distance_z = zscore(dl_distance),
                             n_words_z = zscore(n_words))),
  "detection")

cat("\n---- detection ~ DL distance x length ----\n")
print(round(summary(detection_model)$coefficients, 4))


## ---- 4.2 Are acceptability ratings independent of error awareness? ----
cat("\n===== 4.2 Ratings and error awareness =====\n")

## Figure: item-mean acceptability against sentence length, coloured by
## condition, split into the trials whose error judgment was right and the
## trials whose judgment was wrong.
##
## The shuffled condition is dropped here (only here -- it is still shown in the
## 4.1 and 4.3 figures). Shuffled sentences are detected on 99% of trials, so the
## "judged incorrectly" panel holds only a handful of them spread over a few
## items, and the fit through them is noise.
p_awareness <- ea_trials %>%
  filter(cond != "random") %>%
  mutate(condition_label = droplevels(condition_label)) %>%
  item_means(acceptability, item, n_words, condition_label, detection_label) %>%
  ggplot(aes(n_words, value, color = condition_label, fill = condition_label)) +
  geom_point(alpha = 0.35) +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(~ detection_label, nrow = 1) +
  scale_color_manual(name = "Condition", values = error_colors,
                     aesthetics = c("colour", "fill")) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Acceptability", breaks = seq(1, 7, 2)) +
  coord_cartesian(ylim = c(1, 7))
p_awareness
ggsave("figures/error_awareness_acceptability_by_detection.png", p_awareness,
       width = 11, height = 5.5)

## Model: acceptability from the communication predictors (DL distance per token,
## plus sentence length) crossed with whether the error judgment was correct.
## Correctness is sum-coded so the predictor main effects are averages over the
## two judgment outcomes. Sentence length is between items here, so it gets a
## by-participant but no by-item slope (as in main_analysis.R).
awareness_model <- fit_maximal(
  list(
    acceptability ~ (weighted_dl_distance_z + n_words_z) * detected_c +
      (weighted_dl_distance_z + n_words_z + detected_c | participant_id) +
      (weighted_dl_distance_z + detected_c | item),
    acceptability ~ (weighted_dl_distance_z + n_words_z) * detected_c +
      (weighted_dl_distance_z + n_words_z + detected_c | participant_id) +
      (weighted_dl_distance_z | item),
    acceptability ~ (weighted_dl_distance_z + n_words_z) * detected_c +
      (weighted_dl_distance_z + n_words_z + detected_c || participant_id) +
      (weighted_dl_distance_z || item),
    acceptability ~ (weighted_dl_distance_z + n_words_z) * detected_c +
      (weighted_dl_distance_z + n_words_z + detected_c || participant_id) +
      (1 | item)),
  function(f) lmer(f, control = lmer_ctrl,
                   data = ea_model_trials %>%
                     mutate(weighted_dl_distance_z = zscore(weighted_dl_distance),
                            n_words_z = zscore(n_words),
                            detected_c = ifelse(detected, 0.5, -0.5))),
  "awareness")

cat("\n---- acceptability ~ (DL/length + length) x judged correctly ----\n")
print(round(summary(awareness_model)$coefficients, 4))


## ---- 4.3 Ratings on correctly-judged trials vs. the original experiment ----
cat("\n===== 4.3 Correctly-judged trials vs. the original experiment =====\n")

## Figure: the same by-item picture drawn twice -- for the original English
## experiment (which asked meaning-comprehension questions, so every trial is
## included) and for the error-awareness experiment restricted to the trials
## whose error judgment was correct.
version_panels <- bind_rows(
  crit %>%
    filter(language == "en") %>%
    transmute(panel = "English\n(comprehension questions)",
              item = as.integer(str_extract(as.character(item), "[0-9]+")),
              cond = condition, n_words, acceptability = response_num),
  ea_trials %>%
    filter(detected) %>%
    transmute(panel = "Error awareness\n(errors correctly judged)",
              item, cond, n_words, acceptability)) %>%
  mutate(
    panel = factor(panel, levels = c("English\n(comprehension questions)",
                                     "Error awareness\n(errors correctly judged)")),
    condition_label = factor(cond,
      levels = c("original", "swap_1", "swap_3", "random"),
      labels = c("Original", "1 swap", "3 swaps", "Shuffled")))

p_version <- ggplot(item_means(version_panels, acceptability, panel, item,
                               n_words, condition_label),
                    aes(n_words, value, color = condition_label,
                        fill = condition_label)) +
  geom_point(alpha = 0.35) +
  geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(~ panel, nrow = 1) +
  scale_color_manual(name = "Condition", values = error_colors,
                     aesthetics = c("colour", "fill")) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Acceptability", breaks = seq(1, 7, 2)) +
  coord_cartesian(ylim = c(1, 7))
p_version
ggsave("figures/error_awareness_acceptability_by_version.png", p_version,
       width = 11, height = 5.5)

## The communication account's second prediction, isolated: among sentences that
## the participant correctly judged to be error-free, longer should still be less
## acceptable. Fit on the original-condition trials only.
cat("\n---- error-free (correctly judged) acceptability ~ length ----\n")
print(round(summary(lmer(
  acceptability ~ n_words_c + (n_words_c | participant_id) + (1 | item),
  data = ea_model_trials %>%
    filter(cond == "original", detected) %>% droplevels() %>%
    mutate(n_words_c = n_words - mean(n_words)),
  control = lmer_ctrl))$coefficients, 4))


## ---- 4.4/4.5 The paper's three models on error-awareness trials ----
cat("\n===== 4.4/4.5 Model comparison =====\n")

## The three competing accounts, in both forms the paper uses them: the shape-
## constrained GAMs over by-item means (as in the main text and Sections 1-3) and
## the parametric linear mixed models over trials (as in the linear-model checks
## further below). Single-language versions, so the language random effects drop
## out. The lmers are fit with REML = FALSE because the three models differ in
## their fixed effects and REML deviances are not comparable across them.
fit_ea_model_set <- function(d) {
  d_items <- d %>%
    droplevels() %>%
    group_by(item, cond, n_words, dl_distance, weighted_dl_distance, grammatical) %>%
    summarize(response_num = mean(acceptability), .groups = "drop") %>%
    mutate(item = factor(item),
           grammatical_numeric = ifelse(grammatical == "grammatical", .5, -.5))
  d_trials <- d %>%
    droplevels() %>%
    mutate(dl_distance_z = zscore(dl_distance),
           weighted_dl_distance_z = zscore(weighted_dl_distance),
           n_words_z = zscore(n_words))
  list(
    scam = list(
      binary = scam(response_num ~ grammatical_numeric + s(item, bs = "re"),
                    data = d_items),
      gradient = scam(response_num ~ s(dl_distance, bs = "mpd") + s(item, bs = "re"),
                      data = d_items),
      comm = scam(response_num ~ s(weighted_dl_distance, bs = "mpd") +
                    s(n_words, bs = "mpd") + s(item, bs = "re"),
                  data = d_items)),
    lmer = list(
      binary = lmer(acceptability ~ grammatical +
                      (grammatical | participant_id) + (grammatical | item),
                    data = d_trials, REML = FALSE, control = lmer_ctrl),
      gradient = lmer(acceptability ~ dl_distance_z +
                        (dl_distance_z | participant_id) + (dl_distance_z | item),
                      data = d_trials, REML = FALSE, control = lmer_ctrl),
      comm = lmer(acceptability ~ weighted_dl_distance_z + n_words_z +
                    (weighted_dl_distance_z + n_words_z | participant_id) +
                    (weighted_dl_distance_z | item),
                  data = d_trials, REML = FALSE, control = lmer_ctrl)))
}

## Three nested views of the same trials:
##   All trials          -- everything, regardless of the error judgment
##   Detected only       -- 4.4: only trials whose error judgment was correct
##   + flagged originals -- 4.5: adds back the original sentences that were
##                          wrongly judged to contain an error
## plus the original English experiment (meaning-comprehension questions, no
## error judgment at all) as the reference the section is arguing against. The
## "All trials" and reference rows are not part of the questions as posed, but
## without them the 4.4/4.5 rankings cannot be read: restricting to correctly
## judged trials is a selection on a variable that is itself correlated with the
## rating, and the baselines are what show how much of the result that produces.
ea_reference <- crit %>%
  filter(language == "en", dl_distance <= 3) %>%
  transmute(participant_id = factor(participant_id),
            item = factor(as.integer(str_extract(as.character(item), "[0-9]+"))),
            cond = condition, n_words, dl_distance, weighted_dl_distance,
            grammatical, acceptability = response_num)

ea_datasets <- list(
  "English\n(comprehension questions)" = ea_reference,
  "All trials"          = ea_model_trials,
  "Detected only"       = ea_model_trials %>% filter(detected),
  "+ flagged originals" = ea_model_trials %>% filter(detected | cond == "original"))

if (!file.exists("models/ea_section4_models.RData")) {
  ea_model_sets <- lapply(names(ea_datasets), function(nm) {
    cat("Fitting model set for", gsub("\n", " ", nm), "\n")
    fit_ea_model_set(ea_datasets[[nm]])
  })
  names(ea_model_sets) <- names(ea_datasets)
  save(ea_model_sets, file = "models/ea_section4_models.RData")
} else {
  load("models/ea_section4_models.RData")
}

## What each dataset retains
tibble(dataset = gsub("\n", " ", names(ea_datasets)),
       trials = vapply(ea_datasets, nrow, integer(1)),
       items = vapply(ea_datasets, function(d) n_distinct(d$item), integer(1)),
       originals = vapply(ea_datasets, function(d) sum(d$cond == "original"),
                          integer(1))) %>%
  print(n = Inf)

## Convergence: scam reports its own flag; for the lmers we report genuine
## optimizer failures rather than silently trusting the likelihoods. Boundary
## (singular) fits are not reported -- see converged() above.
for (nm in names(ea_model_sets)) {
  cat("\n---- convergence --", gsub("\n", " ", nm), "----\n")
  cat("  scam converged:",
      paste(names(ea_model_sets[[nm]]$scam),
            vapply(ea_model_sets[[nm]]$scam, function(m) m$conv, logical(1)),
            collapse = "  "), "\n")
  for (mn in names(ea_model_sets[[nm]]$lmer)) {
    cat("  lmer", mn, "- converged:",
        converged(ea_model_sets[[nm]]$lmer[[mn]]), "\n")
  }
}

## Log-likelihoods. scam's own logLik() is conditional on the fitted random
## effects, so marginal_logLik() is used there; lmer's logLik() is already
## marginal. The two families are fit to different data (item means vs. trials)
## and each dataset has a different number of observations, so likelihoods are
## only ever compared WITHIN a family and WITHIN a row.
ea_likelihoods <- bind_rows(lapply(names(ea_model_sets), function(nm) {
  bind_rows(
    tibble(dataset = nm, family = "Shape-constrained GAM\n(item means)",
           model = names(ea_model_sets[[nm]]$scam),
           log_likelihood = vapply(ea_model_sets[[nm]]$scam, marginal_logLik,
                                   numeric(1))),
    tibble(dataset = nm, family = "Linear mixed model\n(trials)",
           model = names(ea_model_sets[[nm]]$lmer),
           log_likelihood = vapply(ea_model_sets[[nm]]$lmer,
                                   function(m) as.numeric(logLik(m)), numeric(1))))
})) %>%
  mutate(
    dataset = factor(dataset, levels = names(ea_datasets)),
    model = factor(model, levels = c("binary", "gradient", "comm"),
                   labels = c("Binary\nGrammaticality", "Gradient\nGrammaticality",
                              "Communication-based")),
    family = factor(family, levels = c("Linear mixed model\n(trials)",
                                       "Shape-constrained GAM\n(item means)"))) %>%
  group_by(dataset, family) %>%
  mutate(likelihood_diff = log_likelihood - max(log_likelihood)) %>%
  ungroup()

ea_likelihoods %>%
  mutate(dataset = gsub("\n", " ", dataset)) %>%
  select(family, dataset, model, log_likelihood, likelihood_diff) %>%
  arrange(family, dataset, model) %>%
  print(n = Inf)

## Which model wins in each cell
ea_likelihoods %>%
  group_by(family, dataset) %>%
  slice_max(log_likelihood, n = 1) %>%
  ungroup() %>%
  mutate(dataset = gsub("\n", " ", dataset)) %>%
  select(family, dataset, best_model = model, log_likelihood) %>%
  print(n = Inf)

## Figure: log-likelihood difference from the best model within each row, one row
## per dataset, one panel per model family. Same idiom and colours as the
## likelihood_differences_* figures elsewhere in the SI.
p_ea_likelihoods <- ggplot(ea_likelihoods) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
  geom_point(aes(y = dataset, x = likelihood_diff, color = model, shape = model),
             size = 3) +
  facet_wrap(~ family, scales = "free_x") +
  scale_color_manual(name = "", values = c("#EB7AA3", "#1E88E5", "#004D40")) +
  scale_shape_manual(name = "", values = c(16, 17, 15)) +
  scale_x_continuous(name = "Log-Likelihood Difference") +
  scale_y_discrete(name = "", limits = rev) +
  theme(legend.position = "top")
p_ea_likelihoods
ggsave("figures/error_awareness_likelihood_differences.png", p_ea_likelihoods,
       width = 12, height = 6)


## #############################################################################
## DIAGNOSTIC AND ADDITIONAL CODE
## Not part of the SI document. Materials & participant descriptives, parametric
## (full-random-effects) linear models, the maximal-exclusion robustness check,
## and diagnostic comparison figures.
## #############################################################################

## ---- Materials Information ----
### ---- Sentence Length Across Languages ----
lengths <- crit %>%
  select(language, item, n_words) %>%
  distinct() %>%
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
      levels=c("Danish", "English", "French", "German", 
               "Hindi", "Korean", "Mandarin", "Russian", 
               "English\n(controlled)")
      )
  )


lengths %>%
  group_by(language) %>% 
  summarize(
    median_length = median(n_words),
    min_length = min(n_words),
    max_length = max(n_words),
    sd_length = sd(n_words),
  )

ggplot(lengths) +
  geom_histogram(
    aes(x=n_words)
  ) + 
  facet_wrap(~language, scales="free_y") + 
  scale_y_continuous(breaks = scales::breaks_width(5))

ggsave("figures/sentence_length_histograms.png", width=10, height=6)

### ---- Comprehension Accuracy by Edit Distance ----
## Descriptive check: how accurately did participants answer the comprehension
## question that followed each sentence, as a function of how corrupted that
## sentence was? Restricted to the three error conditions in the design (0, 1, and
## 3 adjacent-token exchanges); the rare intermediate DL distance of 2 is dropped
## as it is only sparsely populated. Scored the same way as the exclusion
## criterion above (correct = literal_response == response), on the same retained
## sample (participants below the 80% threshold already excluded). Not part of the
## SI document.
comprehension_accuracy <- comprehension %>%
  filter(!participant_id %in% exclusions,
         dl_distance %in% c(0, 1, 3)) %>%
  mutate(correct = literal_response == response) %>%
  group_by(language, dl_distance) %>%
  summarize(n = n(), k = sum(correct, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    accuracy = k / n,
    ## Wilson 95% confidence interval for a binomial proportion (well behaved
    ## near the ceiling, where a normal-approximation interval would exceed 1)
    z = 1.96,
    wilson_center = (accuracy + z^2 / (2 * n)) / (1 + z^2 / n),
    wilson_half   = z * sqrt(accuracy * (1 - accuracy) / n + z^2 / (4 * n^2)) /
                    (1 + z^2 / n),
    lower = wilson_center - wilson_half,
    upper = wilson_center + wilson_half,
    language = label_language(language)
  )

## Pooled accuracy per DL distance, printed for reference
comprehension %>%
  filter(!participant_id %in% exclusions, dl_distance %in% c(0, 1, 3)) %>%
  mutate(correct = literal_response == response) %>%
  group_by(dl_distance) %>%
  summarize(accuracy = mean(correct, na.rm = TRUE), n = n(), .groups = "drop")

ggplot(comprehension_accuracy, aes(x = dl_distance, y = accuracy)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey70") +
  geom_line(color = "grey60") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.15) +
  geom_point(size = 2) +
  facet_wrap(~language) +
  scale_x_continuous(breaks = c(0, 1, 3), name = "Damerau-Levenshtein distance") +
  scale_y_continuous(name = "Comprehension accuracy", limits = c(0.5, 1)) +
  theme(panel.grid.major.y = element_line(color = "grey90"))

ggsave("figures/comprehension_accuracy.png", width = 10, height = 6)

### ---- Count Damerau-Levenshtein distance for each language ----
distance_counts <- crit %>%
  select(language, item, condition, dl_distance, n_words) %>%
  distinct() %>% mutate(
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
      levels=c("Danish", "English", "French", "German", 
               "Hindi", "Korean", "Mandarin", "Russian", 
               "English\n(controlled)")
    )
  )

table(distance_counts$language)
table(distance_counts$item, distance_counts$language)

sink("outputs/dl_distance_counts.txt")
print("Original Sentences")
table(
  distance_counts$dl_distance[distance_counts$condition == "original"], 
  distance_counts$language[distance_counts$condition == "original"]
) %>%
  print(., n=10^4)

print("1 Exchange")
table(
  distance_counts$dl_distance[distance_counts$condition == "swap_1"], 
  distance_counts$language[distance_counts$condition == "swap_1"]
) %>%
  print(., n=10^4)

print("3 Exchanges")
table(
  distance_counts$dl_distance[distance_counts$condition == "swap_3"], 
  distance_counts$language[distance_counts$condition == "swap_3"]
) %>%
  print(., n=10^4)

print("Shuffle")
table(
  distance_counts$dl_distance[distance_counts$condition == "random"], 
  distance_counts$language[distance_counts$condition == "random"]
) %>%
  print(., n=10^4)

sink()

## ---- Demographic Information ----
demographics <- dat %>%
  select(
    participant_id, country_birth, country_birth_text, first_language, 
    first_language_text,second_language , second_language_text, 
    language_contexts, english_first_language, country_of_birth, 
    french_first_language,  french_contexts, language
  ) %>%
  distinct

## Danish
danish_demographics <- demographics %>% filter(language == "da")
table(danish_demographics$first_language, danish_demographics$country_birth)

danish_demographics %>%
  filter(country_birth == "Andet")

## English
english_demographics <- demographics %>% filter(language == "en")

table(english_demographics$english_first_language, 
      english_demographics$country_of_birth)

## English (Length-Controlled Experiment)
english_exp_demographics <- demographics %>% filter(language == "en_exp")

table(english_exp_demographics$english_first_language, 
      english_exp_demographics$country_of_birth,
      useNA = "always"
)

## French
french_demographics <- demographics %>% filter(language == "fr")

table(
  french_demographics$french_first_language,
  french_demographics$country_of_birth
)

french_demographics %>%
  filter(country_of_birth == "Ailleurs") %>%
  pull(french_contexts)

## German
german_demographics <- demographics %>% filter(language == "de")

table(
  german_demographics$first_language,
  german_demographics$country_birth
)

german_demographics %>% 
  filter(country_birth == "Anderes") %>%
  pull(country_birth_text)

## Hindi
hindi_demographics <- demographics %>% filter(language == "hi")

table(
  hindi_demographics$first_language,
  hindi_demographics$country_birth
)

table(hindi_demographics$first_language_text)
table(hindi_demographics$country_birth_text, hindi_demographics$first_language)

## Korean
korean_demographics <- demographics %>% filter(language == "ko")

table(
  korean_demographics$first_language,
  korean_demographics$country_of_birth
)

## Chinese
chinese_demographics <- demographics %>% filter(language == "zh")

table(chinese_demographics$english_first_language)
table(chinese_demographics$country_of_birth)

table(
  chinese_demographics$english_first_language,
  chinese_demographics$country_of_birth
)

## Russian
russian_demographics <- demographics %>% filter(language == "ru")

table(
  russian_demographics$first_language,
  russian_demographics$country_birth
)

table(russian_demographics$country_birth,
      russian_demographics$country_birth_text)



## ---- Linear Models ----
### ---- All Materials ----
### ---- Collapsing Across Language ----
if (file.exists("models/linear_model_fits.Rdata")) {
  load("models/linear_model_fits.Rdata")
} else {
  
  model_grammatical_linear <- lmer(
    response_num ~ grammatical + 
      (grammatical|participant_id) + (grammatical|item) + (grammatical|language),
    data=crit %>%
      filter(
        language != "en_exp"
      ),
    control=lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e6)
    )
  )
  
  model_distance_linear <- lmer(
    response_num ~ dl_distance + 
      (dl_distance|participant_id) + (dl_distance|item) + (dl_distance|language),
    data=crit %>%
      filter(
        language != "en_exp"
      ),
    control=lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e6)
    )
  )
  
  model_weighted_distance_linear <- lmer(
    response_num ~ weighted_dl_distance + n_words + 
      (weighted_dl_distance + n_words|participant_id) + (weighted_dl_distance|item) +
      (weighted_dl_distance + n_words|language),
    data=crit %>%
      filter(
        language != "en_exp"
      ),
    control=lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e6)
    )
  )
  
  save(
    model_grammatical_linear, 
    model_distance_linear, 
    model_weighted_distance_linear,
    file="models/linear_model_fits.Rdata"
  )
}

## Extract log likelihoods
model_grammatical_linear
model_distance_linear
model_weighted_distance_linear

logLik(model_grammatical_linear) %>% as.numeric()
logLik(model_distance_linear) %>% as.numeric()
logLik(model_weighted_distance_linear) %>% as.numeric()

AIC(model_grammatical_linear)
AIC(model_distance_linear)
AIC(model_weighted_distance_linear)

BIC(model_grammatical_linear)
BIC(model_distance_linear)
BIC(model_weighted_distance_linear)


main_likelihoods <- tibble(
  model = c("grammatical", "distance", "weighted_distance"),
  log_likelihood = c(
    logLik(model_grammatical_linear) %>% as.numeric(),
    logLik(model_distance_linear) %>% as.numeric(),
    logLik(model_weighted_distance_linear) %>% as.numeric()
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
    labels=c("Grammatical", "Distance", "Weighted Distance")
  ) +
  scale_x_continuous(
    name = ""
  ) +
  guides(shape = "none") +
  theme(legend.position="top") +
  scale_y_discrete(name="")
p_main_likelihoods


### ---- Analysis By Language ----
if (file.exists("models/models_by_language_linear.RData")) {
  load("models/models_by_language_linear.RData")
} else {
  
  ## Fit a model for each language
  stats_models_by_language <- tibble(
    language = character(),
    model = character(),
    log_likelihood = double()
  )
  
  models_by_language <- list()
  for (language_iter in levels(crit$language)) {
    ## Get specific data
    crit_language_specific <- crit %>%
      filter(language == language_iter)
    
    ## Collapse to the 72 base items for English (controlled): the item id
    ## encodes the length version (en_experimental_62m), so strip the trailing
    ## s/m/l. Length then varies within item.
    if (language_iter == "en_exp") {
      crit_language_specific <- crit_language_specific %>%
        mutate(
          item = as.factor(gsub("[sml]$", "", item))
        )
      
      ## Fit models
      grammatical <- lmer(
        response_num ~ grammatical + 
          (grammatical|participant_id) + (grammatical|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
      
      distance <- lmer(
        response_num ~ dl_distance + 
          (dl_distance|participant_id) + (dl_distance|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
      
      weighted_distance <- lmer(
        response_num ~ weighted_dl_distance + n_words + 
          (weighted_dl_distance + n_words|participant_id) + 
          (weighted_dl_distance + n_words|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
    } else {
      ## Fit models
      grammatical <- lmer(
        response_num ~ grammatical + 
          (grammatical|participant_id) + (grammatical|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
      
      distance <- lmer(
        response_num ~ dl_distance + 
          (dl_distance|participant_id) + (dl_distance|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
      
      weighted_distance <- lmer(
        response_num ~ weighted_dl_distance + n_words + 
          (weighted_dl_distance + n_words|participant_id) + (weighted_dl_distance|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
    }
    
    models_by_language[[language_iter]] <- list(
      grammatical = grammatical,
      distance = distance,
      weighted_distance = weighted_distance
    )
    
    stats_models_by_language <- stats_models_by_language %>% 
      add_row(
        language = language_iter,
        model = "grammatical",
        log_likelihood = grammatical %>% logLik() %>% as.numeric()
      ) %>%
      add_row(
        language = language_iter,
        model = "distance",
        log_likelihood = distance %>% logLik() %>% as.numeric()
      ) %>%
      add_row(
        language = language_iter,
        model = "weighted_distance",
        log_likelihood = weighted_distance %>% logLik() %>% as.numeric()
      )
    
  }
  
  stats_models_by_language <- stats_models_by_language %>%
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
        levels=c("English\n(controlled)", "Russian", "Mandarin",
                 "Korean", "Hindi", "German", "French", "English", "Danish")
      ),
      model = case_when(
        model == "grammatical" ~ "Grammatical",
        model == "distance" ~ "Distance",
        model == "weighted_distance" ~ "Weighted Distance"
      ) %>%
        factor(levels=c("Grammatical", "Distance", "Weighted Distance"))
    )
  
  save(
    stats_models_by_language, 
    models_by_language,
    file="models/models_by_language_linear.RData"
  )
}

### ---- Figure based on individual-language likelihoods  
stats_models_by_language <- stats_models_by_language %>%
  group_by(language) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood
  )

p_language_likelihoods <- ggplot(stats_models_by_language) +
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
  scale_y_discrete(name="")

p_language_likelihoods


p_main_likelihoods / 
  p_language_likelihoods +
  plot_layout(heights=c(1,3))

ggsave("figures/likelihood_differences_linear.png", width=8, height=8)

### ---- DL <= 3 ----
if (file.exists("models/linear_model_fits_dl3_interactions.RData")) {
  load("models/linear_model_fits_dl3_interactions.RData")
} else {
  fit_lmer_interaction <- lmer(
    response_num ~ dl_distance * n_words +
      (dl_distance * n_words|participant_id) + (dl_distance|item) +
      (dl_distance * n_words|language),
    data=crit %>%
      filter(
        language != "en_exp",
        dl_distance <= 3
      ) %>%
      mutate(
        n_words = scale(n_words)[,1],
        dl_distance = scale(dl_distance)[,1]
      ),
    control=lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e6)
    )
  )
  
  fit_lmers_by_language <- lapply(
    levels(crit$language),
    function(lang) {
      print(paste0("Fitting for language: ", lang))
      lmer(
        response_num ~ dl_distance * n_words +
          (dl_distance * n_words|participant_id) + (dl_distance|item),
        data=crit %>%
          filter(
            language == lang,
            dl_distance <= 3
          ) %>%
          mutate(
            n_words = scale(n_words)[,1],
            dl_distance = scale(dl_distance)[,1]
          ),
        control=lmerControl(
          optimizer = "bobyqa",
          optCtrl = list(maxfun = 1e8)
        )
      )
    }
  )
  
  names(fit_lmers_by_language) <- levels(crit$language)
  
  ## Collapse the item id to the 72 base items (see Section 3): length varies
  ## within item only after the trailing s/m/l is stripped.
  fit_lmers_by_language$en_exp <- lmer(
    response_num ~ dl_distance * n_words +
      (dl_distance * n_words|participant_id) + (dl_distance * n_words|item),
    data=crit %>%
      filter(
        language == "en_exp",
        dl_distance <= 3
      ) %>%
      mutate(
        item = as.factor(gsub("[sml]$", "", item)),
        n_words = scale(n_words)[,1],
        dl_distance = scale(dl_distance)[,1]
      ),
    control=lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e8)
    )
  ) 
  
  lapply(fit_lmers_by_language, function(langfit)
    summary(langfit)$coef
         )
  
  save(
    fit_lmer_interaction,
    fit_lmers_by_language,
    file="models/linear_model_fits_dl3_interactions.RData"
  )  
}

summary(fit_lmer_interaction)


if (file.exists("models/linear_model_fits_dl3.Rdata")) {
  load("models/linear_model_fits_dl3.Rdata")
} else {
  
  model_grammatical_linear_dl3 <- lmer(
    response_num ~ grammatical + 
      (grammatical|participant_id) + (grammatical|item) + (grammatical|language),
    data=crit %>%
      filter(
        language != "en_exp",
        dl_distance <= 3
      ),
    control=lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e6)
    )
  )
  
  model_distance_linear_dl3 <- lmer(
    response_num ~ dl_distance + 
      (dl_distance|participant_id) + (dl_distance|item) + (dl_distance|language),
    data=crit %>%
      filter(
        language != "en_exp",
        dl_distance <= 3
      ),
    control=lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e6)
    )
  )
  
  model_weighted_distance_linear_dl3 <- lmer(
    response_num ~ weighted_dl_distance + n_words + 
      (weighted_dl_distance + n_words|participant_id) + (weighted_dl_distance|item) +
      (weighted_dl_distance + n_words|language),
    data=crit %>%
      filter(
        language != "en_exp",
        dl_distance <= 3
      ),
    control=lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 1e6)
    )
  )
  
  save(
    model_grammatical_linear_dl3,
    model_distance_linear_dl3,
    model_weighted_distance_linear_dl3,
    file="models/linear_model_fits_dl3.Rdata"
  )
}

## Extract log likelihoods
model_grammatical_linear_dl3
model_distance_linear_dl3
model_weighted_distance_linear_dl3

logLik(model_grammatical_linear_dl3) %>% as.numeric()
logLik(model_distance_linear_dl3) %>% as.numeric()
logLik(model_weighted_distance_linear_dl3) %>% as.numeric()

AIC(model_grammatical_linear_dl3)
AIC(model_distance_linear_dl3)
AIC(model_weighted_distance_linear_dl3)

BIC(model_grammatical_linear_dl3)
BIC(model_distance_linear_dl3)
BIC(model_weighted_distance_linear_dl3)


main_likelihoods <- tibble(
  model = c("grammatical", "distance", "weighted_distance"),
  log_likelihood = c(
    logLik(model_grammatical_linear_dl3) %>% as.numeric(),
    logLik(model_distance_linear_dl3) %>% as.numeric(),
    logLik(model_weighted_distance_linear_dl3) %>% as.numeric()
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
    labels=c("Grammatical", "Distance", "Weighted Distance")
  ) +
  scale_x_continuous(
    name = ""
  ) +
  guides(shape = "none") +
  theme(legend.position="top") +
  scale_y_discrete(name="")
p_main_likelihoods


### ---- Analysis By Language ----
if (file.exists("models/models_by_language_linear_dl3.RData")) {
  load("models/models_by_language_linear_dl3.RData")
} else {
  
  ## Fit a model for each language
  stats_models_by_language <- tibble(
    language = character(),
    model = character(),
    log_likelihood = double()
  )
  
  models_by_language <- list()
  for (language_iter in levels(crit_by_items_full$language)) {
    ## Get specific data
    crit_language_specific <- crit %>%
      filter(language == language_iter,
             dl_distance <= 3)
    
    ## Collapse to the 72 base items for English (controlled): the item id
    ## encodes the length version (en_experimental_62m), so strip the trailing
    ## s/m/l. Length then varies within item.
    if (language_iter == "en_exp") {
      crit_language_specific <- crit_language_specific %>%
        mutate(
          item = as.factor(gsub("[sml]$", "", item))
        )
      
      ## Fit models
      grammatical <- lmer(
        response_num ~ grammatical + 
          (grammatical|participant_id) + (grammatical|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
      
      distance <- lmer(
        response_num ~ dl_distance + 
          (dl_distance|participant_id) + (dl_distance|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
      
      weighted_distance <- lmer(
        response_num ~ weighted_dl_distance + n_words + 
          (weighted_dl_distance + n_words|participant_id) + 
          (weighted_dl_distance + n_words|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
    } else {
      ## Fit models
      grammatical <- lmer(
        response_num ~ grammatical + 
          (grammatical|participant_id) + (grammatical|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
      
      distance <- lmer(
        response_num ~ dl_distance + 
          (dl_distance|participant_id) + (dl_distance|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
      
      weighted_distance <- lmer(
        response_num ~ weighted_dl_distance + n_words + 
          (weighted_dl_distance + n_words|participant_id) + (weighted_dl_distance|item),
        data=crit_language_specific,
        control=lmerControl(optimizer = "bobyqa")
      )
    }
    
    models_by_language[[language_iter]] <- list(
      grammatical = grammatical,
      distance = distance,
      weighted_distance = weighted_distance
    )
    
    stats_models_by_language <- stats_models_by_language %>% 
      add_row(
        language = language_iter,
        model = "grammatical",
        log_likelihood = grammatical %>% logLik() %>% as.numeric()
      ) %>%
      add_row(
        language = language_iter,
        model = "distance",
        log_likelihood = distance %>% logLik() %>% as.numeric()
      ) %>%
      add_row(
        language = language_iter,
        model = "weighted_distance",
        log_likelihood = weighted_distance %>% logLik() %>% as.numeric()
      )
    
  }
  
  stats_models_by_language <- stats_models_by_language %>%
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
        levels=c("English\n(controlled)", "Russian", "Mandarin",
                 "Korean", "Hindi", "German", "French", "English", "Danish")
      ),
      model = case_when(
        model == "grammatical" ~ "Grammatical",
        model == "distance" ~ "Distance",
        model == "weighted_distance" ~ "Weighted Distance"
      ) %>%
        factor(levels=c("Grammatical", "Distance", "Weighted Distance"))
    )
  
  save(
    stats_models_by_language, 
    models_by_language,
    file="models/models_by_language_linear_dl3.RData"
  )
}

### ---- Figure based on individual-language likelihoods  
stats_models_by_language <- stats_models_by_language %>%
  group_by(language) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood
  )

p_language_likelihoods <- ggplot(stats_models_by_language) +
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
  scale_y_discrete(name="")

p_language_likelihoods


p_main_likelihoods / 
  p_language_likelihoods +
  plot_layout(heights=c(1,3))

ggsave("figures/likelihood_differences_linear_dl3.png", width=8, height=8)


## ---- Main Analysis with Maximal Exclusions ----
crit_exclusions <- crit %>%
  filter(
    first_language %in% c(
      "Dansk",
      "Deutsch",
      "हिंदी",
      "한국어",
      "Русский"
    ) |
      english_first_language %in% c("Yes", "是") |
      french_first_language == "Oui"
  )

table(crit_exclusions$language)

crit_exclusions_by_items <- crit_exclusions %>% 
  mutate(
    weighted_dl_distance = crit_exclusions$dl_distance / crit_exclusions$n_words,
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
           grammatical, sum_surprisal, condition, weighted_dl_distance_surprisal) %>%
  summarize(
    response_num = mean(response_num)
  ) %>%
  mutate(
    grammatical_numeric = ifelse(grammatical == "grammatical", .5, -.5)
  )

### ---- Across Languages ----
if (file.exists("models/main_model_fits_maximal_exclusions.RData")) {
  load("models/main_model_fits_maximal_exclusions.RData")
} else {
  
  mod_grammatical <- scam(
    response_num ~ grammatical_numeric +# + # Fixed effects
      s(language, bs="re") + s(language, grammatical_numeric, bs="re"),# RE
    data=crit_exclusions_by_items %>% 
      filter(language != "en_exp"),
    control = scam.control(trace = TRUE)
  )
  
  mod_distance <- scam(
    response_num ~ s(dl_distance, bs="mpd") + # Fixed effects
      s(language, bs="re") + s(language, dl_distance, bs="re"), # RE
    data=crit_exclusions_by_items %>% 
      filter(language != "en_exp"),
    control = scam.control(trace = TRUE)
  )
  
  mod_weighted_distance <- scam(
    response_num ~ s(weighted_dl_distance, bs="mpd") + s(n_words, bs="mpd") + # Fixed effects
      s(language, bs="re") + s(language, weighted_dl_distance, bs="re") + s(language, n_words, bs="re"), # RE
    data=crit_exclusions_by_items %>% 
      filter(language != "en_exp"),
    control = scam.control(trace = TRUE)
  )
  
  save(
    mod_grammatical, 
    mod_distance, 
    mod_weighted_distance,
    file="models/main_model_fits_maximal_exclusions.RData"
  )
}

## Inspect models
mod_grammatical$conv
mod_distance$conv
mod_weighted_distance$conv

## Extract the marginal log-likelihood of each model (see the note above: scam's
## own logLik()/AIC()/BIC() are conditional on the fitted random effects).
marginal_logLik(mod_grammatical)
marginal_logLik(mod_distance)
marginal_logLik(mod_weighted_distance)

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
    labels=c("Grammatical", "Distance", "Weighted Distance")
  ) +
  scale_x_continuous(
    name = ""
  ) +
  guides(shape = "none") +
  theme(legend.position="top") +
  scale_y_discrete(name="")
p_main_likelihoods

### ---- Analysis By Language: Length ----
if (file.exists("models/models_by_language_maximal_exclusions.RData")) {
  load("models/models_by_language_maximal_exclusions.RData")
} else {
  
  ## Fit a model for each language
  stats_models_by_language <- tibble(
    language = character(),
    model = character(),
    log_likelihood = double()
  )
  
  models_by_language <- list()
  for (language_iter in levels(crit$language)) {
    print(
      paste0("Fitting models for ", language_iter)
    )
    ## Get specific data
    crit_language_specific <- crit_exclusions_by_items %>%
      filter(language == language_iter)
    
    ## Fit models
    grammatical <- scam(
      response_num ~ grammatical,
      data=crit_language_specific
    )
    
    distance <- scam(
      response_num ~ s(dl_distance, bs="mpd"),
      data=crit_language_specific
    )
    
    weighted_distance <- scam(
      response_num ~ s(weighted_dl_distance, bs="mpd") + s(n_words, bs="mpd"),
      data=crit_language_specific
    )

    models_by_language[[language_iter]] <- list(
      grammatical = grammatical,
      distance = distance,
      weighted_distance = weighted_distance
    )

    ## Unlike the other by-language blocks, this tibble is NOT rebuilt after the
    ## cache guard, so it must carry the marginal value directly. These models
    ## carry no random effects, so marginal_logLik() returns the same number as
    ## logLik() -- but using it keeps the figure correct if REs are ever added.
    stats_models_by_language <- stats_models_by_language %>%
      add_row(
        language = language_iter,
        model = "grammatical",
        log_likelihood = marginal_logLik(grammatical)
      ) %>%
      add_row(
        language = language_iter,
        model = "distance",
        log_likelihood = marginal_logLik(distance)
      ) %>%
      add_row(
        language = language_iter,
        model = "weighted_distance",
        log_likelihood = marginal_logLik(weighted_distance)
      )
    
  }
  
  stats_models_by_language <- stats_models_by_language %>%
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
        levels=c("English\n(controlled)", "Russian", "Mandarin",
                 "Korean", "Hindi", "German", "French", "English", "Danish")
      ),
      model = case_when(
        model == "grammatical" ~ "Grammatical",
        model == "distance" ~ "Distance",
        model == "weighted_distance" ~ "Weighted Distance"
      ) %>%
        factor(levels=c("Grammatical", "Distance", "Weighted Distance"))
    )
  
  save(
    stats_models_by_language, 
    models_by_language,
    file="models/models_by_language_maximal_exclusions.RData"
  )
}

stats_models_by_language <- stats_models_by_language %>%
  group_by(language) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood
  )

p_language_likelihoods <- ggplot(stats_models_by_language) +
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
  scale_y_discrete(name="")

p_language_likelihoods


p_main_likelihoods / 
  p_language_likelihoods +
  plot_layout(heights=c(1,3))

ggsave("figures/likelihood_differences_max_exclusions.png", width=8, height=8)


## ---- Figures ----
#### ---- Helper Functions 
YLIM <- c(1, 7)

DL_COLORS <- c("0"="#EB7AA3", "1-2"="#1E88E5", "3"="#FFC107", "4+"="#BBBBBB")
DL_BUCKET_LEVELS <- c("0", "1-2", "3", "4+")

bucket_dl <- function(dl_distance) {
  dplyr::case_when(
    dl_distance == 0 ~ "0",
    dl_distance %in% 1:2 ~ "1-2",
    dl_distance == 3 ~ "3",
    dl_distance >= 4 ~ "4+"
  ) %>%
    factor(levels = DL_BUCKET_LEVELS)
}

scale_dl_color <- function() {
  ggplot2::scale_color_manual(
    name   = "Damerau-Levenshtein Distance",
    values = DL_COLORS,
    breaks = DL_BUCKET_LEVELS,
    limits = DL_BUCKET_LEVELS,
    drop   = FALSE,
    guide  = ggplot2::guide_legend(
      override.aes = list(alpha = 1, linewidth = 0, size = 3)
    )
  )
}


get_prediction_grid_scam <- function(
    lang,
    model_type,
    n_words=4:43, 
    dl_distance=c(0, 1, 3)
) {
  expand.grid(
    n_words = n_words,
    dl_distance = dl_distance
  ) %>%
    mutate(
      weighted_dl_distance = dl_distance / n_words,
      grammatical_numeric = ifelse(dl_distance == 0, .5, -.5),
      item=factor(models_by_language[[lang]][[model_type]]$model$item[1], levels=levels(models_by_language[[lang]][[model_type]]$model$item))
    )
}

get_prediction_grid_lmer <- function(
    lang,
    model_type,
    n_words=4:43, 
    dl_distance=c(0, 1, 3)
) {
  expand.grid(
    n_words = n_words,
    dl_distance = dl_distance
  ) %>%
    mutate(
      weighted_dl_distance = dl_distance / n_words,
      grammatical_numeric = ifelse(dl_distance == 0, .5, -.5),
      item=factor(models_by_language[[lang]][[model_type]]@flist[1]$item[1], levels=levels(models_by_language[[lang]][[model_type]]@flist[1]$item)),
      participant_id=factor(models_by_language[[lang]][[model_type]]@flist[2]$participant_id[1], levels=levels(models_by_language[[lang]][[model_type]]@flist[2]$participant_id))
    )
}

### ---- Baseline data for plotting (item means) ----
plot_data_by_item <- crit %>% 
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


### ---- Model-based plots ----
#### ---- lmer predictions ----
##### ----- Communication, full data ----
load("models/models_by_language_linear.RData")

## Communication
fulldata_communication_predictions <- setNames(
  lapply(names(models_by_language), function(lang) {
    predict(
      models_by_language[[lang]]$weighted_distance,
      newdata = get_prediction_grid_lmer(lang, "weighted_distance"),
      type = "response",
      se.fit = TRUE,
      re.form = NA
    )
  }
  ),
  names(models_by_language)
)

fulldata_communication_plot_data <- bind_rows(
  lapply(names(fulldata_communication_predictions), function(lang) {
    preds <- fulldata_communication_predictions[[lang]]
    grid  <- get_prediction_grid_lmer(lang, "weighted_distance")
    
    tibble(
      language = label_language(lang),
      n_words = grid$n_words,
      dl_distance_bucket = bucket_dl(grid$dl_distance),
      response = preds$fit,
      se = preds$se.fit
    )
  })
)

p_full_communication <- p_baseline +
  geom_ribbon(
    data=fulldata_communication_plot_data %>%
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
    data=fulldata_communication_plot_data %>%
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
    legend.position = "top",
    plot.title = element_text(hjust = 0.5),
    
    # put facet strips on the outside + rotate text
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 90)
  ) +
  ggtitle("Communication Model") +
  scale_y_continuous(
    limits = YLIM,
    breaks = seq(YLIM[1], YLIM[2], by=2),
    labels = seq(YLIM[1], YLIM[2], by=2)
  ) +
  ggplot2::guides(color = "none")

p_full_communication


## Gradient
fulldata_gradient_predictions <- setNames(
  lapply(names(models_by_language), function(lang) {
    predict(
      models_by_language[[lang]]$distance,
      newdata = get_prediction_grid_lmer(lang, "distance"),
      type = "response",
      se.fit = TRUE,
      re.form = NA
    )
  }
  ),
  names(models_by_language)
)

fulldata_gradient_predictions_plot_data <- bind_rows(
  lapply(names(fulldata_gradient_predictions), function(lang) {
    preds <- fulldata_gradient_predictions[[lang]]
    grid  <- get_prediction_grid_lmer(lang, "distance")
    
    tibble(
      language = label_language(lang),
      n_words = grid$n_words,
      dl_distance_bucket = bucket_dl(grid$dl_distance),
      response = preds$fit,
      se = preds$se.fit
    )
  })
)

p_full_gradient <- p_baseline +
  geom_ribbon(
    data=fulldata_gradient_predictions_plot_data %>%
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
    data=fulldata_gradient_predictions_plot_data %>%
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
  ggtitle("Gradient Model") +
  scale_y_continuous(
    limits = YLIM,
    breaks = NULL,
    labels = NULL
  ) +
  ggplot2::guides(color = "none")

p_full_gradient

## with_shuffled
load("models/models_by_language_linear_dl3.RData")
with_shuffled_communication_predictions <- setNames(
  lapply(names(models_by_language), function(lang) {
    predict(
      models_by_language[[lang]]$weighted_distance,
      newdata = get_prediction_grid_lmer(lang, "weighted_distance"),
      type = "response",
      se.fit = TRUE,
      re.form = NA
    )
  }
  ),
  names(models_by_language)
)

with_shuffled_communication_plot_data <- bind_rows(
  lapply(names(with_shuffled_communication_predictions), function(lang) {
    preds <- with_shuffled_communication_predictions[[lang]]
    grid  <- get_prediction_grid_lmer(lang, "weighted_distance")
    
    tibble(
      language = label_language(lang),
      n_words = grid$n_words,
      dl_distance_bucket = bucket_dl(grid$dl_distance),
      response = preds$fit,
      se = preds$se.fit
    )
  })
)

p_with_shuffled_communication <- p_baseline +
  geom_ribbon(
    data=with_shuffled_communication_plot_data,
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
    data=with_shuffled_communication_plot_data,
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
  ggtitle("Communication Model (DL ≤ 3)") +
  scale_y_continuous(
    limits = YLIM,
    breaks = NULL,
    labels = NULL
  ) +
  ggplot2::guides(color = "none")

p_with_shuffled_communication

## Gradient
with_shuffled_gradient_predictions <- setNames(
  lapply(names(models_by_language), function(lang) {
    predict(
      models_by_language[[lang]]$distance,
      newdata = get_prediction_grid_lmer(lang, "distance"),
      type = "response",
      se.fit = TRUE,
      re.form = NA
    )
  }
  ),
  names(models_by_language)
)

with_shuffled_gradient_plot_data <- bind_rows(
  lapply(names(with_shuffled_gradient_predictions), function(lang) {
    preds <- with_shuffled_gradient_predictions[[lang]]
    grid  <- get_prediction_grid_lmer(lang, "distance")
    
    tibble(
      language = label_language(lang),
      n_words = grid$n_words,
      dl_distance_bucket = bucket_dl(grid$dl_distance),
      response = preds$fit,
      se = preds$se.fit
    )
  })
)

p_with_shuffled_gradient <- p_baseline +
  geom_ribbon(
    data=with_shuffled_gradient_plot_data,
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
    data=with_shuffled_gradient_plot_data,
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
  ggtitle("Gradient Model (DL ≤ 3)") +
  scale_y_continuous(
    limits = YLIM,
    breaks = NULL,
    labels = NULL
  ) +
  ggplot2::guides(color = "none")

p_with_shuffled_gradient

## Combine and save linear figures
linear_combined_plot <- (p_full_communication |
  p_full_gradient |
  p_with_shuffled_communication |
  p_with_shuffled_gradient) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
  )

## Save as pdf
ggsave(
  filename = "figures/main_figure_linear_models.pdf",
  plot     = linear_combined_plot,
  width    = 25,
  height   = 35,
  units    = "in"
)

