library(tidyverse)
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

## ---- Marginal log-likelihood for scam/gam models with random effects ----
## scam/gam's logLik() is CONDITIONAL on the fitted random effects b; this returns
## the MARGINAL log-likelihood (integrating the b's out), i.e. the quantity that
## lmer() reports, so that model comparison does not credit a model for its fitted
## item- or language-specific random effects. Copied unchanged from
## analysis/main_analysis.R so that the numbers here are directly comparable.
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

## ---- Read and Prepare Data ----
## Read full data
dat <- read_csv("../data/full_data.csv")

## Examine comprehension questions
comprehension <- dat %>% filter(grepl("comprehension", experiment))

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
  filter(experiment == "exp_b")

## Remove nonsense conditions
crit <- crit %>% 
  filter(!grepl("non_sense", condition))

## Convert responses to numbers
crit$response <- as.numeric(substr(crit$response,1,1))

## ---- Add the SLOR and MORCELA Scores ----
## get_slor_morcela_scores.py wrote one row per sentence per unigram model. The
## item ids there are bare numbers, while full_data.csv prefixes them with the
## folder name, so rebuild the same prefix before joining. Note that the
## controlled English experiment is called en_exp in the data but lives in the
## en_experimental folder.
scores <- read_csv("../data/slor_morcela_scores.csv") %>%
  mutate(
    item = paste(
      ifelse(language == "en_exp", "en_experimental", language),
      item,
      sep = "_"
    )
  )

## We report the model-initial unigram throughout. wordfreq is kept in the file
## as a robustness check: it gives the same answer everywhere, but it does not
## know 9% of the Danish word tokens or 10% of the Korean ones, so it cannot be
## the default in a cross-linguistic analysis.
UNIGRAM_MODEL <- "model_initial"

crit <- crit %>%
  left_join(
    scores %>%
      filter(unigram_model == UNIGRAM_MODEL) %>%
      select(language, item, condition, log_p_sentence, log_p_unigram),
    by = c("language", "item", "condition")
  )

## ---- Aggregate to item means ----
## The same aggregation main_analysis.R uses for its scam models: one row per
## item and condition, ratings averaged over participants, distances above 3
## dropped because those materials sit at floor.
crit_by_items <- crit %>% 
  mutate(
    weighted_dl_distance = dl_distance / n_words,
    grammatical = factor(
      ifelse(dl_distance == 0, "grammatical", "ungrammatical"),
      levels=c("ungrammatical", "grammatical")),
    item = as.factor(item),
    language = as.factor(language),
    response_num = response
  ) %>%
  group_by(language, item, n_words, dl_distance, weighted_dl_distance, 
           grammatical, condition, log_p_sentence, log_p_unigram) %>%
  summarize(
    response_num = mean(response_num)
  ) %>%
  mutate(
    grammatical_numeric = ifelse(grammatical == "grammatical", .5, -.5)
  ) %>%
  ungroup() %>%
  filter(
    dl_distance <= 3
  )

## ---- Compute SLOR ----
## SLOR (Lau, Clark & Lappin 2017) is the sentence's log probability, corrected
## for how frequent its words are, divided by how long it is:
##
##     SLOR = ( log p_model(sentence) - log p_unigram(sentence) ) / length
##
## Subtracting the unigram term removes the penalty a sentence pays simply for
## containing rare words; dividing by length removes the penalty it pays for
## being long. We measure length in words, so that the denominator is the same
## n_words the communication-based model uses.
crit_by_items <- crit_by_items %>%
  mutate(
    slor = (log_p_sentence - log_p_unigram) / n_words
  )

## ---- Compute MORCELA ----
## MORCELA (Tjuatja, Neubig, Linzen & Hao 2025) is SLOR with the two corrections
## allowed to be weaker or stronger than SLOR assumes:
##
##     MORCELA = ( log p_model(sentence) - beta * log p_unigram(sentence) + gamma ) / length
##
## SLOR is the special case beta = 1, gamma = 0. The two weights are learned from
## the acceptability ratings themselves, by the recipe in that paper: regress the
## ratings on the three terms below and read the weights off the coefficients.
##
## Note that this fits beta and gamma on the very ratings the score is then used
## to predict. That is what the MORCELA paper does, and it is deliberately
## generous to MORCELA -- we want to compare against the strongest version of it.
morcela_score <- function(model_data) {
  morcela_regression <- lm(
    response_num ~ log_p_sentence_over_length + log_p_unigram_over_length + one_over_length,
    data = model_data
  )
  weights <- coef(morcela_regression)
  beta  <- -weights[["log_p_unigram_over_length"]] / weights[["log_p_sentence_over_length"]]
  gamma <-  weights[["one_over_length"]]           / weights[["log_p_sentence_over_length"]]

  model_data$log_p_sentence_over_length -
    beta * model_data$log_p_unigram_over_length +
    gamma * model_data$one_over_length
}

crit_by_items <- crit_by_items %>%
  mutate(
    log_p_sentence_over_length = log_p_sentence / n_words,
    log_p_unigram_over_length  = log_p_unigram / n_words,
    one_over_length            = 1 / n_words
  )

## ---- Prepare One Language's Data ----
## Used both when fitting below and when plotting further down, so that the two
## can never drift apart. MORCELA's weights are learned separately for each
## language, since there is no reason a single beta and gamma should suit all of
## them.
language_data <- function(language_iter) {
  language_specific <- crit_by_items %>%
    filter(language == language_iter) %>% droplevels()

  language_specific$morcela <- morcela_score(language_specific)

  if (language_iter == "en_exp") {
    ## In this experiment each of the 72 base items was written in 3 lengths, so
    ## the item id encodes the length version (e.g. en_experimental_62m).
    ## Stripping the trailing s/m/l recovers the 72 base items, and length then
    ## varies WITHIN item, which is what makes the by-item slopes identifiable.
    language_specific <- language_specific %>%
      mutate(
        item = as.factor(gsub("[slm]$", "", item))
      )
  }

  language_specific
}

## ---- Analysis By Language ----
if (file.exists("models/models_by_language_slor_morcela.RData")) {
  load("models/models_by_language_slor_morcela.RData")
} else {

  ## Fit a model for each language. The three models from the paper are fit
  ## exactly as in analysis/main_analysis.R, so that their panels in the figure
  ## below are identical to the published ones. SLOR and MORCELA are given the
  ## same random-effects structure as those three.
  ##
  ## Both SLOR and MORCELA get a monotonically INCREASING smooth: a higher score
  ## means the sentence is more probable than its words alone would predict, and
  ## so should be judged more acceptable. We checked this rather than assumed it
  ## -- fitting the decreasing version instead collapses the smooth to a flat
  ## line in every language, because the relationship really does go up.
  models_by_language <- list()
  for (language_iter in levels(crit_by_items$language)) {
    print(
      paste0("Fitting models for ", language_iter)
    )

    ## Get specific data
    crit_by_items_language_specific <- language_data(language_iter)

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

    models_by_language[[language_iter]] <- list(
      grammatical = grammatical,
      distance = distance,
      weighted_distance = weighted_distance,
      slor = slor,
      morcela = morcela
    )

  }

  save(
    models_by_language,
    file="models/models_by_language_slor_morcela.RData"
  )
}

## Build the by-language log-likelihood table using the MARGINAL log-likelihood
## (integrating the item random effects out) rather than the conditional value
## logLik() returns. Done here, after the cache block, so it applies whether the
## models were just fit or loaded from disk.
stats_models_by_language <- do.call(rbind, lapply(names(models_by_language), function(lg) {
  M <- models_by_language[[lg]]
  tibble(language = lg, model = names(M),
         log_likelihood = vapply(M, marginal_logLik, numeric(1)))
})) %>%
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
      model == "grammatical" ~ "Binary",
      model == "distance" ~ "Gradient",
      model == "weighted_distance" ~ "Communication",
      model == "slor" ~ "SLOR",
      model == "morcela" ~ "MORCELA"
    ) %>%
      factor(levels=c("Communication", "Gradient", "Binary", "SLOR", "MORCELA"))
  )

write_csv(stats_models_by_language, "figures/likelihoods_by_language.csv")

### ---- Figure based on individual-language likelihoods ----
stats_models_by_language <- stats_models_by_language %>%
  group_by(language) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood
  )

## Filled shapes (21 to 25 take both a colour and a fill), so every model is a
## solid marker. The three models from the paper keep the colours they have
## there, and their shapes are the filled versions of the defaults ggplot gave
## them: Binary a circle, Gradient a triangle, Communication a square.
MODEL_COLORS <- c(
  "Communication" = "#004D40",
  "Gradient"      = "#1E88E5",
  "Binary"        = "#EB7AA3",
  "SLOR"          = "#E08A2E",
  "MORCELA"       = "#A31F34"
)
MODEL_SHAPES <- c(
  "Communication" = 22,
  "Gradient"      = 24,
  "Binary"        = 21,
  "SLOR"          = 23,
  "MORCELA"       = 25
)

p_language_likelihoods <- ggplot(stats_models_by_language) +
  geom_point(
    aes(
      y=language,
      x=likelihood_diff,
      color=model,
      fill=model,
      shape=model
    ),
    size=4
  ) +
  scale_color_manual(name="", values=MODEL_COLORS) +
  scale_fill_manual(name="", values=MODEL_COLORS) +
  scale_shape_manual(name="", values=MODEL_SHAPES) +
  scale_x_continuous(
    name = "Log-Likelihood Difference"
  ) +
  theme(legend.position="top") +
  scale_y_discrete(name="")

ggsave("figures/likelihood_differences_slor_morcela.png",
       plot = p_language_likelihoods, width = 10, height = 7)

## ---- Main figure ----
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
  levels=c("Danish", "English", "French", "German", "Hindi", "Korean", 
           "Mandarin", "Russian", "English\n(controlled)")
  )
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

## The three models from the paper are functions of length and of the number of
## errors, so we can ask them what they predict for any combination of the two,
## whether or not a sentence like that was in the materials. That is what this
## grid is for, and it is the same grid main_analysis.R uses.
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
      item=factor(models_by_language[[lang]][[model_type]]$model$item[1], levels=levels(models_by_language[[lang]][[model_type]]$model$item))
    )
}

### ---- Baseline data for plotting (item means) ----
## Same aggregation as crit_by_items above, but deliberately WITHOUT the
## filter(dl_distance <= 3): the scatter shows the shuffled items as well, in the
## grey "4+" bucket, as context. The models (and hence the fitted curves drawn
## over these points) are still fit on distances up to 3 only.
plot_data_by_item <- crit %>%
  mutate(
    item = as.factor(item),
    language = as.factor(language),
    response_num = response
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

### ---- Model-based plots ----
## One function for the three models from the paper, since they differ only in
## which model object is asked for predictions.
grid_model_plot_data <- function(model_type) {
  bind_rows(
    lapply(names(models_by_language), function(lang) {
      grid <- get_prediction_grid(
        lang, model_type,
        n_words = range(crit[crit$language == lang,]$n_words)
      )
      preds <- predict(
        models_by_language[[lang]][[model_type]],
        newdata = grid,
        type = "response",
        se.fit = TRUE,
        exclude = "s(item)"
      )
      tibble(
        language = label_language(lang),
        n_words = grid$n_words,
        dl_distance = grid$dl_distance,
        dl_distance_bucket = bucket_dl(grid$dl_distance),
        response = preds$fit,
        se = preds$se.fit
      )
    })
  )
}

## The three models above give one smooth curve each because they are smooth
## functions of length to begin with. SLOR and MORCELA give one prediction per
## sentence, and sentences of the same length differ, so their raw predictions
## zig-zag. We therefore smooth them over length in exactly the way the Human
## Judgments panel smooths the ratings themselves, so that the six columns are
## drawn on the same terms.
smooth_over_length <- function(model_predictions) {
  if (n_distinct(model_predictions$n_words) < 4 ||
      sd(model_predictions$response) < 1e-8) {
    return(tibble(
      n_words = range(model_predictions$n_words),
      response = mean(model_predictions$response),
      se = mean(model_predictions$se)
    ))
  }
  smooth_fit <- gam(response ~ s(n_words, bs = "tp"), data = model_predictions)
  lengths_to_plot <- seq(
    min(model_predictions$n_words),
    max(model_predictions$n_words),
    length.out = 120
  )
  smooth_predictions <- predict(
    smooth_fit,
    newdata = data.frame(n_words = lengths_to_plot),
    se.fit = TRUE
  )
  tibble(
    n_words = lengths_to_plot,
    response = as.numeric(smooth_predictions$fit),
    se = as.numeric(smooth_predictions$se.fit)
  )
}

## SLOR and MORCELA are NOT functions of length and error count: many different
## sentences share a length and a number of errors, and they have quite different
## scores. So there is no grid to ask them about. Instead we ask them about the
## sentences that actually exist, and average their answers at each length. This
## is the honest version of the same question -- what these models predict, on
## average, for real sentences of a given length and error count.
observed_model_plot_data <- function(model_type) {
  bind_rows(
    lapply(names(models_by_language), function(lang) {
      model_data <- language_data(lang)
      preds <- predict(
        models_by_language[[lang]][[model_type]],
        newdata = model_data,
        type = "response",
        se.fit = TRUE,
        exclude = "s(item)"
      )
      tibble(
        language = label_language(lang),
        n_words = model_data$n_words,
        dl_distance = model_data$dl_distance,
        dl_distance_bucket = bucket_dl(model_data$dl_distance),
        response = preds$fit,
        se = preds$se.fit
      ) %>%
        filter(dl_distance %in% c(0, 1, 3)) %>%
        group_by(language, dl_distance, dl_distance_bucket) %>%
        group_modify(~ smooth_over_length(.x)) %>%
        ungroup()
    })
  )
}

## One function to draw a model panel, so the five look alike. The colours
## argument exists for the binary model: it predicts the same rating for every
## ungrammatical sentence, so drawing its 1-2 and 3 curves in different colours
## would hide one line underneath the other. Painting them the same colour shows
## honestly that the model makes one prediction there.
model_panel <- function(plot_data, panel_title, colours = NULL) {
  p <- p_baseline +
    ggnewscale::new_scale_color() +
    geom_ribbon(
      data=plot_data,
      aes(
        x = n_words,
        ymin = response - 2 * se,
        ymax = response + 2 * se,
        group = dl_distance
      ),
      alpha = 0.25,
      color = NA
    ) +
    geom_line(
      data=plot_data,
      aes(
        x = n_words,
        y = response,
        color = dl_distance_bucket,
        group = dl_distance
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
    ggtitle(panel_title) +
    scale_y_continuous(
      limits = YLIM,
      breaks = NULL,
      labels = NULL
    ) +
    ggplot2::guides(color = "none")

  if (is.null(colours)) {
    p + scale_color_manual(values = DL_COLORS)
  } else {
    p + scale_color_manual(values = colours)
  }
}

p_grammatical <- model_panel(
  grid_model_plot_data("grammatical"),
  "Binary Grammaticality Model",
  colours = c("0" = "#EB7AA3", "1-2" = "#06402B", "3" = "#06402B", "4+" = "#06402B")
)

p_distance <- model_panel(
  grid_model_plot_data("distance"),
  "Gradient Grammaticality Model"
)

p_weighted_distance <- model_panel(
  grid_model_plot_data("weighted_distance") %>%
    mutate(
      response = ifelse(response > 7.5, 7.5, response)
    ),
  "Communication Model"
)

p_slor <- model_panel(
  observed_model_plot_data("slor"),
  "SLOR"
)

p_morcela <- model_panel(
  observed_model_plot_data("morcela"),
  "MORCELA"
)

### ---- Combine plots into main figure ----
main_figure = (p_gam | p_morcela | p_slor | p_weighted_distance | p_distance | p_grammatical) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)
        )

ggsave(
  filename = "figures/main_figure_slor_morcela.png",
  plot     = main_figure,
  width    = 28,
  height   = 20,
  units    = "in",
  dpi      = 200
)

## Save as pdf
ggsave(
  filename = "figures/main_figure_slor_morcela.pdf",
  plot     = main_figure,
  width    = 28,
  height   = 20,
  units    = "in"
)

print("Done.")
