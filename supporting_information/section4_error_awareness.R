## =============================================================================
## SI Section 4: Error awareness (English)
##
## The English experiment in which participants also judged whether each
## sentence contained an error. Produces the four error_awareness_*.png figures.
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

## A separate English experiment with 80 participants, on the same materials as
## the original English experiment. On every trial participants both rated the
## sentence's acceptability (exp_b) and answered "Are there any errors in this
## sentence?" (exp_b_comprehension). `literal_response` holds the correct Yes/No
## answer -- "No" for the original condition, "Yes" once any adjacent-token
## exchanges were introduced -- so a trial counts as correctly judged exactly
## when response == literal_response.
##
## Trial order: the error question came FIRST and the rating second. In the
## Surveyor materials the error question is the anchor trial (it carries the
## sentence) and the rating is its coordinated follow-up, so the rating always
## landed on the next screen -- verified in the data, where the rating's
## presentation_order is exactly one greater on all 5,686 trials. The original
## English experiment is the other way round: it rated first and asked its
## meaning-comprehension question second. So the two experiments differ in the
## order of the two questions as well as in what the second one asked, which the
## 4.3 and 4.4/4.5 comparisons below have to be read in light of.
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

## ---- Shared setup ----
## Colours matching the main-text DL-distance buckets (0 / 1-2 / 3 / 4+), mapped
## to the four error conditions so this section matches the rest of the study.
error_colors <- c("Original" = "#EB7AA3", "1 swap" = "#1E88E5",
                  "3 swaps" = "#FFC107", "Shuffled" = "#BBBBBB")

condition_labels <- function(x) factor(
  x, levels = c("original", "swap_1", "swap_3", "random"),
  labels = c("Original", "1 swap", "3 swaps", "Shuffled"))

## Every figure in this section plots item means, so collapse them in one place.
item_means <- function(d, value, ...) {
  d %>%
    group_by(...) %>%
    summarize(value = mean({{ value }}), n_trials = n(), .groups = "drop")
}

## Predictors are z-scored inside each fitting subset: it helps the maximal
## random-effect structures converge, and a linear rescaling of a predictor
## leaves the log-likelihood unchanged (the random-effect covariances are
## unstructured), so the 4.4/4.5 comparisons are unaffected.
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
## `formulas` only if the optimizer actually fails to converge on it.
fit_maximal <- function(formulas, fitter, label) {
  for (i in seq_along(formulas)) {
    fit <- suppressWarnings(fitter(as.formula(formulas[[i]])))
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


## ---- The two error-awareness experiments ----
## The same design was run twice: once on the original English materials, where
## sentence length varies BETWEEN items, and once on the length-controlled
## materials, where each of the 72 base items exists in a short, medium and long
## version so length varies WITHIN item. That difference is what `slopes_by_item`
## records -- it is the only thing that changes in the random-effect structures.
##
## Both are prepared to the same shape: one row per participant x sentence,
## carrying the acceptability rating, that trial's Yes/No error judgment, and the
## stimulus properties joined in from the corresponding original experiment.

## English: DL distance and length come from `crit`, which already applies the
## >= 80% comprehension exclusion. Item ids are "en_17" there and "17" here.
prepare_en <- function() {
  key <- crit %>%
    filter(language == "en") %>%
    distinct(item, condition, dl_distance, n_words) %>%
    mutate(item = as.integer(str_extract(as.character(item), "[0-9]+"))) %>%
    rename(cond_family = condition)

  read_csv("../data/en_error_awareness/data_en_error_awareness.csv",
           show_col_types = FALSE) %>%
    mutate(item = as.integer(item), cond_family = cond,
           item_version = as.character(item)) %>%
    left_join(key, by = c("item", "cond_family"))
}

## Controlled: the conditions carry the length ("swap-3_medium"), and the
## stimulus properties come from the controlled experiment's stimulus table.
## That table has all 864 (item, condition) cells; full_data.csv is response data
## and only covers the cells its own participants happened to see.
prepare_en_exp <- function() {
  key <- read_csv("../data/en_experimental/stimuli_data.csv", show_col_types = FALSE) %>%
    transmute(item_version = item, cond_family = condition, dl_distance, n_words)

  read_csv("../data/en_experimental_error_awareness/data_en_exp_error_awareness.csv",
           show_col_types = FALSE) %>%
    mutate(
      item = as.integer(item),
      sentence_length = str_extract(cond, "(short|medium|long)$"),
      cond_family = str_replace(str_remove(cond, "_(short|medium|long)$"), "-", "_"),
      item_version = paste0(item, case_when(sentence_length == "short"  ~ "s",
                                            sentence_length == "medium" ~ "m",
                                            sentence_length == "long"   ~ "l"))
    ) %>%
    left_join(key, by = c("item_version", "cond_family"))
}

## Join each rating to the same sentence's error judgment. `literal_response` is
## the correct Yes/No answer -- "No" for the original condition, "Yes" once any
## adjacent-token exchanges were introduced -- so a trial counts as correctly
## judged exactly when response == literal_response.
build_ea_trials <- function(raw) {
  raw %>%
    filter(experiment == "exp_b") %>%
    transmute(participant_id, item, item_version, cond_family, dl_distance, n_words,
              acceptability = as.numeric(substr(response, 1, 1))) %>%
    inner_join(
      raw %>%
        filter(experiment == "exp_b_comprehension") %>%
        transmute(participant_id, item, item_version, cond_family,
                  said_error = response == "Yes",
                  detected = response == literal_response),
      by = c("participant_id", "item", "item_version", "cond_family")) %>%
    mutate(
      condition_label = condition_labels(cond_family),
      detection_label = factor(
        ifelse(detected, "Error judged correctly", "Error judged incorrectly"),
        levels = c("Error judged correctly", "Error judged incorrectly")),
      weighted_dl_distance = dl_distance / n_words,
      grammatical = factor(ifelse(dl_distance == 0, "grammatical", "ungrammatical"),
                           levels = c("ungrammatical", "grammatical")))
}

## The comprehension-question experiment each one is compared against in 4.3.
reference_trials <- function(language_code) {
  crit %>%
    filter(language == language_code) %>%
    transmute(item_version = str_remove(as.character(item), "^en(_experimental)?_"),
              cond_family = condition, n_words, dl_distance,
              acceptability = response_num,
              condition_label = condition_labels(condition))
}

ea_experiments <- list(
  en = list(
    label = "English", slug = "en", reference = "en",
    prepare = prepare_en, slopes_by_item = FALSE),
  en_exp = list(
    label = "English (controlled)", slug = "en_exp", reference = "en_exp",
    prepare = prepare_en_exp, slopes_by_item = TRUE)
)

ea_likelihoods_all <- NULL
ea_trials_all <- NULL

for (spec in ea_experiments) {
  cat("\n\n#############################################################\n")
  cat("## Error awareness --", spec$label, "\n")
  cat("#############################################################\n")

  trials <- build_ea_trials(spec$prepare())
  ea_trials_all <- bind_rows(ea_trials_all,
                             trials %>% mutate(experiment_label = spec$label))

  ## The models restrict to dl_distance <= 3, as the main analyses do.
  model_trials <- trials %>%
    filter(dl_distance <= 3) %>%
    mutate(participant_id = factor(participant_id), item = factor(item),
           item_version = factor(item_version)) %>%
    droplevels()

  ## Length is between items in the English experiment and within item in the
  ## controlled one, so only the latter supports a by-item length slope.
  item_slope <- function(terms) {
    terms <- if (spec$slopes_by_item) terms else setdiff(terms, "n_words_z")
    paste0("(", paste(terms, collapse = " + "), " | item)")
  }


  ## ---- 4.1 Error detection: descriptives, figure, and model ----
  cat("\n===== 4.1 Error detection --", spec$label, "=====\n")

  trials %>%
    summarize(participants = n_distinct(participant_id), items = n_distinct(item),
              sentences = n_distinct(item_version), trials = n()) %>%
    print()

  trials %>%
    group_by(condition_label) %>%
    summarize(n = n(), accuracy = mean(detected), said_error = mean(said_error),
              mean_rating = mean(acceptability), .groups = "drop") %>%
    print(n = Inf)

  trials %>%
    summarize(overall_accuracy = mean(detected),
              false_alarm_rate = mean(said_error[cond_family == "original"]),
              hit_rate = mean(said_error[cond_family != "original"])) %>%
    print()

  ## Figure: by-item detection accuracy against sentence length, one curve per
  ## condition. Dots are item means; the curve is a binomial GAM over the
  ## trial-level 0/1 data, with a deliberately small shrinkage-penalised basis
  ## (k = 4, bs = "cs") so it describes the trend rather than tracking items --
  ## the shuffled condition sits near ceiling and an unconstrained smooth there
  ## is all wiggle and no signal.
  p_detection <- ggplot(mapping = aes(x = n_words, color = condition_label,
                                      fill = condition_label)) +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey70") +
    geom_point(data = item_means(trials, detected, item_version, n_words,
                                 condition_label),
               aes(y = value), alpha = 0.35) +
    geom_smooth(data = trials %>% mutate(correct = as.integer(detected)),
                aes(y = correct), method = "gam",
                method.args = list(family = "binomial"),
                formula = y ~ s(x, bs = "cs", k = 4)) +
    scale_color_manual(name = "Condition", values = error_colors,
                       aesthetics = c("colour", "fill")) +
    scale_x_continuous(name = "Sentence length (tokens)") +
    scale_y_continuous(name = "Probability of correct judgment", labels = scales::percent,
                       breaks = seq(0, 1, 0.25)) +
    coord_cartesian(ylim = c(0, 1)) +
    ggtitle(spec$label)
  print(p_detection)
  ggsave(sprintf("figures/error_awareness_detection_by_length_%s.png", spec$slug),
         p_detection, width = 8, height = 5.5)

  ## Model: does noticing an error depend on how much was displaced and on how
  ## long the sentence is? Same structure as the main-analysis acceptability
  ## model, in binomial form.
  detection_model <- fit_maximal(
    list(
      paste("detected ~ dl_distance_z * n_words_z +",
            "(dl_distance_z * n_words_z | participant_id) +", item_slope(c("dl_distance_z", "n_words_z"))),
      paste("detected ~ dl_distance_z * n_words_z +",
            "(dl_distance_z + n_words_z | participant_id) +", item_slope("dl_distance_z")),
      paste("detected ~ dl_distance_z * n_words_z +",
            "(dl_distance_z + n_words_z || participant_id) + (1 | item)")),
    function(f) glmer(f, family = binomial, control = glmer_ctrl,
                      data = model_trials %>%
                        mutate(dl_distance_z = zscore(dl_distance),
                               n_words_z = zscore(n_words))),
    paste("detection", spec$slug))
  cat("\n---- detection ~ DL distance x length --", spec$label, "----\n")
  print(round(summary(detection_model)$coefficients, 4))


  ## ---- 4.2 Are acceptability ratings independent of error awareness? ----
  cat("\n===== 4.2 Ratings and error awareness --", spec$label, "=====\n")

  ## Shuffled is dropped from this figure only: it is detected on nearly every
  ## trial, so the "judged incorrectly" panel holds too few of them to fit.
  p_awareness <- trials %>%
    filter(cond_family != "random") %>%
    mutate(condition_label = droplevels(condition_label)) %>%
    item_means(acceptability, item_version, n_words, condition_label, detection_label) %>%
    ggplot(aes(n_words, value, color = condition_label, fill = condition_label)) +
    geom_point(alpha = 0.35) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "tp")) +
    facet_wrap(~ detection_label, nrow = 1) +
    scale_color_manual(name = "Condition", values = error_colors,
                       aesthetics = c("colour", "fill")) +
    scale_x_continuous(name = "Sentence length (tokens)") +
    scale_y_continuous(name = "Acceptability", breaks = seq(1, 7, 2)) +
    coord_cartesian(ylim = c(1, 7)) +
    ggtitle(spec$label)
  print(p_awareness)
  ggsave(sprintf("figures/error_awareness_acceptability_by_detection_%s.png", spec$slug),
         p_awareness, width = 11, height = 5.5)

  ## Model: acceptability from the communication predictors (DL distance per
  ## token, plus length) crossed with whether the error judgment was correct.
  ## Correctness is sum-coded so the predictor main effects average over it.
  awareness_model <- fit_maximal(
    list(
      paste("acceptability ~ (weighted_dl_distance_z + n_words_z) * detected_c +",
            "(weighted_dl_distance_z + n_words_z + detected_c | participant_id) +",
            item_slope(c("weighted_dl_distance_z", "n_words_z", "detected_c"))),
      paste("acceptability ~ (weighted_dl_distance_z + n_words_z) * detected_c +",
            "(weighted_dl_distance_z + n_words_z + detected_c | participant_id) +",
            item_slope("weighted_dl_distance_z")),
      paste("acceptability ~ (weighted_dl_distance_z + n_words_z) * detected_c +",
            "(weighted_dl_distance_z + n_words_z + detected_c || participant_id) +",
            "(1 | item)")),
    function(f) lmer(f, control = lmer_ctrl,
                     data = model_trials %>%
                       mutate(weighted_dl_distance_z = zscore(weighted_dl_distance),
                              n_words_z = zscore(n_words),
                              detected_c = ifelse(detected, 0.5, -0.5))),
    paste("awareness", spec$slug))
  cat("\n---- acceptability ~ (DL/length + length) x judged correctly --",
      spec$label, "----\n")
  print(round(summary(awareness_model)$coefficients, 4))


  ## ---- 4.3 Correctly-judged trials vs. the comprehension-question version ----
  cat("\n===== 4.3 Versions --", spec$label, "=====\n")

  ## Two panels: the same materials run with meaning-comprehension questions
  ## (all trials, since there is no error judgment to condition on) and run with
  ## the error question (restricted to the trials judged correctly).
  ## Three panels: the same materials with meaning-comprehension questions, then
  ## with the error question over ALL trials, then restricted to the trials whose
  ## error judgment was correct. The middle panel is the one that keeps the
  ## comparison honest -- it shows the error-detection task's effect on ratings
  ## without the selection that restricting to caught errors introduces.
  ea_version_panel <- function(d, label) {
    d %>% transmute(item_version, cond_family, n_words, dl_distance, acceptability,
                    condition_label, panel = label)
  }
  version_panels <- bind_rows(
    reference_trials(spec$reference) %>%
      mutate(panel = "Comprehension questions"),
    ea_version_panel(trials, "Error detection\n(all trials)"),
    ea_version_panel(trials %>% filter(detected), "Error detection\n(judged correctly)")) %>%
    mutate(panel = factor(panel, levels = c("Comprehension questions",
                                            "Error detection\n(all trials)",
                                            "Error detection\n(judged correctly)")))

  p_version <- ggplot(item_means(version_panels, acceptability, panel,
                                 item_version, n_words, condition_label),
                      aes(n_words, value, color = condition_label,
                          fill = condition_label)) +
    geom_point(alpha = 0.35) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "tp")) +
    facet_wrap(~ panel, nrow = 1) +
    scale_color_manual(name = "Condition", values = error_colors,
                       aesthetics = c("colour", "fill")) +
    scale_x_continuous(name = "Sentence length (tokens)") +
    scale_y_continuous(name = "Acceptability", breaks = seq(1, 7, 2)) +
    coord_cartesian(ylim = c(1, 7)) +
    ggtitle(spec$label)
  print(p_version)
  ggsave(sprintf("figures/error_awareness_acceptability_by_version_%s.png", spec$slug),
         p_version, width = 15, height = 5.5)

  ## The communication account's second prediction, isolated: among sentences
  ## correctly judged error-free, longer should still be less acceptable.
  cat("\n---- error-free (correctly judged) acceptability ~ length --",
      spec$label, "----\n")
  ef_item <- if (spec$slopes_by_item) "(n_words_c | item)" else "(1 | item)"
  print(round(summary(lmer(
    as.formula(paste("acceptability ~ n_words_c + (n_words_c | participant_id) +", ef_item)),
    data = model_trials %>%
      filter(cond_family == "original", detected) %>% droplevels() %>%
      mutate(n_words_c = n_words - mean(n_words)),
    control = lmer_ctrl))$coefficients, 4))


  ## ---- 4.4/4.5 The paper's three models ----
  cat("\n===== 4.4/4.5 Model comparison --", spec$label, "=====\n")

  ## The three competing accounts in both forms the paper uses: shape-constrained
  ## GAMs over by-item means, and parametric linear mixed models over trials.
  ## Single-language versions, so the language random effects drop out. The lmers
  ## use REML = FALSE because the three models differ in their fixed effects and
  ## REML deviances are not comparable across them.
  scam_item_re <- function(x) if (spec$slopes_by_item)
    sprintf("s(item, bs='re') + s(item, %s, bs='re')", x) else "s(item, bs='re')"

  fit_ea_model_set <- function(d) {
    d_items <- d %>%
      droplevels() %>%
      group_by(item, item_version, cond_family, n_words, dl_distance,
               weighted_dl_distance, grammatical) %>%
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
        binary = scam(as.formula(paste("response_num ~ grammatical_numeric +",
                                       scam_item_re("grammatical_numeric"))), data = d_items),
        gradient = scam(as.formula(paste("response_num ~ s(dl_distance, bs='mpd') +",
                                         scam_item_re("dl_distance"))), data = d_items),
        comm = scam(as.formula(paste("response_num ~ s(weighted_dl_distance, bs='mpd') +",
                                     "s(n_words, bs='mpd') +",
                                     scam_item_re("weighted_dl_distance"))), data = d_items)),
      lmer = list(
        binary = lmer(as.formula(paste("acceptability ~ grammatical +",
                                       "(grammatical | participant_id) + (grammatical | item)")),
                      data = d_trials, REML = FALSE, control = lmer_ctrl),
        gradient = lmer(as.formula(paste("acceptability ~ dl_distance_z +",
                                         "(dl_distance_z | participant_id) + (dl_distance_z | item)")),
                        data = d_trials, REML = FALSE, control = lmer_ctrl),
        comm = lmer(as.formula(paste("acceptability ~ weighted_dl_distance_z + n_words_z +",
                                     "(weighted_dl_distance_z + n_words_z | participant_id) +",
                                     item_slope(c("weighted_dl_distance_z", "n_words_z")))),
                    data = d_trials, REML = FALSE, control = lmer_ctrl)))
  }

  ## Three nested views of the same trials, plus the comprehension-question
  ## experiment as the reference. The "All trials" and reference rows are not
  ## part of the questions as posed, but without them 4.4/4.5 cannot be read:
  ## restricting to correctly judged trials selects on a variable that is itself
  ## correlated with the rating.
  reference_set <- crit %>%
    filter(language == spec$reference, dl_distance <= 3) %>%
    mutate(item_id = str_remove(as.character(item), "^en(_experimental)?_")) %>%
    transmute(participant_id = factor(participant_id),
              item = factor(str_remove(item_id, "[sml]$")),
              item_version = item_id,
              cond_family = condition, n_words, dl_distance, weighted_dl_distance,
              grammatical, acceptability = response_num)

  ea_datasets <- list(
    "Comprehension questions" = reference_set,
    "All trials"              = model_trials,
    "Detected only"           = model_trials %>% filter(detected),
    "+ flagged originals"     = model_trials %>% filter(detected | cond_family == "original"))

  cache <- sprintf("models/ea_section4_models_%s.RData", spec$slug)
  if (!file.exists(cache)) {
    ea_model_sets <- lapply(names(ea_datasets), function(nm) {
      cat("Fitting model set for", nm, "\n"); fit_ea_model_set(ea_datasets[[nm]])
    })
    names(ea_model_sets) <- names(ea_datasets)
    save(ea_model_sets, file = cache)
  } else {
    load(cache)
  }

  tibble(dataset = names(ea_datasets),
         trials = vapply(ea_datasets, nrow, integer(1)),
         items = vapply(ea_datasets, function(d) n_distinct(d$item), integer(1)),
         originals = vapply(ea_datasets,
                            function(d) sum(d$cond_family == "original"), integer(1))) %>%
    print(n = Inf)

  for (nm in names(ea_model_sets)) {
    cat("\n---- convergence --", nm, "----\n")
    cat("  scam converged:",
        paste(names(ea_model_sets[[nm]]$scam),
              vapply(ea_model_sets[[nm]]$scam, function(m) m$conv, logical(1)),
              collapse = "  "), "\n")
    for (mn in names(ea_model_sets[[nm]]$lmer))
      cat("  lmer", mn, "- converged:", converged(ea_model_sets[[nm]]$lmer[[mn]]), "\n")
  }

  ## scam's own logLik() is conditional on the fitted random effects, so
  ## marginal_logLik() is used there; lmer's logLik() is already marginal. The
  ## two families are fit to different data, so likelihoods are only ever
  ## compared WITHIN a family and WITHIN a row.
  lik <- bind_rows(lapply(names(ea_model_sets), function(nm) bind_rows(
    tibble(dataset = nm, family = "Shape-constrained GAM\n(item means)",
           model = names(ea_model_sets[[nm]]$scam),
           log_likelihood = vapply(ea_model_sets[[nm]]$scam, marginal_logLik, numeric(1))),
    tibble(dataset = nm, family = "Linear mixed model\n(trials)",
           model = names(ea_model_sets[[nm]]$lmer),
           log_likelihood = vapply(ea_model_sets[[nm]]$lmer,
                                   function(m) as.numeric(logLik(m)), numeric(1)))))) %>%
    mutate(experiment = spec$label,
           dataset = factor(dataset, levels = names(ea_datasets)),
           model = factor(model, levels = c("binary", "gradient", "comm"),
                          labels = c("Binary grammaticality", "Gradient grammaticality",
                                     "Communication-based")),
           family = factor(family, levels = c("Linear mixed model\n(trials)",
                                              "Shape-constrained GAM\n(item means)"))) %>%
    group_by(dataset, family) %>%
    mutate(likelihood_diff = log_likelihood - max(log_likelihood)) %>%
    ungroup()

  lik %>% select(family, dataset, model, log_likelihood, likelihood_diff) %>%
    arrange(family, dataset, model) %>% print(n = Inf)
  lik %>% group_by(family, dataset) %>% slice_max(log_likelihood, n = 1) %>%
    ungroup() %>% select(family, dataset, best_model = model, log_likelihood) %>%
    print(n = Inf)

  ea_likelihoods_all <- bind_rows(ea_likelihoods_all, lik)
}


## ---- Detection figure, both experiments side by side ----
## The same by-item detection accuracy as the per-experiment figures, but with
## the two experiments in one frame so the length effect can be compared. Length
## is confounded with item in the English experiment and varies within item in
## the controlled one, which is the comparison this figure exists to show.
ea_trials_all <- ea_trials_all %>%
  mutate(experiment_label = factor(experiment_label,
                                   levels = c("English", "English (controlled)")))

p_detection_both <- ggplot(mapping = aes(x = n_words, color = condition_label,
                                         fill = condition_label)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey70") +
  geom_point(data = item_means(ea_trials_all, detected, experiment_label,
                               item_version, n_words, condition_label),
             aes(y = value), alpha = 0.3) +
  geom_smooth(data = ea_trials_all %>% mutate(correct = as.integer(detected)),
              aes(y = correct), method = "gam",
              method.args = list(family = "binomial"),
              formula = y ~ s(x, bs = "cs", k = 4)) +
  facet_wrap(~ experiment_label, nrow = 1) +
  scale_color_manual(name = "Condition", values = error_colors,
                     aesthetics = c("colour", "fill")) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Probability of correct judgment",
                     labels = scales::percent, breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1))
print(p_detection_both)
ggsave("figures/error_awareness_detection_by_experiment.png", p_detection_both,
       width = 11, height = 5.5)


## ---- Likelihood-difference figure, both experiments ----
## Shape-constrained GAMs only. The linear mixed models are still fitted and
## reported above, but they are left off this figure: they are fit to trials
## rather than item means, so their likelihoods are not on the same footing, and
## the paper's comparison is the shape-constrained one. The "+ flagged
## originals" robustness dataset is likewise reported above but not plotted.
ea_row_labels <- c("Comprehension questions" = "comprehension",
                   "All trials"              = "error awareness\n(all responses)",
                   "Detected only"           = "error awareness\n(judged correctly)")

ea_likelihoods_plot <- ea_likelihoods_all %>%
  filter(family == "Shape-constrained GAM\n(item means)",
         dataset %in% names(ea_row_labels)) %>%
  mutate(row = paste0(experiment, "\n", ea_row_labels[as.character(dataset)]),
         row = factor(row, levels = unlist(lapply(
           c("English", "English (controlled)"),
           function(e) paste0(e, "\n", ea_row_labels)))))

p_ea_likelihoods <- ggplot(ea_likelihoods_plot) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
  geom_point(aes(y = row, x = likelihood_diff, color = model, shape = model), size = 3) +
  scale_color_manual(name = "", values = c("#EB7AA3", "#1E88E5", "#004D40")) +
  scale_shape_manual(name = "", values = c(16, 17, 15)) +
  scale_x_continuous(name = "Log-likelihood difference from the best model") +
  scale_y_discrete(name = "", limits = rev) +
  ## Legend across the top, as in the other likelihood figures in the SI. The
  ## model names are kept on one line each -- the two-line form made three
  ## entries read like six.
  guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size = 15),
        axis.text.y = element_text(size = 13, lineheight = 0.9))
print(p_ea_likelihoods)
ggsave("figures/error_awareness_likelihood_differences.png", p_ea_likelihoods,
       width = 11, height = 6)
