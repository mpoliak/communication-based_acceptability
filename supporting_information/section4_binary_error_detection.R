## =============================================================================
## SI Section 4: Binary error detection (English)
##
## The English experiments in which participants also judged whether each
## sentence contained an error (../data/en_binary_error_detection,
## ../data/en_experimental_binary_error_detection), and the English exchange
## identification task (../data/en_exchange_identification). Figures:
## binary_error_detection_*.png, en_exchange_identification_*.png and
## en_swap_penalty_by_length.png.
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
## shuffled trials). Produces binary_error_detection_accuracy_by_length_*.png,
## binary_error_detection_acceptability_by_detection_*.png and
## binary_error_detection_acceptability_by_version_*.png (one per experiment,
## en and en_exp), binary_error_detection_accuracy_by_experiment.png and
## binary_error_detection_likelihood_differences.png.

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


## ---- The two binary error detection experiments ----
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

  read_csv("../data/en_binary_error_detection/data_en_binary_error_detection.csv",
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

  read_csv("../data/en_experimental_binary_error_detection/data_en_exp_binary_error_detection.csv",
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
  cat("## Binary error detection --", spec$label, "\n")
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
  ggsave(sprintf("figures/binary_error_detection_accuracy_by_length_%s.png", spec$slug),
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


  ## ---- 4.2 Are acceptability ratings independent of binary error detection? ----
  cat("\n===== 4.2 Ratings and binary error detection --", spec$label, "=====\n")

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
  ggsave(sprintf("figures/binary_error_detection_acceptability_by_detection_%s.png", spec$slug),
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
  ggsave(sprintf("figures/binary_error_detection_acceptability_by_version_%s.png", spec$slug),
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
ggsave("figures/binary_error_detection_accuracy_by_experiment.png", p_detection_both,
       width = 11, height = 5.5)


## ---- Likelihood-difference figure, both experiments ----
## Shape-constrained GAMs only. The linear mixed models are still fitted and
## reported above, but they are left off this figure: they are fit to trials
## rather than item means, so their likelihoods are not on the same footing, and
## the paper's comparison is the shape-constrained one. The "+ flagged
## originals" robustness dataset is likewise reported above but not plotted.
ea_row_labels <- c("Comprehension questions" = "comprehension",
                   "All trials"              = "binary error detection\n(all responses)",
                   "Detected only"           = "binary error detection\n(judged correctly)")

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
ggsave("figures/binary_error_detection_likelihood_differences.png", p_ea_likelihoods,
       width = 11, height = 6)


## ---- 4.7 Word identification (en_exchange_identification) ----
## A different task on the ORIGINAL English materials: every word is clickable
## and the participant clicks the ones they think are misplaced. No rating, no
## Yes/No question. Only `original` and `swap_1`, so exactly two words are
## displaced when anything is (verified below). Surveyor's
## per_word_multiple_choice always writes one row per trial -- `NONE` when
## nothing was clicked -- so a missing row means the trial was never reached,
## which makes the completeness check below meaningful.
cat("\n\n#############################################################\n")
cat("## Word identification (en_exchange_identification)\n")
cat("#############################################################\n")

## Degree of agreement about WHERE the error is: mean pairwise Jaccard overlap
## between the sets of positions different participants clicked. 1 = everyone
## picked exactly the same words, 0 = no overlap at all. NA when fewer than two
## participants clicked anything, since agreement is undefined for one set.
## For swap_1, marking one of the two displaced words and marking both are the
## same event -- either way the participant has pointed at the error. Collapse
## the two target positions into a single unit (-1L, which cannot collide with a
## real 1-based position) before measuring overlap, so those two responses count
## as agreeing rather than as half-agreeing. Does not apply to the original
## condition, where no word was displaced and there is nothing to collapse.
## A response counts as correctly located if it touches at least one displaced
## word AND everything it marks sits in one unbroken run containing them --
## marking the two swapped words, just one of them, or either plus the words
## immediately beside them all point at the same disfluent region. Verified:
## every swap_1 target is exactly two ADJACENT positions, so the union is
## contiguous exactly when the participant marked no distant extras.
marked_target_region <- function(clicked, target) {
  if (!length(intersect(clicked, target))) return(FALSE)
  span <- sort(union(clicked, target))
  all(diff(span) == 1)
}

collapse_targets <- function(clicked, target) {
  c(setdiff(clicked, target),
    if (length(intersect(clicked, target))) -1L else integer(0))
}

agreement_of <- function(click_sets) {
  click_sets <- click_sets[lengths(click_sets) > 0]
  if (length(click_sets) < 2) return(NA_real_)
  pairs <- combn(length(click_sets), 2)
  mean(apply(pairs, 2, function(ij) {
    a <- click_sets[[ij[1]]]; b <- click_sets[[ij[2]]]
    length(intersect(a, b)) / length(union(a, b))
  }))
}

id_materials <- read_csv("../data/en_exchange_identification/materials_en_exchange_identification.csv",
                         show_col_types = FALSE) %>%
  transmute(item = as.character(item), condition, sentence)

## The displaced positions for each item, recovered by diffing swap_1 against
## its own original. literal_response is NA in the stimuli file, so this is the
## ground truth for scoring. 1-based, to match the stored response format.
id_targets <- id_materials %>%
  pivot_wider(names_from = condition, values_from = sentence) %>%
  rowwise() %>%
  mutate(target = list(which(str_split(original, " ")[[1]] !=
                             str_split(swap_1, " ")[[1]]))) %>%
  ungroup() %>%
  select(item, target)

id_raw <- read_csv("../data/en_exchange_identification/en_exchange_identification_data.csv",
                   show_col_types = FALSE) %>%
  mutate(item = as.character(item))

## ---- did it save properly? ----
cat("\n---- export integrity ----\n")
cat("participants in file:", n_distinct(id_raw$participant_id), "\n")
id_raw %>%
  count(participant_id, name = "trials") %>%
  count(trials, name = "n_participants") %>%
  print()
cat("(72 trials each = complete; per_word_multiple_choice always emits a row)\n")

## Every position:word pair should name the word actually at that position in
## that participant's version of the sentence. Words are percent-escaped for
## % , and : so they survive the encoding, hence URLdecode.
id_clicks <- id_raw %>%
  filter(response != "NONE") %>%
  separate_rows(response, sep = ",") %>%
  separate(response, into = c("position", "word"), sep = ":", extra = "merge") %>%
  mutate(position = as.integer(position),
         word = map_chr(word, URLdecode)) %>%
  left_join(id_materials, by = c("item", "cond" = "condition")) %>%
  mutate(word_at_position = map2_chr(sentence, position,
                                     ~ str_split(.x, " ")[[1]][.y]))
cat("position:word pairs:", nrow(id_clicks),
    "| matching the stimulus:", sum(id_clicks$word == id_clicks$word_at_position),
    "| out of range:", sum(is.na(id_clicks$word_at_position)), "\n")

## `test` is the experimenter's own run-through, not a participant.
id_trials <- id_raw %>%
  filter(participant_id != "test") %>%
  left_join(id_targets, by = "item") %>%
  mutate(clicked = map(response, function(resp) {
           if (resp == "NONE") return(integer(0))
           as.integer(str_extract(str_split(resp, ",")[[1]], "^[0-9]+"))
         }),
         n_clicked = lengths(clicked),
         flagged = n_clicked > 0,
         ## correct location = at least one of the two displaced words clicked
         hit = map2_lgl(clicked, target, ~ length(intersect(.x, .y)) > 0),
         both = map2_lgl(clicked, target, ~ length(intersect(.x, .y)) == 2))

cat("\nreal participants:", n_distinct(id_trials$participant_id),
    "| trials:", nrow(id_trials), "\n")
id_trials %>% count(cond, flagged) %>% print()

## ---- (1) original sentences: false alarms, and where they land ----
cat("\n---- (1) ORIGINAL sentences: is anything flagged, and do people agree where? ----\n")
id_orig <- id_trials %>% filter(cond == "original")
cat("false-alarm rate (clicked at least one word in a sentence with no error):",
    round(mean(id_orig$flagged), 3), "\n")
cat("mean words clicked when a false alarm happened:",
    round(mean(id_orig$n_clicked[id_orig$flagged]), 2), "\n")

id_orig_by_item <- id_orig %>%
  group_by(item) %>%
  summarize(n = n(), n_flagged = sum(flagged), flag_rate = mean(flagged),
            agreement_where = agreement_of(clicked[flagged]), .groups = "drop")
cat("\nitems by how often they drew a false alarm:\n")
id_orig_by_item %>% count(n_flagged, name = "items") %>% print()
cat("\nthe items most often flagged despite having no error:\n")
id_orig_by_item %>% arrange(desc(flag_rate), desc(n_flagged)) %>%
  filter(n_flagged > 0) %>% head(10) %>% print()

## ---- (2) swap_1 sentences: noticed, located, agreed ----
cat("\n---- (2) SWAP_1 sentences: noticed, located correctly, agreement ----\n")
id_swap <- id_trials %>% filter(cond == "swap_1")
id_swap %>%
  summarize(trials = n(),
            noticed = mean(flagged),
            located_correctly = mean(hit),
            located_given_noticed = mean(hit[flagged]),
            both_words = mean(both),
            both_given_noticed = mean(both[flagged]),
            mean_words_clicked = mean(n_clicked[flagged])) %>%
  print(width = Inf)

id_swap_by_item <- id_swap %>%
  mutate(clicked_collapsed = map2(clicked, target, collapse_targets)) %>%
  group_by(item) %>%
  summarize(n = n(), n_noticed = sum(flagged), noticed = mean(flagged),
            n_hit = sum(hit), located = mean(hit),
            agreement_where = agreement_of(clicked_collapsed[flagged]),
            agreement_raw = agreement_of(clicked[flagged]), .groups = "drop")
cat("\nper-item distribution of correct localisation:\n")
id_swap_by_item %>% count(n_hit, n, name = "items") %>% print(n = Inf)
cat("\nitems where the displaced words were never located:\n")
id_swap_by_item %>% filter(n_hit == 0) %>%
  left_join(id_materials %>% filter(condition == "swap_1") %>% select(item, sentence),
            by = "item") %>%
  mutate(sentence = str_trunc(sentence, 70)) %>%
  select(item, n, n_noticed, sentence) %>% print(n = Inf)

cat("\n---- agreement summary (mean pairwise Jaccard of clicked positions) ----\n")
cat("swap_1 uses the collapsed measure: hitting either or both displaced words\n")
cat("counts as the same event. Raw (each word separate) shown for comparison.\n")
id_swap_by_item %>%
  summarize(collapsed = mean(agreement_where, na.rm = TRUE),
            raw = mean(agreement_raw, na.rm = TRUE)) %>%
  print()
bind_rows(
  id_orig_by_item %>% transmute(cond = "original", agreement_where),
  id_swap_by_item %>% transmute(cond = "swap_1", agreement_where)) %>%
  filter(!is.na(agreement_where)) %>%
  group_by(cond) %>%
  summarize(items_with_2plus_flaggers = n(),
            mean_agreement = mean(agreement_where), .groups = "drop") %>%
  print()

## ---- 4.7 figures ----
## Sentence length comes from `crit` (the same English experiment these
## materials are drawn from), so lengths match the rest of section 4.
id_lengths <- crit %>%
  filter(language == "en") %>%
  distinct(item, n_words) %>%
  transmute(item = str_remove(as.character(item), "^en_"), n_words) %>%
  distinct(item, .keep_all = TRUE)

id_plot_trials <- id_trials %>%
  left_join(id_lengths, by = "item") %>%
  mutate(condition_label = condition_labels(cond))

## Figure 1: does noticing the swap depend on sentence length, and does
## localisation track it? Points are item means, curves are binomial GAMs over
## the trial-level 0/1 data, same idiom as the detection figures above. The
## false-alarm line is the original condition, where any click is an error.
id_by_length <- bind_rows(
  id_plot_trials %>% filter(cond == "swap_1") %>%
    transmute(item, n_words, measure = "Noticed the swap", y = as.integer(flagged)),
  id_plot_trials %>% filter(cond == "swap_1") %>%
    transmute(item, n_words, measure = "Located it correctly", y = as.integer(hit)),
  id_plot_trials %>% filter(cond == "original") %>%
    transmute(item, n_words, measure = "False alarm (no error present)",
              y = as.integer(flagged))) %>%
  mutate(measure = factor(measure, levels = c("Noticed the swap",
                                              "Located it correctly",
                                              "False alarm (no error present)")))

id_measure_colors <- c("Noticed the swap" = "#1E88E5",
                       "Located it correctly" = "#004D40",
                       "False alarm (no error present)" = "#EB7AA3")

p_id_length <- ggplot(mapping = aes(x = n_words, color = measure, fill = measure)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey80") +
  geom_point(data = id_by_length %>%
               group_by(measure, item, n_words) %>%
               summarize(value = mean(y), .groups = "drop"),
             aes(y = value), alpha = 0.3) +
  geom_smooth(data = id_by_length, aes(y = y), method = "gam",
              method.args = list(family = "binomial"),
              formula = y ~ s(x, bs = "cs", k = 4)) +
  scale_color_manual(name = "", values = id_measure_colors,
                     aesthetics = c("colour", "fill")) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Proportion of trials", labels = scales::percent,
                     breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1)) +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
  theme(legend.position = "top", legend.text = element_text(size = 13))
print(p_id_length)
ggsave("figures/en_exchange_identification_by_length.png", p_id_length, width = 9, height = 5.5)

## Figure 2: what actually got clicked on a swap_1 trial. "Correct" means every
## word clicked was one of the two displaced ones -- clicking just one of the
## pair still points at the right place, so it counts. "Mixed" is the partly
## right case: at least one displaced word, but irrelevant words too.
id_outcomes <- id_plot_trials %>%
  filter(cond == "swap_1") %>%
  mutate(n_extra = map2_int(clicked, target, ~ length(setdiff(.x, .y))),
         n_on_target = map2_int(clicked, target, ~ length(intersect(.x, .y))),
         in_region = map2_lgl(clicked, target, marked_target_region),
         outcome = case_when(
           !flagged          ~ "Wrong\n(clicked nothing)",
           in_region         ~ "Correct\n(target \u00b1 adjacent words)",
           n_on_target == 0  ~ "Wrong\n(other words only)",
           TRUE              ~ "Mixed\n(target + distant words)"),
         outcome = factor(outcome, levels = c(
           "Correct\n(target \u00b1 adjacent words)", "Mixed\n(target + distant words)",
           "Wrong\n(other words only)", "Wrong\n(clicked nothing)")))

## setNames, not c("..." = ...): a name written as an argument tag is converted
## to the native encoding, which in a non-UTF-8 locale turns the "\u00b1" into
## the literal text "<U+00B1>". The palette then no longer matches the data and
## the Correct level silently falls back to grey.
id_outcome_colors <- setNames(
  c("#004D40", "#FFC107", "#EB7AA3", "#BBBBBB"),
  c("Correct\n(target \u00b1 adjacent words)", "Mixed\n(target + distant words)",
    "Wrong\n(other words only)", "Wrong\n(clicked nothing)"))

p_id_outcome <- id_outcomes %>%
  count(outcome) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(outcome, prop, fill = outcome)) +
  geom_col(alpha = 0.9, show.legend = FALSE) +
  scale_fill_manual(values = id_outcome_colors) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", 100 * prop, n)), vjust = -0.3, size = 4.5) +
  scale_y_continuous(name = "Proportion of swap_1 trials", labels = scales::percent,
                     expand = expansion(mult = c(0, 0.18))) +
  scale_x_discrete(name = "") +
  theme(axis.text.x = element_text(size = 12))
print(p_id_outcome)
ggsave("figures/en_exchange_identification_outcomes.png", p_id_outcome, width = 10, height = 5.5)

## Figure 2b: the same four outcomes as a function of sentence length. Each
## outcome is its own 0/1 indicator, so the four curves sum to 1 at any length.
## Points are item means, curves are binomial GAMs over the trial-level data --
## same idiom as the other by-length figures in this section.
id_outcome_long <- id_outcomes %>%
  select(item, n_words, outcome) %>%
  tidyr::expand_grid(level = levels(id_outcomes$outcome)) %>%
  mutate(y = as.integer(as.character(outcome) == level),
         level = factor(level, levels = levels(id_outcomes$outcome)))

p_id_outcome_length <- ggplot(mapping = aes(x = n_words, color = level, fill = level)) +
  geom_point(data = id_outcome_long %>%
               group_by(level, item, n_words) %>%
               summarize(value = mean(y), .groups = "drop"),
             aes(y = value), alpha = 0.25) +
  geom_smooth(data = id_outcome_long, aes(y = y), method = "gam",
              method.args = list(family = "binomial"),
              formula = y ~ s(x, bs = "cs", k = 4)) +
  scale_color_manual(name = "", values = id_outcome_colors,
                     aesthetics = c("colour", "fill")) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Proportion of swap_1 trials", labels = scales::percent,
                     breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1)) +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
  theme(legend.position = "top", legend.text = element_text(size = 11))
print(p_id_outcome_length)
ggsave("figures/en_exchange_identification_outcomes_by_length.png", p_id_outcome_length,
       width = 11, height = 5.5)

## Figures 2c/2d: the same treatment for the ORIGINAL sentences. There are no
## displaced words, so the only correct response is to click nothing and any
## click at all is a false alarm -- two outcomes, not graded by how many words
## were marked.
id_orig_outcomes <- id_plot_trials %>%
  filter(cond == "original") %>%
  mutate(outcome = factor(ifelse(n_clicked == 0,
                                 "Correct\n(clicked nothing)",
                                 "False alarm\n(1+ words marked)"),
                          levels = c("Correct\n(clicked nothing)",
                                     "False alarm\n(1+ words marked)")))

id_orig_colors <- c("Correct\n(clicked nothing)"      = "#004D40",
                    "False alarm\n(1+ words marked)"  = "#EB7AA3")

p_id_orig_outcome <- id_orig_outcomes %>%
  count(outcome) %>%
  mutate(prop = n / sum(n)) %>%
  ggplot(aes(outcome, prop, fill = outcome)) +
  geom_col(alpha = 0.9, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", 100 * prop, n)), vjust = -0.3, size = 4.5) +
  scale_fill_manual(values = id_orig_colors) +
  scale_y_continuous(name = "Proportion of original trials", labels = scales::percent,
                     expand = expansion(mult = c(0, 0.18))) +
  scale_x_discrete(name = "") +
  theme(axis.text.x = element_text(size = 12))
print(p_id_orig_outcome)
ggsave("figures/en_exchange_identification_outcomes_original.png", p_id_orig_outcome,
       width = 10, height = 5.5)

id_orig_outcome_long <- id_orig_outcomes %>%
  select(item, n_words, outcome) %>%
  tidyr::expand_grid(level = levels(id_orig_outcomes$outcome)) %>%
  mutate(y = as.integer(as.character(outcome) == level),
         level = factor(level, levels = levels(id_orig_outcomes$outcome)))

p_id_orig_length <- ggplot(mapping = aes(x = n_words, color = level, fill = level)) +
  geom_point(data = id_orig_outcome_long %>%
               group_by(level, item, n_words) %>%
               summarize(value = mean(y), .groups = "drop"),
             aes(y = value), alpha = 0.25) +
  geom_smooth(data = id_orig_outcome_long, aes(y = y), method = "gam",
              method.args = list(family = "binomial"),
              formula = y ~ s(x, bs = "cs", k = 4)) +
  scale_color_manual(name = "", values = id_orig_colors,
                     aesthetics = c("colour", "fill")) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Proportion of original trials", labels = scales::percent,
                     breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1)) +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
  theme(legend.position = "top", legend.text = element_text(size = 11))
print(p_id_orig_length)
ggsave("figures/en_exchange_identification_outcomes_original_by_length.png", p_id_orig_length,
       width = 11, height = 5.5)

cat("\n---- original-condition outcome breakdown ----\n")
id_orig_outcomes %>% count(outcome) %>% mutate(prop = round(n / sum(n), 3)) %>% print()

## Figure 3: how much do participants agree with each other about where the
## error is? One point per item, mean pairwise Jaccard of the clicked position
## sets among participants who clicked anything.
id_agreement <- bind_rows(
  id_orig_by_item %>% transmute(item, condition_label = condition_labels("original"),
                                agreement_where, n_flagged),
  id_swap_by_item %>% transmute(item, condition_label = condition_labels("swap_1"),
                                agreement_where, n_flagged = n_noticed)) %>%
  filter(!is.na(agreement_where))

p_id_agreement <- ggplot(id_agreement, aes(agreement_where, fill = condition_label)) +
  geom_histogram(binwidth = 0.1, boundary = 0, color = "white", alpha = 0.9) +
  facet_wrap(~ condition_label, nrow = 1, scales = "free_y") +
  scale_fill_manual(name = "Condition", values = error_colors) +
  ## coord_cartesian, not scale limits: a hard limit clips the bin containing
  ## perfect agreement (1.0) and silently drops those items from the figure.
  scale_x_continuous(name = "Agreement about which words are misplaced\n(mean pairwise Jaccard)",
                     breaks = seq(0, 1, 0.25)) +
  coord_cartesian(xlim = c(-0.05, 1.05)) +
  scale_y_continuous(name = "Items") +
  theme(legend.position = "none")
print(p_id_agreement)
ggsave("figures/en_exchange_identification_agreement.png", p_id_agreement, width = 10, height = 5)

## Figure 4: among swap_1 trials where the participant marked SOMETHING, how
## much do they agree about where? Plotted against sentence length. One point
## per item, sized by how many participants flagged it -- an item with only two
## flaggers gives a very noisy agreement estimate, so the smooth is weighted by
## that count. Uses the collapsed measure (one or both displaced words = same
## event).
##
## Raw Jaccard is not comparable across lengths: two people each clicking two
## words overlap by luck far more often in a 6-word sentence than a 30-word one,
## so part of the apparent decline is just a shrinking chance baseline. The
## second line corrects for it. Per item, chance is found by permutation --
## keep each participant's number of clicks, randomise WHICH words within that
## sentence's length, collapse the targets the same way, recompute. Corrected
## agreement is (observed - chance)/(1 - chance): 0 is chance, 1 is perfect,
## on the same scale at every length.
set.seed(2026)
id_chance_agreement <- function(click_sets, sentence_len, target_positions,
                                reps = 200) {
  sizes <- lengths(click_sets)
  sizes <- sizes[sizes > 0]
  if (length(sizes) < 2) return(NA_real_)
  mean(vapply(seq_len(reps), function(i) {
    perm <- lapply(sizes, function(k)
      collapse_targets(sample.int(sentence_len, k), target_positions))
    agreement_of(perm)
  }, numeric(1)), na.rm = TRUE)
}

id_agree_length <- id_swap %>%
  left_join(id_lengths, by = "item") %>%
  mutate(clicked_collapsed = map2(clicked, target, collapse_targets)) %>%
  group_by(item, n_words) %>%
  summarize(n_noticed = sum(flagged),
            observed = agreement_of(clicked_collapsed[flagged]),
            chance = id_chance_agreement(clicked[flagged], first(n_words),
                                         first(target)),
            .groups = "drop") %>%
  filter(!is.na(observed)) %>%
  mutate(corrected = (observed - chance) / (1 - chance))

## Only the observed agreement is plotted. The chance-corrected version is still
## computed and printed below -- it is the check that the decline with length is
## not just a shrinking chance baseline -- but it tracks the observed line so
## closely that drawing both only crowded the figure.
p_id_agree_length <- ggplot(id_agree_length, aes(n_words, observed)) +
  geom_point(aes(size = n_noticed), alpha = 0.45, color = "#1E88E5") +
  geom_smooth(aes(weight = n_noticed), method = "gam",
              formula = y ~ s(x, bs = "cs", k = 4),
              color = "#1E88E5", fill = "#1E88E5") +
  scale_size_continuous(name = "participants who\nmarked something",
                        range = c(1, 6)) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Agreement about which words are misplaced\n(mean pairwise Jaccard)",
                     breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1)) +
  ggtitle("swap_1, trials where an error was marked") +
  theme(legend.position = "right", plot.title = element_text(size = 15))
print(p_id_agree_length)
ggsave("figures/en_exchange_identification_agreement_by_length.png", p_id_agree_length,
       width = 9, height = 5.5)

cat("\n---- agreement vs length (swap_1, flagged trials only) ----\n")
id_agree_length %>%
  mutate(bin = cut(n_words, breaks = c(0, 10, 15, 20, 25, Inf),
                   labels = c("<=10", "11-15", "16-20", "21-25", "26+"))) %>%
  group_by(bin) %>%
  summarize(items = n(), mean_flaggers = mean(n_noticed),
            observed = mean(observed), chance = mean(chance),
            corrected = mean(corrected), .groups = "drop") %>%
  print()

cat("\n---- figure data check ----\n")
id_outcomes %>% count(outcome) %>% mutate(prop = round(n / sum(n), 3)) %>% print()
id_agreement %>% group_by(condition_label) %>%
  summarize(items = n(), mean_agreement = mean(agreement_where),
            median = median(agreement_where), .groups = "drop") %>% print()


## ---- 4.8 The original English experiment, restricted to validated items ----
## en_exchange_identification tells us, item by item, whether the swap is actually findable.
## Keep only items where more than THRESH of swap_1 responses located the error
## (target plus adjacent words) AND more than THRESH of original responses
## correctly clicked nothing, then refit the paper's three accounts on the
## ORIGINAL English acceptability experiment restricted to those items. Run at
## two thresholds so the conclusion is not an artefact of where the line is put.
##
## READ THE LENGTH STATISTICS BELOW BEFORE INTERPRETING ANY OF THIS. The filter
## is not independent of length: detection fails on long sentences, so screening
## on detection throws the long items out. The surviving set tops out well short
## of the full range, which is a problem for models whose whole claim is about
## length.
cat("\n\n#############################################################\n")
cat("## 4.8 Original English experiment, validated items only\n")
cat("#############################################################\n")

id_item_quality <- id_trials %>%
  mutate(in_region = map2_lgl(clicked, target, marked_target_region)) %>%
  group_by(item, cond) %>%
  summarize(rate = mean(if (first(cond) == "swap_1") in_region else !flagged),
            n = n(), .groups = "drop") %>%
  pivot_wider(names_from = cond, values_from = c(rate, n))

validated_analysis <- function(thresh) {
  tag <- sprintf("%d", round(thresh * 100))
  cat("\n\n=============================================================\n")
  cat(sprintf("== threshold: >%s%% hits AND >%s%% correct rejections\n", tag, tag))
  cat("=============================================================\n")

  good <- id_item_quality %>%
    filter(rate_swap_1 > thresh, rate_original > thresh) %>%
    pull(item)
  cat("\nitems passing both thresholds:", length(good), "of",
      n_distinct(id_item_quality$item), "\n")

  ## What do the surviving materials look like? This is the key diagnostic.
  cat("\n---- sentence length of kept vs dropped items ----\n")
  id_lengths %>%
    mutate(kept = item %in% good) %>%
    group_by(kept) %>%
    summarize(items = n(), min = min(n_words), max = max(n_words),
              mean = mean(n_words), median = median(n_words), .groups = "drop") %>%
    print()
  cat("\nlength distribution of kept items:\n")
  id_lengths %>%
    filter(item %in% good) %>%
    mutate(bin = cut(n_words, c(0, 10, 15, 20, 25, Inf),
                     labels = c("<=10", "11-15", "16-20", "21-25", "26+"))) %>%
    count(bin, .drop = FALSE) %>% print()

  ## Refit the three accounts on the original English acceptability data,
  ## restricted to these items. en varies length BETWEEN items, so no by-item
  ## length slope (same convention as the main en analysis above).
  dat <- crit %>%
    filter(language == "en", dl_distance <= 3) %>%
    mutate(item_num = str_remove(as.character(item), "^en_")) %>%
    filter(item_num %in% good) %>%
    transmute(participant_id = factor(participant_id), item = factor(item_num),
              item_version = item_num, cond_family = condition, n_words,
              dl_distance, weighted_dl_distance, grammatical,
              acceptability = response_num)

  cat("\n---- data going into the models ----\n")
  dat %>%
    summarize(trials = n(), participants = n_distinct(participant_id),
              items = n_distinct(item), length_min = min(n_words),
              length_max = max(n_words)) %>% print()
  dat %>% count(cond_family) %>% print()

  items_df <- dat %>%
    group_by(item, cond_family, n_words, dl_distance, weighted_dl_distance, grammatical) %>%
    summarize(response_num = mean(acceptability), .groups = "drop") %>%
    mutate(item = factor(item),
           grammatical_numeric = ifelse(grammatical == "grammatical", .5, -.5))
  trials_df <- dat %>%
    mutate(dl_distance_z = zscore(dl_distance),
           weighted_dl_distance_z = zscore(weighted_dl_distance),
           n_words_z = zscore(n_words))

  cache <- sprintf("models/ea_section4_models_validated_%s.RData", tag)
  if (!file.exists(cache)) {
    cat("\nFitting model set on validated items\n")
    models <- list(
      scam = list(
        binary = scam(response_num ~ grammatical_numeric + s(item, bs = "re"),
                      data = items_df),
        gradient = scam(response_num ~ s(dl_distance, bs = "mpd") + s(item, bs = "re"),
                        data = items_df),
        comm = scam(response_num ~ s(weighted_dl_distance, bs = "mpd") +
                      s(n_words, bs = "mpd") + s(item, bs = "re"), data = items_df)),
      lmer = list(
        binary = lmer(acceptability ~ grammatical + (grammatical | participant_id) +
                        (grammatical | item), data = trials_df,
                      REML = FALSE, control = lmer_ctrl),
        gradient = lmer(acceptability ~ dl_distance_z + (dl_distance_z | participant_id) +
                          (dl_distance_z | item), data = trials_df,
                        REML = FALSE, control = lmer_ctrl),
        comm = lmer(acceptability ~ weighted_dl_distance_z + n_words_z +
                      (weighted_dl_distance_z + n_words_z | participant_id) +
                      (weighted_dl_distance_z | item), data = trials_df,
                    REML = FALSE, control = lmer_ctrl)))
    save(models, file = cache)
  } else {
    load(cache)
  }

  cat("\n---- convergence ----\n")
  cat("  scam:", paste(names(models$scam),
      vapply(models$scam, function(m) m$conv, logical(1)), collapse = "  "), "\n")
  for (mn in names(models$lmer))
    cat("  lmer", mn, "- converged:", converged(models$lmer[[mn]]), "\n")

  lik <- bind_rows(
    tibble(family = "Shape-constrained GAM (item means)", model = names(models$scam),
           log_likelihood = vapply(models$scam, marginal_logLik, numeric(1))),
    tibble(family = "Linear mixed model (trials)", model = names(models$lmer),
           log_likelihood = vapply(models$lmer,
                                   function(m) as.numeric(logLik(m)), numeric(1)))) %>%
    group_by(family) %>%
    mutate(likelihood_diff = log_likelihood - max(log_likelihood)) %>%
    ungroup()
  cat("\n---- model comparison on validated items ----\n")
  print(lik, n = Inf)

  ## Figure: the by-item acceptability picture on this validated subset.
  p <- dat %>%
    mutate(condition_label = condition_labels(cond_family)) %>%
    item_means(acceptability, item_version, n_words, condition_label) %>%
    ggplot(aes(n_words, value, color = condition_label, fill = condition_label)) +
    geom_point(alpha = 0.4) +
    geom_smooth(method = "gam", formula = y ~ s(x, bs = "tp")) +
    scale_color_manual(name = "Condition", values = error_colors,
                       aesthetics = c("colour", "fill")) +
    scale_x_continuous(name = "Sentence length (tokens)") +
    scale_y_continuous(name = "Acceptability", breaks = seq(1, 7, 2)) +
    coord_cartesian(ylim = c(1, 7)) +
    ggtitle(sprintf("English, %d items with accuracy > %s%%", length(good), tag))
  print(p)
  ggsave(sprintf("figures/en_exchange_identification_validated_items_acceptability_%s.png", tag), p,
         width = 9, height = 5.5)

  list(items = good, lik = lik)
}

id_validated_runs <- lapply(c(0.80, 0.75), validated_analysis)
names(id_validated_runs) <- c("80", "75")
id_good_items <- id_validated_runs[["80"]]$items

cat("\n---- GAM comparison side by side, both thresholds ----\n")
imap_dfr(id_validated_runs, function(r, tag)
  r$lik %>% filter(family == "Shape-constrained GAM (item means)") %>%
    transmute(threshold = paste0(">", tag, "%"), model,
              log_likelihood = round(log_likelihood, 1),
              likelihood_diff = round(likelihood_diff, 1))) %>%
  print(n = Inf)


## ---- 4.9 Is the communication result just dislike of long sentences? ----
## The worry: the communication model is the only one carrying a length term, so
## maybe it wins simply because people rate long sentences lower, and the
## dl/length weighting does no real work.
##
## Two things discriminate the accounts:
##   (a) If it were a plain length preference, adding length as a separate
##       predictor to the gradient model should capture it. Fit that model, plus
##       a length-only baseline, and see.
##   (b) Length aversion is a MAIN EFFECT: every condition should decline with
##       length together. The communication account instead says the SAME error
##       costs less in a longer sentence -- an interaction, with swapped
##       sentences declining less or improving. The sign of the per-condition
##       length slope settles it.
cat("\n\n#############################################################\n")
cat("## 4.9 Length preference vs communication\n")
cat("#############################################################\n")

id_en_base <- crit %>%
  filter(language == "en", dl_distance <= 3) %>%
  mutate(item_num = str_remove(as.character(item), "^en_"))

id_thresholds <- list(
  "All 72 items" = unique(id_en_base$item_num),
  "Validated >75%" = id_item_quality %>% filter(rate_swap_1 > 0.75, rate_original > 0.75) %>% pull(item),
  "Validated >80%" = id_good_items)

## (a) five models, including a length-only baseline and gradient + length
id_five_models <- function(d) {
  di <- d %>%
    group_by(item_num, condition, n_words, dl_distance, weighted_dl_distance, grammatical) %>%
    summarize(response_num = mean(response_num), .groups = "drop") %>%
    mutate(item = factor(item_num),
           grammatical_numeric = ifelse(grammatical == "grammatical", .5, -.5))
  m <- list(
    "Communication\n(dl/length + length)" =
      scam(response_num ~ s(weighted_dl_distance, bs="mpd") + s(n_words, bs="mpd") + s(item, bs="re"), data = di),
    "Gradient + length" =
      scam(response_num ~ s(dl_distance, bs="mpd") + s(n_words, bs="mpd") + s(item, bs="re"), data = di),
    "Gradient" = scam(response_num ~ s(dl_distance, bs="mpd") + s(item, bs="re"), data = di),
    "Binary" = scam(response_num ~ grammatical_numeric + s(item, bs="re"), data = di),
    "Length only" = scam(response_num ~ s(n_words, bs="mpd") + s(item, bs="re"), data = di))
  tibble(model = names(m), ll = vapply(m, marginal_logLik, numeric(1)))
}

id_five_lik <- imap_dfr(id_thresholds, function(items, label) {
  id_five_models(id_en_base %>% filter(item_num %in% items)) %>% mutate(subset = label)
}) %>%
  group_by(subset) %>% mutate(diff = ll - max(ll)) %>% ungroup() %>%
  mutate(subset = factor(subset, levels = names(id_thresholds)),
         model = factor(model, levels = rev(c("Communication\n(dl/length + length)",
                        "Gradient + length", "Gradient", "Binary", "Length only"))))
cat("\n---- five-model comparison ----\n")
id_five_lik %>% arrange(subset, desc(ll)) %>%
  mutate(ll = round(ll, 1), diff = round(diff, 1)) %>% print(n = Inf)

p_id_five <- ggplot(id_five_lik, aes(diff, model, color = subset, shape = subset)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey80") +
  geom_point(size = 3, position = position_dodge(width = 0.4)) +
  scale_color_manual(name = "", values = c("#004D40", "#1E88E5", "#EB7AA3")) +
  scale_shape_manual(name = "", values = c(16, 17, 15)) +
  scale_x_continuous(name = "Log-likelihood difference from the best model") +
  scale_y_discrete(name = "") +
  guides(color = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  theme(legend.position = "top", axis.text.y = element_text(size = 12, lineheight = 0.9))
print(p_id_five)
ggsave("figures/en_exchange_identification_length_vs_communication.png", p_id_five, width = 10, height = 5.5)

## (b) the per-condition length slope -- the sign is the argument
id_slopes <- imap_dfr(id_thresholds, function(items, label) {
  id_en_base %>%
    filter(item_num %in% items, condition %in% c("original", "swap_1", "swap_3")) %>%
    group_by(condition) %>%
    group_modify(~ {
      m <- lmer(response_num ~ n_words + (1 | participant_id) + (1 | item),
                data = .x, control = lmer_ctrl)
      s <- summary(m)$coefficients
      tibble(slope = s[2, 1], se = s[2, 2], p = s[2, 5])
    }) %>% mutate(subset = label)
}) %>%
  mutate(subset = factor(subset, levels = names(id_thresholds)),
         condition_label = condition_labels(condition))
cat("\n---- acceptability change per extra token, within condition ----\n")
cat("length aversion predicts every slope < 0; communication predicts a sign flip\n")
id_slopes %>% mutate(across(c(slope, se), ~ round(.x, 4)), p = signif(p, 3)) %>%
  select(subset, condition_label, slope, se, p) %>% print(n = Inf)

p_id_slopes <- ggplot(id_slopes, aes(condition_label, slope, color = subset)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(ymin = slope - 1.96 * se, ymax = slope + 1.96 * se),
                  position = position_dodge(width = 0.5), size = 0.7) +
  scale_color_manual(name = "", values = c("#004D40", "#1E88E5", "#EB7AA3")) +
  scale_x_discrete(name = "") +
  scale_y_continuous(name = "Acceptability change per extra token") +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "top")
print(p_id_slopes)
ggsave("figures/en_exchange_identification_length_slopes.png", p_id_slopes, width = 9, height = 5.5)


## ---- 4.10 The within-item test: does the SAME error cost less in a longer sentence? ----
## 4.9 rules out a length main effect by its sign. This is the sharper version.
##
## Within one item, every condition has the SAME word count -- a swap only
## reorders. So take the acceptability penalty
##     penalty = mean(original) - mean(swapped)
## item by item. Any main effect of length subtracts the same amount from both
## terms and cancels exactly. So:
##   * length aversion       -> penalty flat in length (slope 0)
##   * binary grammaticality -> penalty flat in length (slope 0)
##   * gradient (raw dl)     -> penalty flat in length (slope 0)
##   * communication         -> penalty SHRINKS with length (slope < 0): the same
##                              displacement is a smaller share of the message.
## Every account except communication makes the same prediction here, so a
## non-zero slope is evidence for communication specifically and against all
## three alternatives at once.
cat("\n\n#############################################################\n")
cat("## 4.10 Within-item swap penalty as a function of length\n")
cat("#############################################################\n")

penalty_of <- function(d, item_col) {
  item_col <- rlang::ensym(item_col)
  d %>%
    filter(condition %in% c("original", "swap_1", "swap_3")) %>%
    group_by(!!item_col, n_words, condition) %>%
    summarize(m = mean(response_num), .groups = "drop") %>%
    pivot_wider(names_from = condition, values_from = m) %>%
    filter(!is.na(original), !is.na(swap_1) | !is.na(swap_3)) %>%
    transmute(item = !!item_col, n_words,
              `1 swap` = original - swap_1,
              `3 swaps` = original - swap_3) %>%
    pivot_longer(c(`1 swap`, `3 swaps`), names_to = "condition_label",
                 values_to = "penalty") %>%
    filter(!is.na(penalty))
}

## Sanity check: length really is constant within an item.
cat("\nitems whose conditions disagree on word count (should be 0):\n")
crit %>% filter(language %in% c("en", "en_exp"), dl_distance <= 3) %>%
  group_by(language, item) %>% summarize(k = n_distinct(n_words), .groups = "drop") %>%
  filter(k > 1) %>% nrow() %>% cat("  ", ., "\n")

## (a) The original English experiment.
en_pen <- crit %>%
  filter(language == "en", dl_distance <= 3) %>%
  penalty_of(item)

cat("\n---- en: penalty vs length ----\n")
en_pen %>% group_by(condition_label) %>%
  group_modify(~ {
    m <- lm(penalty ~ n_words, data = .x)
    s <- summary(m)$coefficients
    tibble(items = nrow(.x), intercept = round(s[1, 1], 3),
           slope = round(s[2, 1], 4), se = round(s[2, 2], 4),
           p = signif(s[2, 4], 3))
  }) %>% print()

## Trial-level confirmation: the condition x length interaction is the same
## quantity, estimated with participant and item random effects.
en_int_dat <- crit %>%
  filter(language == "en", dl_distance <= 3,
         condition %in% c("original", "swap_1", "swap_3")) %>%
  mutate(condition = factor(condition, levels = c("original", "swap_1", "swap_3")),
         n_words_c = n_words - mean(n_words))
m_en_int <- lmer(response_num ~ condition * n_words_c +
                   (1 | participant_id) + (1 | item),
                 data = en_int_dat, REML = FALSE, control = lmer_ctrl)
cat("\n---- en: condition x length interaction (trial level) ----\n")
cat("positive interaction = the swap penalty shrinks as the sentence gets longer\n")
print(round(summary(m_en_int)$coefficients, 4))

## (b) The controlled English experiment. Here length is manipulated WITHIN a
## base item -- short, medium and long versions of the same content -- so this
## repeats the test with the topic held constant as well as the word count.
exp_pen <- crit %>%
  filter(language == "en_exp", dl_distance <= 3) %>%
  penalty_of(item)

cat("\n---- en_exp (controlled): penalty vs length ----\n")
exp_pen %>% group_by(condition_label) %>%
  group_modify(~ {
    m <- lm(penalty ~ n_words, data = .x)
    s <- summary(m)$coefficients
    tibble(items = nrow(.x), intercept = round(s[1, 1], 3),
           slope = round(s[2, 1], 4), se = round(s[2, 2], 4),
           p = signif(s[2, 4], 3))
  }) %>% print()

exp_int_dat <- crit %>%
  filter(language == "en_exp", dl_distance <= 3,
         condition %in% c("original", "swap_1", "swap_3")) %>%
  mutate(base_item = factor(str_remove(as.character(item), "[sml]$")),
         condition = factor(condition, levels = c("original", "swap_1", "swap_3")),
         n_words_c = n_words - mean(n_words))
m_exp_int <- lmer(response_num ~ condition * n_words_c +
                    (1 | participant_id) + (1 | base_item),
                  data = exp_int_dat, REML = FALSE, control = lmer_ctrl)
cat("\n---- en_exp: condition x length interaction, content held constant ----\n")
cat("random effect is the BASE item, so short/medium/long versions of the same\n")
cat("passage are compared against each other\n")
print(round(summary(m_exp_int)$coefficients, 4))

## Figure: the penalty against length in both experiments.
p_penalty <- bind_rows(
  en_pen %>% mutate(exp = "English (original)"),
  exp_pen %>% mutate(exp = "English (controlled)")) %>%
  mutate(exp = factor(exp, levels = c("English (original)", "English (controlled)"))) %>%
  ggplot(aes(n_words, penalty, color = condition_label, fill = condition_label)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(alpha = 0.35) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "tp")) +
  facet_wrap(~ exp) +
  scale_color_manual(name = "Condition", values = error_colors[c("1 swap", "3 swaps")],
                     aesthetics = c("colour", "fill")) +
  scale_x_continuous(name = "Sentence length (tokens)") +
  scale_y_continuous(name = "Acceptability cost of the swap\n(original - swapped, same item)") +
  guides(color = guide_legend(nrow = 1), fill = guide_legend(nrow = 1)) +
  theme(legend.position = "top") +
  ggtitle("The same error costs less in a longer sentence")
print(p_penalty)
ggsave("figures/en_swap_penalty_by_length.png", p_penalty, width = 10, height = 5)
