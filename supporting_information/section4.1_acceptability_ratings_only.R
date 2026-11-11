## =============================================================================
## SI Section 4.1: Acceptability ratings only (constructed sentences)
##
## en_artificial_acceptability_ratings. The 243 constructed sentences (27 base
## items x short, medium and long x original, 1 swap and 3 swaps; see
## section4_setup.R), and every participant sees every base item once.
## Participants only rate each sentence (1-7), with the prompt and scale of the
## main experiments.
##
## The main text's analyses, with sentence surprisal as Information(s) and the
## specifications of the main analysis's controlled English experiment:
##   (1) the five models (Binary, Gradient, Communication, SLOR, MORCELA) as
##       shape-constrained GAMs on item x condition means, compared on their
##       marginal log-likelihood;
##   (2) the linear mixed model of acceptability on DL distance x surprisal,
##       over trials;
##   (3) the surprisal slope among the original (error-free) sentences;
##   (4) the length check, in the CSV only: Communication with length in place
##       of surprisal.
## There is no second task in this study, so no trials can be singled out as
## correct: every analysis uses all trials.
##
## Mixed models have the maximal random-effects structure and are fit with
## bobyqa; singular fits are kept. ../analysis/main_analysis.R and the other
## section 4 studies fit the linear mixed model of (2) with nlminb, but here
## nlminb stops short (see (2)). For every mixed model the CSVs record the
## optimizer, every message and warning printed while fitting it (lme4's and
## lmerTest's) and lme4's gradient and Hessian checks; see (2).
##
## Figure: figures/s4_acceptability_ratings_only_all_trials.png/.pdf, the main
## text's human-judgment curves for this study
## Numbers: outputs/s4_acceptability_ratings_only_*.csv
## The log-likelihoods are drawn together with the other two studies' by
## section4_model_comparison_figure.R.
## =============================================================================

## Print warnings as they happen, so that they land next to the step that raised
## them in the log
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

## Libraries, plotting theme, the main analysis's helper functions
## (marginal_logLik(), morcela_weights(), morcela_score(), from setup.R), the
## figures' colours, and the table of the 243 sentences
source("section4_setup.R")

task_name <- "acceptability_ratings_only"

## ---- Responses ----
responses <- read_csv("../data/en_artificial_acceptability_ratings/en_artificial_acceptability_ratings_data.csv",
                      show_col_types = FALSE)

cat("\n---- export integrity ----\n")
cat("participants in the file:", n_distinct(responses$participant_id), "\n")
responses %>% count(participant_id, name = "trials") %>% count(trials, name = "participants") %>% print()
cat("(27 trials = complete; Surveyor writes no row for a skipped question)\n")
responses %>% distinct(participant_id, native_language) %>% count(native_language) %>% print()

## The study is of native English speakers: drop anyone who answered No.
## A blank answer is kept, as in the other two studies.
non_native <- responses %>%
  distinct(participant_id, native_language) %>%
  filter(native_language %in% "No") %>%
  pull(participant_id)
cat("excluded as non-native:", length(non_native), "participant(s)\n")

trials <- responses %>%
  filter(!participant_id %in% non_native, experiment == "exp_b") %>%
  transmute(participant_id, item = as.integer(item), cond,
            acceptability = as.numeric(substr(response, 1, 1))) %>%
  left_join(sentences, by = c("item", "cond"))

stopifnot(!anyNA(trials))
cat("participants analysed:", n_distinct(trials$participant_id), "| trials:", nrow(trials), "\n")

cells <- trials %>% count(item, cond)
cat("item x condition cells seen:", nrow(cells), "of", nrow(sentences),
    "| ratings per cell: median", median(cells$n),
    "range", min(cells$n), "-", max(cells$n), "\n")

cat("\n---- mean rating by condition and length ----\n")
trials %>%
  group_by(condition_label, length_level) %>%
  summarize(mean_rating = round(mean(acceptability), 2), .groups = "drop") %>%
  pivot_wider(names_from = length_level, values_from = mean_rating) %>%
  print()

## ---- The trials analysed ----
## One long table, with every analysis below run on each subset of it. This
## study has one subset: all trials.
analysis_trials <- bind_rows("All trials" = trials, .id = "subset")

sample_sizes <- analysis_trials %>%
  group_by(subset) %>%
  summarize(n_participants = n_distinct(participant_id), n_trials = n())
sample_sizes

## ---- (1) The five models on item x condition means ----
## As in the main analysis, the models are fit to the mean rating of each
## sentence (one base item at one length in one version) over the trials
## analysed. SLOR is computed as in the main analysis.
item_means <- analysis_trials %>%
  group_by(subset, item, cond, length_level, condition_label, n_words,
           sum_surprisal, dl_distance, weighted_dl_distance,
           weighted_dl_distance_surprisal, grammatical_numeric,
           log_p_sentence, log_p_unigram) %>%
  summarize(response_num = mean(acceptability), n_ratings = n(), .groups = "drop") %>%
  mutate(
    item = as.factor(item),
    slor = (log_p_sentence - log_p_unigram) / n_words,
    morcela = NA_real_
  )

## MORCELA's two weights are estimated on these same means, separately for each
## subset, by the main analysis's recipe
morcela_weights_by_subset <- tibble()
for (subset_iter in unique(item_means$subset)) {
  in_subset <- item_means$subset == subset_iter
  weights <- morcela_weights(item_means[in_subset, ])
  item_means$morcela[in_subset] <- morcela_score(item_means[in_subset, ], weights)
  morcela_weights_by_subset <- bind_rows(
    morcela_weights_by_subset,
    tibble(task = task_name, subset = subset_iter,
           beta = weights$beta, gamma = weights$gamma)
  )
}
morcela_weights_by_subset

## The models, with the specifications of the main analysis's controlled
## English branch: the base item takes a random intercept and a random slope
## for every predictor.
##   - Binary: grammatical vs. ungrammatical
##   - Gradient: a decreasing smooth of the number of errors (DL distance)
##   - Communication: a decreasing smooth of the errors per bit of information
##     (DL distance / surprisal of the intended sentence) plus a decreasing
##     smooth of that surprisal
##   - SLOR and MORCELA: an increasing smooth of the score
## The last model is the length check: Communication with length (in words) as
## the information measure, with the random effects of SI section 3's
## controlled English fit (section3_length.R).
models_by_subset <- list()
for (subset_iter in unique(item_means$subset)) {
  print(paste0("Fitting models for ", subset_iter))

  item_means_subset <- item_means %>% filter(subset == subset_iter)

  grammatical <- scam(
    response_num ~ grammatical_numeric +
      s(item, bs="re") + s(item, grammatical_numeric, bs="re"),
    data=item_means_subset
  )

  distance <- scam(
    response_num ~ s(dl_distance, bs="mpd") +
      s(item, bs="re") + s(item, dl_distance, bs="re"),
    data=item_means_subset
  )

  weighted_distance <- scam(
    response_num ~ s(weighted_dl_distance_surprisal, bs="mpd") + s(sum_surprisal, bs="mpd") +
      s(item, bs="re") + s(item, weighted_dl_distance_surprisal, bs="re") +
      s(item, sum_surprisal, bs="re"),
    data=item_means_subset
  )

  slor <- scam(
    response_num ~ s(slor, bs="mpi") +
      s(item, bs="re") + s(item, slor, bs="re"),
    data=item_means_subset
  )

  morcela <- scam(
    response_num ~ s(morcela, bs="mpi") +
      s(item, bs="re") + s(item, morcela, bs="re"),
    data=item_means_subset
  )

  weighted_distance_length <- scam(
    response_num ~ s(weighted_dl_distance, bs="mpd") + s(n_words, bs="mpd") +
      s(item, bs="re") + s(item, weighted_dl_distance, bs="re") +
      s(item, n_words, bs="re"),
    data=item_means_subset
  )

  models_by_subset[[subset_iter]] <- list(
    grammatical = grammatical,
    distance = distance,
    weighted_distance = weighted_distance,
    slor = slor,
    morcela = morcela,
    weighted_distance_length = weighted_distance_length
  )
}

## Evaluate convergence (one column per subset)
sapply(models_by_subset, sapply, `[[`, "conv")

## Marginal log-likelihoods, with the random effects integrated out, as in the
## main analysis
model_fits <- tibble()
for (subset_iter in names(models_by_subset)) {
  model_fits <- bind_rows(model_fits, tibble(
    subset = subset_iter,
    fit = names(models_by_subset[[subset_iter]]),
    log_likelihood = vapply(models_by_subset[[subset_iter]], marginal_logLik, numeric(1)),
    converged = sapply(models_by_subset[[subset_iter]], `[[`, "conv")
  ))
}

## The comparison with surprisal as the information measure is the five models.
## The length check puts Communication's length fit in place of its surprisal
## fit; the other four models have no information term, so they are the same
## fits in both.
likelihoods <- bind_rows(
  Surprisal = model_fits %>% filter(fit != "weighted_distance_length"),
  Length = model_fits %>%
    filter(fit != "weighted_distance") %>%
    mutate(fit = recode(fit, weighted_distance_length = "weighted_distance")),
  .id = "information"
) %>%
  mutate(
    task = task_name,
    model = case_when(
      fit == "grammatical" ~ "Binary",
      fit == "distance" ~ "Gradient",
      fit == "weighted_distance" ~ "Communication",
      fit == "slor" ~ "SLOR",
      fit == "morcela" ~ "MORCELA"
    ) %>%
      factor(levels = MODEL_LEVELS)
  ) %>%
  group_by(subset, information) %>%
  mutate(
    max_likelihood = max(log_likelihood),
    likelihood_diff = log_likelihood - max_likelihood
  ) %>%
  ungroup() %>%
  left_join(sample_sizes, by = "subset") %>%
  left_join(count(item_means, subset, name = "n_cells"), by = "subset") %>%
  arrange(subset, desc(information), model) %>%
  select(task, subset, information, model, log_likelihood, max_likelihood,
         likelihood_diff, converged, n_participants, n_trials, n_cells)

likelihoods %>%
  select(subset, information, model, log_likelihood, likelihood_diff, converged) %>%
  mutate(log_likelihood = round(log_likelihood, 1),
         likelihood_diff = round(likelihood_diff, 1)) %>%
  print(n = Inf)

## ---- (2) The linear mixed model of the main text ----
## Acceptability on DL distance x surprisal over trials, with the random effects
## of the main analysis's controlled English fit, which are the maximal ones:
## by-participant and by-item (base item) slopes for both predictors and their
## interaction. Both predictors are z-scored within the trials analysed, as in
## the main analysis, and the p-values are lmerTest's (Satterthwaite).
##
## Convergence. Each fit runs inside capture.output(type = "message"), which
## records every message and warning printed while the model is fit (with
## warn = 1, set above, warnings print as they happen): lme4's messages, and
## lmerTest's warning when the Hessian of the variance parameters has a
## negative eigenvalue. They go to the convergence column of the CSV. lme4
## skips its gradient and Hessian checks on singular fits, which all of these
## are, so those two checks are computed by hand below.
lmer_fits <- list()
lmer_messages <- list()
for (subset_iter in unique(analysis_trials$subset)) {
  print(paste0("Fitting the linear mixed model for ", subset_iter))
  lmer_messages[[subset_iter]] <- capture.output(
    lmer_fits[[subset_iter]] <- lmer(
      acceptability ~ dl_distance * sum_surprisal +
        (dl_distance * sum_surprisal|participant_id) +
        (dl_distance * sum_surprisal|item),
      data=analysis_trials %>%
        filter(subset == subset_iter) %>%
        mutate(
          item = as.factor(item),
          sum_surprisal = scale(sum_surprisal)[,1],
          dl_distance = scale(dl_distance)[,1]
        ),
      ## bobyqa rather than the main text's nlminb: in lme4's allFit(), nlminbwrap
      ## stops 26.6 REML log-likelihood units short of the optimum that bobyqa
      ## and NLopt's BOBYQA reach (within 1.2e-4), with a negative eigenvalue of
      ## the Hessian, and its fixed effects differ (surprisal slope .121 against
      ## .185).
      control = lmerControl(
        optimizer ="bobyqa",
        optCtrl = list(maxfun = 1e5)
      )
    ),
    type = "message"
  )
}

lapply(lmer_fits, summary)

## Save the fixed effects, with the optimizer and its return code (0 = OK),
## whether the fit is singular, the messages printed while fitting ("none" if
## there were none), and lme4's gradient and Hessian checks, computed as lme4's
## checkConv() does from the gradient and Hessian that lme4 stores at the
## optimum (fit@optinfo$derivs):
##   - max_gradient: the largest absolute gradient, scaled by the Hessian; lme4
##     warns above 0.002. NA when the Hessian is not positive definite, where
##     lme4 says "unable to evaluate scaled gradient".
##   - negative_hessian_eigenvalues: the Hessian's eigenvalues below -1e-6;
##     lme4 warns of a degenerate Hessian when there are any.
lmm_dl_by_surprisal <- tibble()
for (subset_iter in names(lmer_fits)) {
  coefficients <- summary(lmer_fits[[subset_iter]])$coefficients
  derivs <- lmer_fits[[subset_iter]]@optinfo$derivs
  scaled_gradient <- try(solve(chol(derivs$Hessian), derivs$gradient), silent = TRUE)
  lmm_dl_by_surprisal <- bind_rows(lmm_dl_by_surprisal, tibble(
    task = task_name,
    subset = subset_iter,
    term = rownames(coefficients),
    estimate = coefficients[, "Estimate"],
    se = coefficients[, "Std. Error"],
    df = coefficients[, "df"],
    t = coefficients[, "t value"],
    p = coefficients[, "Pr(>|t|)"],
    singular = isSingular(lmer_fits[[subset_iter]]),
    optimizer = lmer_fits[[subset_iter]]@optinfo$optimizer,
    optimizer_code = lmer_fits[[subset_iter]]@optinfo$conv$opt,
    convergence = ifelse(length(lmer_messages[[subset_iter]]) > 0,
                         paste(lmer_messages[[subset_iter]], collapse = "; "), "none"),
    max_gradient = ifelse(inherits(scaled_gradient, "try-error"), NA_real_,
                          max(pmin(abs(scaled_gradient), abs(derivs$gradient)))),
    negative_hessian_eigenvalues = sum(eigen(derivs$Hessian, symmetric = TRUE,
                                             only.values = TRUE)$values < -1e-6)
  ))
}
lmm_dl_by_surprisal <- lmm_dl_by_surprisal %>%
  left_join(sample_sizes, by = "subset")

lmm_dl_by_surprisal %>%
  select(subset, term, estimate, se, t, p, singular, optimizer, optimizer_code,
         convergence, max_gradient, negative_hessian_eigenvalues) %>%
  print(n = Inf, width = Inf)

## ---- (3) The surprisal slope among the error-free sentences ----
## The main analysis's controlled English fit: the original sentences only,
## surprisal z-scored within them, and by-participant and by-item surprisal
## slopes (the maximal structure). Convergence is recorded as in (2).
error_free_fits <- list()
error_free_messages <- list()
for (subset_iter in unique(analysis_trials$subset)) {
  print(paste0("Fitting the error-free slope for ", subset_iter))
  error_free_messages[[subset_iter]] <- capture.output(
    error_free_fits[[subset_iter]] <- lmer(
      acceptability ~ sum_surprisal +
        (sum_surprisal | participant_id) +
        (sum_surprisal | item),
      data=analysis_trials %>%
        filter(subset == subset_iter, dl_distance == 0) %>%
        mutate(
          item = as.factor(item),
          sum_surprisal = scale(sum_surprisal)[,1]
        ),
      control = lmerControl(
        optimizer ="bobyqa",
        optCtrl = list(maxfun = 1e5)
      )
    ),
    type = "message"
  )
}

lapply(error_free_fits, summary)

## beta is per SD of surprisal; sd_surprisal_bits (the SD among the error-free
## trials analysed) converts it to per bit
error_free_sizes <- analysis_trials %>%
  filter(dl_distance == 0) %>%
  group_by(subset) %>%
  summarize(sd_surprisal_bits = sd(sum_surprisal),
            n_participants = n_distinct(participant_id),
            n_trials = n())

## The convergence columns are those of (2)
error_free_surprisal_slopes <- tibble()
for (subset_iter in names(error_free_fits)) {
  coefficients <- summary(error_free_fits[[subset_iter]])$coefficients
  derivs <- error_free_fits[[subset_iter]]@optinfo$derivs
  scaled_gradient <- try(solve(chol(derivs$Hessian), derivs$gradient), silent = TRUE)
  error_free_surprisal_slopes <- bind_rows(error_free_surprisal_slopes, tibble(
    task = task_name,
    subset = subset_iter,
    beta = coefficients["sum_surprisal", "Estimate"],
    se = coefficients["sum_surprisal", "Std. Error"],
    df = coefficients["sum_surprisal", "df"],
    t = coefficients["sum_surprisal", "t value"],
    p = coefficients["sum_surprisal", "Pr(>|t|)"],
    singular = isSingular(error_free_fits[[subset_iter]]),
    optimizer = error_free_fits[[subset_iter]]@optinfo$optimizer,
    optimizer_code = error_free_fits[[subset_iter]]@optinfo$conv$opt,
    convergence = ifelse(length(error_free_messages[[subset_iter]]) > 0,
                         paste(error_free_messages[[subset_iter]], collapse = "; "), "none"),
    max_gradient = ifelse(inherits(scaled_gradient, "try-error"), NA_real_,
                          max(pmin(abs(scaled_gradient), abs(derivs$gradient)))),
    negative_hessian_eigenvalues = sum(eigen(derivs$Hessian, symmetric = TRUE,
                                             only.values = TRUE)$values < -1e-6)
  ))
}
error_free_surprisal_slopes <- error_free_surprisal_slopes %>%
  left_join(error_free_sizes, by = "subset")

error_free_surprisal_slopes %>%
  select(subset, beta, se, t, p, sd_surprisal_bits, singular, optimizer,
         optimizer_code, convergence, max_gradient, negative_hessian_eigenvalues) %>%
  print(width = Inf)

## ---- Save the tables ----
write_csv(likelihoods, paste0("outputs/s4_", task_name, "_likelihoods.csv"))
write_csv(lmm_dl_by_surprisal, paste0("outputs/s4_", task_name, "_lmm_dl_by_surprisal.csv"))
write_csv(error_free_surprisal_slopes,
          paste0("outputs/s4_", task_name, "_error_free_surprisal_slopes.csv"))
write_csv(morcela_weights_by_subset, paste0("outputs/s4_", task_name, "_morcela_weights.csv"))

## ---- Figure: acceptability ----
## The human-judgment curves of the main text's qualitative-predictions figure,
## with its plotting code (p_baseline and p_gam in analysis/main_analysis.R),
## for all trials of this study: it has no second task. Points are sentence
## means (one item at one length in one condition); lines are GAM smooths of
## those means over surprisal, one per DL distance, with 95% confidence bands.
## Like the main figure, it is drawn at its print size, one PNAS column (8.7 cm)
## wide, in the main figure's 8-pt print theme (theme_print(8, 0.3) there).
## Where it differs from the main code: it has no "Human Judgments" title (there
## are no model columns to tell it from) and shows its own legend, at the
## bottom like the main figure's; the legend is centred on the whole plot, as
## the main figure's collected legend is, so that it fits at this width; and it
## leaves out "4+", since no sentence here has a DL distance above 3.
plot_data_by_item <- analysis_trials %>%
  filter(subset == "All trials") %>%
  group_by(subset, item, length_level, condition_label, sum_surprisal, dl_distance) %>%
  summarize(response = mean(acceptability), .groups = "drop") %>%
  mutate(dl_distance_bucket = bucket_dl(dl_distance))

p_acceptability <- ggplot(
  plot_data_by_item,
  aes(x = sum_surprisal, y = response, color = dl_distance_bucket)
) +
  geom_point(alpha = 0.3, size = 0.6, stroke = 0) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, bs = "tp"),
    se = TRUE,
    linewidth = 0.5,
    show.legend = FALSE
  ) +
  facet_wrap(~subset, ncol = 1, strip.position = "left") +
  labs(x = "Sentence surprisal (bits)") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        legend.key = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        strip.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5)) +
  theme(
    legend.position = "bottom",
    legend.location = "plot",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.spacing = unit(2, "pt"),
    strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1),
    panel.spacing = unit(4, "pt")
  ) +
  scale_y_continuous(
    limits = YLIM,
    breaks = seq(1, 7, 2),
    labels = seq(1, 7, 2),
    name = "Acceptability Ratings (1-7 scale)"
  ) +
  scale_color_manual(
    name   = "Damerau-Levenshtein Distance",
    values = DL_COLORS,
    breaks = c("0", "1-2", "3"),
    limits = DL_BUCKET_LEVELS,
    drop   = FALSE,
    guide  = guide_legend(
      override.aes = list(alpha = 1, linewidth = 0, size = 1.6)
    )
  )

if (interactive()) print(p_acceptability)

ggsave(paste0("figures/s4_", task_name, "_all_trials.png"),
       p_acceptability, width = 3.42, height = 2.6, units = "in", dpi = 300)

## Save as pdf, with the fonts embedded (cairo_pdf)
ggsave(paste0("figures/s4_", task_name, "_all_trials.pdf"),
       p_acceptability, width = 3.42, height = 2.6, units = "in", device = cairo_pdf)

print("Done.")
