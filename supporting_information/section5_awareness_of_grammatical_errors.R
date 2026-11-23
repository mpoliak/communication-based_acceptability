## =============================================================================
## SI Section 5: Probing participants' awareness of grammatical errors
##
## One script for the four follow-up experiments of SI section 5, in the order
## of the SI:
##   5.1  Binary error detection + acceptability ratings, on the original
##        materials of the English experiment (../data/en_binary_error_detection)
##   5.2  Acceptability ratings only, on the constructed sentences
##   5.3  Binary error detection + acceptability ratings, constructed sentences
##   5.4  Exchange identification + acceptability ratings, constructed sentences
## and then the figure of the five models' log-likelihoods in the three studies
## of constructed sentences (SI Figure S11).
##
## 5.1 starts from setup.R; 5.2-5.4 and the figure share the set-up of the
## constructed sentences, which follows 5.1. Before each part, the objects of
## the previous part are removed, so that no part uses another part's objects.
## Outputs: outputs/s5_* and figures/s5_*.
## Run from inside supporting_information/:
##   mkdir -p ../logs
##   LC_ALL=en_US.UTF-8 Rscript section5_awareness_of_grammatical_errors.R \
##     > ../logs/section5_awareness_of_grammatical_errors.log 2>&1
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


## #############################################################################
## 5.1  Binary error detection + acceptability ratings (original sentences)
## #############################################################################
## =============================================================================
## SI Section 5.1: Binary error detection + acceptability ratings (original
## English materials)
##
## en_binary_error_detection. The 72 items of the main text's English
## experiment (sampled from Universal Dependencies), each in that experiment's
## four conditions: original, 1 swap (swap_1), 3 swaps (swap_3) and shuffled
## (random). The 288 sentences are those of the English experiment, character
## for character (checked below against ../data/en/materials_en.csv). Every
## participant sees every item once, 18 items in each condition. On each trial
## the participant first answers "Are there any errors in this sentence?"
## (Yes/No) and then, on the next screen, rates how acceptable the sentence is
## (1-7). The instructions in demographics_en_binary_error_detection.csv
## describe the opposite order (rate first, then the error question); the data
## show the order above (checked below).
##
## Judged correctly: the answer to the binary error detection question matches
## literal_response, i.e. "No" (no error) for an original and "Yes" for a
## swapped or shuffled sentence.
##
## The main text's analyses, with sentence surprisal as Information(s) and the
## specifications of the main analysis's English experiment (en), each on all
## trials and again on the trials judged correctly only:
##   (1) the five models (Binary, Gradient, Communication, SLOR, MORCELA) as
##       shape-constrained GAMs on item x condition means, compared on their
##       marginal log-likelihood;
##   (2) the linear mixed model of acceptability on DL distance x surprisal,
##       over trials;
##   (3) the surprisal slope among the original (error-free) sentences;
##   (4) the length check, in the CSV only: Communication with length in place
##       of surprisal.
## Restricting to the trials judged correctly selects on a response that is
## correlated with the rating, so all trials are the baseline to read the
## correct-trial results against.
##
## Unlike the constructed sentences of sections 5.2-5.4, these items have no
## length versions: surprisal and length are properties of the item, constant
## across its conditions, as in the main analysis's English experiment. So the
## models take that experiment's random effects: a by-item intercept in the
## GAMs, by-item DL distance slopes (and no surprisal slopes) in (2), and a
## by-item intercept in (3). As in the main analysis, the models use the
## sentences with a DL distance of at most 3, which drops the shuffled
## condition except for its 3 sentences with a DL distance of 2 or 3 (items 15,
## 16 and 23); the five models with the shuffled condition included are SI
## section 2's analysis of the main experiments. Binary error detection is
## reported for all four conditions, shuffled included.
##
## Mixed models have the maximal random-effects structure and are fit with
## bobyqa; singular fits are kept. ../analysis/main_analysis.R fits the linear
## mixed model of (2) with nlminb, which reaches the same optimum here (see
## (2)). For every mixed model the CSVs record the optimizer,
## every message and warning printed while fitting it (lme4's and lmerTest's)
## and lme4's gradient and Hessian checks.
##
## Figures: figures/s5_original_binary_error_detection_yes_responses.png
## (proportion of "yes" answers, over all binary error detection judgments) and
## figures/s5_original_binary_error_detection_correct_trials.png
## (the main text's human-judgment curves, on the trials judged correctly)
## Numbers: outputs/s5_original_binary_error_detection_*.csv
## =============================================================================

## Libraries, plotting theme, the main analysis's helper functions
## (marginal_logLik(), morcela_weights(), morcela_score()), the five models'
## colours, shapes and labels, the DL-distance colours and buckets, and the main
## experiments' data (crit), all from setup.R
source("setup.R")

task_name <- "original_binary_error_detection"

## ---- Colours ----
## The four conditions take the colours of their DL distances in the main
## text's figures: the swaps have DL distance 1 and (mostly) 3, and nearly all
## shuffled sentences fall in the grey "4+" bucket
CONDITION_COLORS <- c("Original" = DL_COLORS[["0"]],
                      "1 swap"   = DL_COLORS[["1-2"]],
                      "3 swaps"  = DL_COLORS[["3"]],
                      "Shuffled" = DL_COLORS[["4+"]])

## ---- The 288 sentences ----
## In the materials file, the sentence is on the row of the trial's first
## screen, the binary error detection question (the rating that follows it has
## an empty sentence)
materials <- read_csv("../data/en_binary_error_detection/materials_en_binary_error_detection.csv",
                      show_col_types = FALSE) %>%
  filter(experiment == "exp_b_comprehension") %>%
  transmute(item = as.integer(item), cond = condition, sentence)

english_materials <- read_csv("../data/en/materials_en.csv", show_col_types = FALSE) %>%
  filter(experiment == "exp_b") %>%
  transmute(item = as.integer(item), cond = condition, english_sentence = sentence)

## Every item appears once in each of the four conditions, with the sentence
## of the English experiment
stopifnot(
  nrow(materials) == 72 * 4,
  nrow(distinct(materials, item, cond)) == 72 * 4,
  all(materials$sentence ==
        left_join(materials, english_materials, by = c("item", "cond"))$english_sentence)
)

## DL distance and length in words, from the stimulus table of the English
## experiment (written by data/preprocess_stimuli.R; full_data.csv takes its
## dl_distance and n_words from it)
english_stimuli <- read_csv("../data/en/stimuli_data.csv", show_col_types = FALSE) %>%
  transmute(item = as.integer(item), cond = condition, dl_distance, n_words)

## The language-model scores, from the one source of the main analysis
## (data/slor_morcela_scores.csv: Goldfish 1000mb, log2, model-initial unigram)
english_scores <- read_csv("../data/slor_morcela_scores.csv",
                           col_types = cols(item = col_character()),
                           show_col_types = FALSE) %>%
  filter(language == "en", unigram_model == "model_initial") %>%
  transmute(item = as.integer(item), cond = condition, log_p_sentence, log_p_unigram)

## Information(s) in the communication account is the surprisal of the intended
## sentence, i.e. of the original: all four conditions of an item share it. SLOR
## and MORCELA instead score the sentence that was actually read, so they take
## each condition's own log probabilities.
original_surprisal <- english_scores %>%
  filter(cond == "original") %>%
  transmute(item, sum_surprisal = -log_p_sentence)

sentences <- materials %>%
  left_join(english_stimuli, by = c("item", "cond")) %>%
  left_join(original_surprisal, by = "item") %>%
  left_join(english_scores, by = c("item", "cond")) %>%
  mutate(
    condition_label = factor(recode(cond,
                                    "original" = "Original",
                                    "swap_1" = "1 swap",
                                    "swap_3" = "3 swaps",
                                    "random" = "Shuffled"),
                             levels = names(CONDITION_COLORS)),
    weighted_dl_distance = dl_distance / n_words,
    weighted_dl_distance_surprisal = dl_distance / sum_surprisal,
    grammatical_numeric = ifelse(dl_distance == 0, .5, -.5)
  ) %>%
  select(item, cond, sentence, condition_label, n_words, sum_surprisal,
         dl_distance, weighted_dl_distance, weighted_dl_distance_surprisal,
         grammatical_numeric, log_p_sentence, log_p_unigram)

## Every sentence found its DL distance, length and scores, and they are the
## values the main analysis uses for the English experiment (crit, from
## setup.R), in every cell that experiment's participants saw
english_crit <- crit %>%
  filter(language == "en") %>%
  transmute(item = as.integer(str_remove(as.character(item), "^en_")), cond = condition,
            dl_distance, n_words, sum_surprisal, log_p_sentence, log_p_unigram) %>%
  distinct()

stopifnot(
  nrow(sentences) == 72 * 4,
  !anyNA(sentences),
  nrow(inner_join(english_crit, sentences, by = names(english_crit))) == nrow(english_crit)
)
cat("cells of the English experiment matched:", nrow(english_crit), "of", nrow(sentences), "\n")

cat("DL distance by condition (number of sentences):\n")
sentences %>% count(condition_label, dl_distance) %>% print(n = Inf)

## The 72 original sentences: their length in words and their surprisal in bits
## (the Information(s) of all four conditions)
cat("Length (words) and surprisal (bits) of the original sentences:\n")
sentences %>%
  filter(condition_label == "Original") %>%
  summarize(n_sentences = n(),
            min_words = min(n_words), median_words = median(n_words),
            max_words = max(n_words),
            min_surprisal = min(sum_surprisal), median_surprisal = median(sum_surprisal),
            max_surprisal = max(sum_surprisal)) %>%
  print(width = Inf)

## ---- Responses: one row per screen: the binary error detection question (exp_b_comprehension), then the rating (exp_b) ----
responses <- read_csv("../data/en_binary_error_detection/data_en_binary_error_detection.csv",
                      show_col_types = FALSE)

cat("\n---- export integrity ----\n")
cat("participants in the file:", n_distinct(responses$participant_id), "\n")
responses %>% distinct(participant_id, native_language) %>% count(native_language) %>% print()
responses %>% distinct(participant_id, country_of_birth) %>% count(country_of_birth) %>% print()
responses %>% distinct(participant_id, other_languages) %>% count(other_languages) %>% print()

## The study is of native English speakers: drop anyone who answered No (nobody
## did here). A blank answer is kept, as in the section 5.2-5.4 studies.
non_native <- responses %>%
  distinct(participant_id, native_language) %>%
  filter(native_language %in% "No") %>%
  pull(participant_id)
cat("excluded as non-native:", length(non_native), "participant(s)\n")
responses <- responses %>% filter(!participant_id %in% non_native)

## Sections 5.3 and 5.4 exclude the people who took both of those studies,
## which show the same sentences: the studies of constructed sentences share
## one set of pseudonyms (en_artificial_<n>), so the same id in two of them is
## the same person. The participant ids here are pseudonyms given within this
## study (en_awareness_<n>), as are those of the English experiment (en_<n>),
## and the raw Prolific ids are not kept, so who else took part in another
## study cannot be checked, and nobody is excluded for it.

judgments <- responses %>%
  filter(experiment == "exp_b_comprehension") %>%
  transmute(participant_id, item = as.integer(item), cond,
            judgment_order = presentation_order,
            ## Judged correctly: "No" for an original, "Yes" for a swapped or
            ## shuffled sentence
            correct = response == literal_response,
            ## Reported an error: answered "Yes"
            reported_error = response == "Yes")
## A "yes" is correct for every sentence except an original
stopifnot(all(judgments$reported_error == (judgments$correct != str_detect(judgments$cond, "^original"))))

ratings <- responses %>%
  filter(experiment == "exp_b") %>%
  transmute(participant_id, item = as.integer(item), cond,
            rating_order = presentation_order,
            acceptability = as.numeric(substr(response, 1, 1)))

responses %>% count(participant_id, name = "rows") %>% count(rows, name = "participants") %>% print()
cat("(144 rows = 72 sentences x 2 questions; Surveyor writes no row for a skipped question)\n")
cat("judgments without a rating:",
    nrow(anti_join(judgments, ratings, by = c("participant_id", "item", "cond"))),
    "| ratings without a judgment:",
    nrow(anti_join(ratings, judgments, by = c("participant_id", "item", "cond"))), "\n")

## Condition assignment: every participant sees each item at most once, in one
## condition, and the 62 who answered all 72 binary error detection questions
## saw 18 items in each condition. Surveyor rotates the item-to-condition
## assignment by a seed from the participant id, so every participant has their
## own assignment and the item x condition cells are not equally filled.
stopifnot(!anyDuplicated(judgments[c("participant_id", "item")]))
judgments %>%
  count(participant_id, cond) %>%
  pivot_wider(names_from = cond, values_from = n, values_fill = 0) %>%
  count(original, swap_1, swap_3, random, name = "participants") %>%
  print()
cat("distinct item-to-condition assignments:",
    judgments %>%
      arrange(item) %>%
      group_by(participant_id) %>%
      summarize(assignment = paste(item, cond, collapse = ",")) %>%
      pull(assignment) %>%
      n_distinct(),
    "over", n_distinct(judgments$participant_id), "participants\n")

## Every binary error detection judgment, rated or not, with its sentence: the
## correctness table and figure use all of them, as in section 5.3
judgment_trials <- judgments %>%
  left_join(sentences, by = c("item", "cond"))

## A trial is a sentence with both its binary error detection judgment and its
## rating: the acceptability analyses use these
trials <- judgments %>%
  inner_join(ratings, by = c("participant_id", "item", "cond")) %>%
  left_join(sentences, by = c("item", "cond"))

## The rating is always the screen right after its binary error detection
## question
stopifnot(
  !anyNA(judgment_trials),
  !anyNA(trials),
  all(trials$rating_order == trials$judgment_order + 1)
)
cat("binary error detection judgments:", nrow(judgment_trials),
    "| participants analysed:", n_distinct(trials$participant_id),
    "| trials (judged and rated):", nrow(trials), "\n")

cells <- trials %>% count(item, cond)
cat("item x condition cells seen:", nrow(cells), "of", nrow(sentences),
    "| responses per cell: median", median(cells$n),
    "range", min(cells$n), "-", max(cells$n), "\n")

## ---- Binary error detection: proportion judged correctly ----
## Over all binary error detection judgments, rated or not (as in section 5.3),
## by condition, and over all conditions
correctness <- bind_rows(
  judgment_trials %>%
    group_by(condition_label = as.character(condition_label)) %>%
    summarize(n_trials = n(), proportion_correct = mean(correct),
              proportion_reported_error = mean(reported_error), .groups = "drop"),
  judgment_trials %>%
    summarize(n_trials = n(), proportion_correct = mean(correct),
              proportion_reported_error = mean(reported_error)) %>%
    mutate(condition_label = "All")
) %>%
  mutate(task = task_name, .before = 1) %>%
  arrange(factor(condition_label, levels = c(names(CONDITION_COLORS), "All")))

print(correctness, n = Inf)
write_csv(correctness, paste0("outputs/s5_", task_name, "_correctness.csv"))

## ---- The trials analysed ----
## One long table, with every analysis below run on each subset of it: all
## trials, and the trials judged correctly. It keeps the shuffled sentences,
## for the counts and the figure; the models below take the sentences with a
## DL distance of at most 3, as in the main analysis.
analysis_trials <- bind_rows(
  "All trials" = trials,
  "Judged correctly" = trials %>% filter(correct),
  .id = "subset"
)

cat("\n---- trials by condition (all conditions) ----\n")
analysis_trials %>%
  count(subset, condition_label) %>%
  pivot_wider(names_from = condition_label, values_from = n) %>%
  print()

## The sample sizes of the models (DL distance up to 3)
sample_sizes <- analysis_trials %>%
  filter(dl_distance <= 3) %>%
  group_by(subset) %>%
  summarize(n_participants = n_distinct(participant_id), n_trials = n())
print(sample_sizes)

cat("\n---- mean rating by condition ----\n")
analysis_trials %>%
  group_by(subset, condition_label) %>%
  summarize(mean_rating = round(mean(acceptability), 2), .groups = "drop") %>%
  pivot_wider(names_from = condition_label, values_from = mean_rating) %>%
  print()

## ---- (1) The five models on item x condition means ----
## As in the main analysis, the models are fit to the mean rating of each
## sentence (one item in one condition) with a DL distance of at most 3, over
## the trials analysed. SLOR is computed as in the main analysis.
item_means <- analysis_trials %>%
  filter(dl_distance <= 3) %>%
  group_by(subset, item, cond, condition_label, n_words, sum_surprisal,
           dl_distance, weighted_dl_distance, weighted_dl_distance_surprisal,
           grammatical_numeric, log_p_sentence, log_p_unigram) %>%
  summarize(response_num = mean(acceptability), .groups = "drop") %>%
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
print(morcela_weights_by_subset)

## The models, with the specifications of the main analysis's English
## experiment, which, like every main experiment without length versions of its
## items, gives the item a random intercept only.
##   - Binary: grammatical vs. ungrammatical
##   - Gradient: a decreasing smooth of the number of errors (DL distance)
##   - Communication: a decreasing smooth of the errors per bit of information
##     (DL distance / surprisal of the intended sentence) plus a decreasing
##     smooth of that surprisal
##   - SLOR and MORCELA: an increasing smooth of the score
## The last model is the length check: Communication with length (in words) as
## the information measure, with the random effects of SI section 3's English
## fit (section3_length.R).
models_by_subset <- list()
for (subset_iter in unique(item_means$subset)) {
  print(paste0("Fitting models for ", subset_iter))

  item_means_subset <- item_means %>% filter(subset == subset_iter)

  grammatical <- scam(
    response_num ~ grammatical_numeric +
      s(item, bs="re"),
    data=item_means_subset
  )

  distance <- scam(
    response_num ~ s(dl_distance, bs="mpd") +
      s(item, bs="re"),
    data=item_means_subset
  )

  weighted_distance <- scam(
    response_num ~ s(weighted_dl_distance_surprisal, bs="mpd") + s(sum_surprisal, bs="mpd") +
      s(item, bs="re"),
    data=item_means_subset
  )

  slor <- scam(
    response_num ~ s(slor, bs="mpi") +
      s(item, bs="re"),
    data=item_means_subset
  )

  morcela <- scam(
    response_num ~ s(morcela, bs="mpi") +
      s(item, bs="re"),
    data=item_means_subset
  )

  weighted_distance_length <- scam(
    response_num ~ s(weighted_dl_distance, bs="mpd") + s(n_words, bs="mpd") +
      s(item, bs="re"),
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

## scam's inner fitting loop converged (one column per subset); the
## smoothing-parameter search is checked below
print(sapply(models_by_subset, sapply, `[[`, "conv"))

## Marginal log-likelihoods, with the random effects integrated out, as in the
## main analysis. `converged` is scam's inner fitting loop ($conv), which is
## TRUE even when scam's search for the smoothing parameters stopped short;
## bfgs_convergence is how that search ended ("Full convergence" when the
## gradient of its smoothness score, GCV/UBRE, became small enough).
model_fits <- tibble()
for (subset_iter in names(models_by_subset)) {
  model_fits <- bind_rows(model_fits, tibble(
    subset = subset_iter,
    fit = names(models_by_subset[[subset_iter]]),
    log_likelihood = vapply(models_by_subset[[subset_iter]], marginal_logLik, numeric(1)),
    converged = sapply(models_by_subset[[subset_iter]], `[[`, "conv"),
    bfgs_convergence = sapply(lapply(models_by_subset[[subset_iter]], `[[`, "bfgs.info"), `[[`, "conv"),
    optimizer = "bfgs"
  ))
}

model_fits %>% count(bfgs_convergence) %>% print()

## A fit whose search did not reach full convergence (it reached its 200-step
## limit, or its last step found no lower score) is refit with scam's two other
## optimizers, nlm and optim, as ../analysis/main_analysis.R checks its fits.
## The fit reported is the one with the lowest smoothness score among the
## original and the refits whose search converged; `optimizer` records which
## (bfgs is scam's default). nlm's two messages of success end in "probably
## solution"; optim's is "Successful completion".
for (row_index in which(model_fits$bfgs_convergence != "Full convergence")) {
  m <- models_by_subset[[model_fits$subset[row_index]]][[model_fits$fit[row_index]]]
  refit_fits <- list(
    bfgs = m,
    nlm = scam(m$formula, data = m$model, optimizer = "nlm"),
    optim = scam(m$formula, data = m$model, optimizer = "optim")
  )
  refits <- tibble(
    subset = model_fits$subset[row_index],
    fit = model_fits$fit[row_index],
    optimizer = names(refit_fits),
    search_convergence = c(m$bfgs.info$conv, refit_fits$nlm$nlm.info$conv,
                           refit_fits$optim$optim.info$conv),
    smoothness_score = sapply(refit_fits, `[[`, "gcv.ubre"),
    log_likelihood = vapply(refit_fits, marginal_logLik, numeric(1))
  )
  print(as.data.frame(refits), digits = 7)

  kept <- refits %>%
    filter(optimizer == "bfgs" |
             str_detect(search_convergence, "probably solution|Successful completion")) %>%
    slice_min(smoothness_score, n = 1, with_ties = FALSE)
  print(paste0("Kept: ", kept$optimizer, "; marginal log-likelihood change ",
               round(kept$log_likelihood - model_fits$log_likelihood[row_index], 3)))

  model_fits$optimizer[row_index] <- kept$optimizer
  model_fits$log_likelihood[row_index] <- kept$log_likelihood
  model_fits$converged[row_index] <- refit_fits[[kept$optimizer]]$conv
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
         likelihood_diff, converged, bfgs_convergence, optimizer,
         n_participants, n_trials, n_cells)

likelihoods %>%
  select(subset, information, model, log_likelihood, likelihood_diff, converged,
         optimizer) %>%
  mutate(log_likelihood = round(log_likelihood, 1),
         likelihood_diff = round(likelihood_diff, 1)) %>%
  print(n = Inf)

## ---- (2) The linear mixed model of the main text ----
## Acceptability on DL distance x surprisal over trials with a DL distance of
## at most 3, with the random effects of the main analysis's English fit, which
## are the maximal ones: by-participant slopes for both predictors and their
## interaction, and by-item DL distance slopes (surprisal does not vary within
## an item). Both predictors are z-scored within the trials analysed, as in the
## main analysis, and the p-values are lmerTest's (Satterthwaite).
##
## Convergence. Each fit runs inside capture.output(type = "message"), which
## records every message and warning printed while the model is fit (with
## warn = 1, set above, warnings print as they happen): lme4's messages, and
## lmerTest's warning when the Hessian of the variance parameters has a
## negative eigenvalue. They go to the convergence column of the CSV. lme4
## skips its gradient and Hessian checks on singular fits, so those two checks
## are computed by hand below.
lmer_fits <- list()
lmer_messages <- list()

## Both subsets: bobyqa. ../analysis/main_analysis.R fits the English experiment
## with nlminb (nlminbwrap); on these data the two reach the same optimum (REML
## criterion and fixed effects equal to at least 4 decimals, in both subsets),
## and the bobyqa fits pass every check.
for (subset_iter in unique(analysis_trials$subset)) {
  print(paste0("Fitting the linear mixed model for ", subset_iter))
  lmer_messages[[subset_iter]] <- capture.output(
    lmer_fits[[subset_iter]] <- lmer(
      acceptability ~ dl_distance * sum_surprisal +
        (dl_distance * sum_surprisal|participant_id) +
        (dl_distance|item),
      data=analysis_trials %>%
        filter(subset == subset_iter, dl_distance <= 3) %>%
        mutate(
          item = as.factor(item),
          sum_surprisal = scale(sum_surprisal)[,1],
          dl_distance = scale(dl_distance)[,1]
        ),
      control = lmerControl(
        optimizer ="bobyqa",
        optCtrl = list(maxfun = 1e5)
      )
    ),
    type = "message"
  )
}

print(lapply(lmer_fits, summary))

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
## The main analysis's English fit: the original sentences only, surprisal
## z-scored within them, a by-participant surprisal slope and a by-item
## intercept (surprisal does not vary within an item, so this is the maximal
## structure). Convergence is recorded as in (2).
error_free_fits <- list()
error_free_messages <- list()
for (subset_iter in unique(analysis_trials$subset)) {
  print(paste0("Fitting the error-free slope for ", subset_iter))
  error_free_messages[[subset_iter]] <- capture.output(
    error_free_fits[[subset_iter]] <- lmer(
      acceptability ~ sum_surprisal +
        (sum_surprisal | participant_id) +
        (1 | item),
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

print(lapply(error_free_fits, summary))

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
write_csv(likelihoods, paste0("outputs/s5_", task_name, "_likelihoods.csv"))
write_csv(lmm_dl_by_surprisal, paste0("outputs/s5_", task_name, "_lmm_dl_by_surprisal.csv"))
write_csv(error_free_surprisal_slopes,
          paste0("outputs/s5_", task_name, "_error_free_surprisal_slopes.csv"))
write_csv(morcela_weights_by_subset, paste0("outputs/s5_", task_name, "_morcela_weights.csv"))


## ---- Figure: "yes" responses in binary error detection ----
## Over all binary error detection judgments, rated or not, like the
## correctness table. The y axis is the proportion of "yes" answers to "Are
## there any errors in this sentence?": false alarms for the originals, and
## detected errors for the swapped and shuffled sentences. Points are sentences (one item in one condition). The four sentences of one item share
## their x (the surprisal of the intended sentence), so the points are
## jittered sideways, by up to 2 bits, and never up or down, so that no point
## is drawn outside 0-100%. Lines are logistic regressions of the trial-level
## outcome on surprisal, one per condition; their 95% confidence bands are
## computed on the logit scale, so they stay within 0-100%. They ignore that
## trials are grouped by participant and item. The legend lists the conditions
## from top to bottom in the order of the curves, and calls the swaps
## "exchanges", as the text does.
p_yes_responses <- ggplot(mapping = aes(x = sum_surprisal, color = condition_label,
                                        fill = condition_label)) +
  geom_point(data = judgment_trials %>%
               group_by(condition_label, item, sum_surprisal) %>%
               summarize(value = mean(reported_error), .groups = "drop"),
             aes(y = value), alpha = 0.4,
             position = position_jitter(width = 2, height = 0, seed = 2026)) +
  geom_smooth(data = judgment_trials, aes(y = as.integer(reported_error)),
              method = "glm", method.args = list(family = binomial), formula = y ~ x) +
  scale_color_manual(name = "Condition", values = CONDITION_COLORS,
                     labels = c("Original" = "Original", "1 swap" = "1 exchange",
                                "3 swaps" = "3 exchanges", "Shuffled" = "Shuffled"),
                     aesthetics = c("colour", "fill"),
                     guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous(name = "Sentence surprisal (bits)") +
  scale_y_continuous(name = "Proportion of \"yes\" responses\n(reported an error)",
                     labels = scales::percent, breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1))

if (interactive()) print(p_yes_responses)

ggsave(paste0("figures/s5_", task_name, "_yes_responses.png"),
       p_yes_responses, width = 8, height = 5.5)


## ---- Figure: acceptability on the trials judged correctly ----
## The human-judgment curves of the main text's qualitative-predictions figure,
## with its plotting code (p_baseline and p_gam in analysis/main_analysis.R),
## for the trials judged correctly. Points are sentence means (one item in one
## condition); lines are GAM smooths of those means over surprisal, one per DL
## distance bucket, with 95% confidence bands. As in the main figure, the
## shuffled sentences are drawn too, nearly all in the grey "4+" bucket, though
## the models are fit on DL distances up to 3 only. Like the main figure, it is
## drawn at its print size, one PNAS column (8.7 cm) wide, in the main figure's
## 8-pt print theme (theme_print_8pt there). Where it differs from the main
## code: it has no "Human Judgments" title (there are no model columns to tell
## it from) and shows its own legend, at the bottom like the main figure's; the
## legend is centred on the whole plot, as the main figure's collected legend
## is, so that it fits at this width.
plot_data_by_item <- analysis_trials %>%
  filter(subset == "Judged correctly") %>%
  group_by(subset, item, condition_label, sum_surprisal, dl_distance) %>%
  summarize(response = mean(acceptability), .groups = "drop") %>%
  mutate(dl_distance_bucket = bucket_dl(dl_distance))

p_acceptability <- ggplot(
  plot_data_by_item,
  aes(x = sum_surprisal, y = response, color = dl_distance_bucket)
) +
  geom_point(alpha = 0.5, size = 1.3, stroke = 0) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, bs = "tp"),
    se = TRUE,
    linewidth = 0.5,
    show.legend = FALSE
  ) +
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
    breaks = DL_BUCKET_LEVELS,
    limits = DL_BUCKET_LEVELS,
    drop   = FALSE,
    guide  = guide_legend(
      override.aes = list(alpha = 1, linewidth = 0, size = 1.6)
    )
  )

if (interactive()) print(p_acceptability)

ggsave(paste0("figures/s5_", task_name, "_correct_trials.png"),
       p_acceptability, width = 3.42, height = 2.6, units = "in", dpi = 300)


print("Done.")


## #############################################################################
## Shared set-up of 5.2-5.4 (the constructed sentences)
## #############################################################################
## Start from a clean slate, as the stand-alone scripts did
rm(list = setdiff(ls(), "get_script_dir"))

## =============================================================================
## SI Section 5: set-up shared by the three studies of constructed sentences
##
## Used by parts 5.2, 5.3 and 5.4 and the figure below. It
##   - sources setup.R, for the libraries, the plotting theme, the main
##     analysis's marginal_logLik(), morcela_weights() and morcela_score(), the
##     five models' colours, shapes and legend labels, and the main text's
##     DL-distance colours and buckets (DL_COLORS, DL_BUCKET_LEVELS,
##     bucket_dl(), YLIM) for the acceptability figures (setup.R also reads the
##     main experiments' data, which parts 5.2-5.4 do not use);
##   - sets the colours of the three conditions, for the correctness figures;
##   - builds the table of the 243 sentences that all three studies show.
##
## The sentences: 27 constructed (artificial) base items, each written in a
## short, medium and long version, and each version shown original, with 1 swap
## or with 3 swaps. That makes 9 conditions and 243 sentences, and every
## participant sees every base item once, in one of the 9 conditions. As in the
## controlled English experiment of the main text, every item comes in length
## versions of one base item, so the analyses group by the base item, with that
## experiment's by-item random effects. Surprisal then varies within item, which
## is what makes the by-item surprisal slopes identifiable.
## =============================================================================

## Libraries, plotting theme, helper functions, and the models' colours, shapes
## and labels
source("setup.R")
library(stringdist)

## ---- Colours ----
## The three conditions have DL distance 0, 1 and 3, and take the colours of
## those distances in the main text's figures
CONDITION_COLORS <- c("Original" = DL_COLORS[["0"]],
                      "1 swap"   = DL_COLORS[["1-2"]],
                      "3 swaps"  = DL_COLORS[["3"]])

## ---- The 243 sentences ----
## In each study's materials file, the sentence is on the row of the trial's
## first screen (a rating that follows another question has an empty sentence)
materials <- bind_rows(
  acceptability_ratings_only = read_csv(
    "../data/en_artificial_acceptability_ratings/materials_en_artificial_acceptability_ratings.csv",
    show_col_types = FALSE) %>%
    filter(experiment == "exp_b"),
  binary_error_detection_acceptability_ratings = read_csv(
    "../data/en_artificial_binary_error_detection/materials_en_artificial_binary_error_detection.csv",
    show_col_types = FALSE) %>%
    filter(experiment == "exp_b_comprehension"),
  exchange_identification_acceptability_ratings = read_csv(
    "../data/en_artificial_exchange_identification_acceptability/materials_en_artificial_exchange_identification_acceptability.csv",
    show_col_types = FALSE) %>%
    filter(experiment == "annot"),
  .id = "task"
) %>%
  transmute(task, item = as.integer(item), cond = condition, sentence)

## Every study shows each of the 243 sentences once, and all three show the
## same sentences
stopifnot(
  all(count(materials, task)$n == 27 * 9),
  nrow(distinct(materials, item, cond, sentence)) == 27 * 9
)

sentences <- materials %>%
  distinct(item, cond, sentence) %>%
  mutate(
    length_level = factor(str_extract(cond, "[^_]+$"),
                          levels = c("short", "medium", "long")),
    condition_label = factor(recode(str_extract(cond, "^[^_]+"),
                                    "original" = "Original",
                                    "swap-1" = "1 swap",
                                    "swap-3" = "3 swaps"),
                             levels = names(CONDITION_COLORS))
  )

## Damerau-Levenshtein distance from each version to its original, computed as
## in data/preprocess_stimuli.R: lowercase, drop punctuation, write each
## distinct word of the item (at that length) as one letter, and take
## stringdist(method = "dl") between the letter strings. Length is the number
## of space-separated words.
sentence_letters <- sentences %>%
  separate_rows(sentence, sep = " ") %>%
  mutate(word = str_remove_all(str_to_lower(sentence), "[^a-z0-9]")) %>%
  group_by(item, length_level) %>%
  mutate(letter = c(letters, LETTERS)[match(word, unique(word))]) %>%
  group_by(item, cond) %>%
  summarize(letter_string = paste(letter, collapse = ""), n_words = n(),
            .groups = "drop")

## The language-model scores, from the one source of the main analysis
## (data/slor_morcela_scores.csv: Goldfish 1000mb, log2, model-initial unigram)
artificial_scores <- read_csv("../data/slor_morcela_scores.csv",
                              col_types = cols(item = col_character()),
                              show_col_types = FALSE) %>%
  filter(language == "en_artificial", unigram_model == "model_initial") %>%
  transmute(item = as.integer(item), cond = condition, log_p_sentence, log_p_unigram)

## Information(s) in the communication account is the surprisal of the intended
## sentence, i.e. of the original at the same length: the original and swapped
## versions of an item at one length share it, and its short, medium and long
## versions differ. SLOR and MORCELA instead score the sentence that was
## actually read, so they take each condition's own log probabilities.
original_surprisal <- artificial_scores %>%
  filter(str_starts(cond, "original")) %>%
  transmute(item,
            length_level = factor(str_extract(cond, "[^_]+$"),
                                  levels = c("short", "medium", "long")),
            sum_surprisal = -log_p_sentence)

sentences <- sentences %>%
  left_join(sentence_letters, by = c("item", "cond")) %>%
  group_by(item, length_level) %>%
  mutate(original_string = letter_string[condition_label == "Original"]) %>%
  ungroup() %>%
  left_join(original_surprisal, by = c("item", "length_level")) %>%
  left_join(artificial_scores, by = c("item", "cond")) %>%
  mutate(
    dl_distance = stringdist(letter_string, original_string, method = "dl"),
    weighted_dl_distance = dl_distance / n_words,
    weighted_dl_distance_surprisal = dl_distance / sum_surprisal,
    grammatical_numeric = ifelse(dl_distance == 0, .5, -.5)
  ) %>%
  select(item, cond, sentence, length_level, condition_label, n_words,
         sum_surprisal, dl_distance, weighted_dl_distance,
         weighted_dl_distance_surprisal, grammatical_numeric,
         log_p_sentence, log_p_unigram)

## Every sentence found its scores, and every condition has the DL distance it
## was built to have: 0, 1 and 3. (So the main analysis's restriction to DL
## distances up to 3, which drops its shuffled condition, would remove nothing
## here, and is not repeated.)
stopifnot(
  nrow(sentences) == 27 * 9,
  !anyNA(sentences),
  all(sentences$dl_distance ==
        c("Original" = 0, "1 swap" = 1, "3 swaps" = 3)[as.character(sentences$condition_label)])
)

cat("DL distance by condition (number of sentences):\n")
sentences %>% count(condition_label, dl_distance) %>% print()

## The 27 original sentences of each length version: their length in words and
## their surprisal in bits (the Information(s) of all three conditions)
cat("Length (words) and surprisal (bits) of the original sentences, by length version:\n")
sentences %>%
  filter(condition_label == "Original") %>%
  group_by(length_level) %>%
  summarize(n_sentences = n(),
            min_words = min(n_words), median_words = median(n_words),
            max_words = max(n_words),
            min_surprisal = min(sum_surprisal), median_surprisal = median(sum_surprisal),
            max_surprisal = max(sum_surprisal)) %>%
  print(width = Inf)

## Everything defined so far is shared by 5.2-5.4 and the figure; each part
## below starts by removing what the previous part defined
shared_objects <- c(ls(), "shared_objects")


## #############################################################################
## 5.2  Acceptability ratings only (constructed sentences)
## #############################################################################
rm(list = setdiff(ls(), shared_objects))

## =============================================================================
## SI Section 5.2: Acceptability ratings only (constructed sentences)
##
## en_artificial_acceptability_ratings. The 243 constructed sentences (27 base
## items x short, medium and long x original, 1 swap and 3 swaps; see
## the shared set-up above), and every participant sees every base item once.
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
## bobyqa; singular fits are kept. ../analysis/main_analysis.R and sections
## 5.3-5.4 fit the linear mixed model of (2) with nlminb, but here
## nlminb stops short (see (2)). For every mixed model the CSVs record the
## optimizer, every message and warning printed while fitting it (lme4's and
## lmerTest's) and lme4's gradient and Hessian checks; see (2).
##
## Figure: figures/s5_acceptability_ratings_only_all_trials.png, the main
## text's human-judgment curves for this study
## Numbers: outputs/s5_acceptability_ratings_only_*.csv
## The log-likelihoods are drawn together with the other two studies' by
## the figure at the end of this script.
## =============================================================================

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

## Condition assignment: every participant rates each base item at most once,
## in one of its 9 conditions, and those who rated all 27 saw each condition 3
## times. The participants by how many items they rated, in how many of the 9
## conditions, and the fewest and most items they rated in one condition:
stopifnot(!anyDuplicated(trials[c("participant_id", "item")]))
trials %>%
  count(participant_id, cond) %>%
  group_by(participant_id) %>%
  summarize(items = sum(n), conditions = n(), fewest_per_condition = min(n),
            most_per_condition = max(n)) %>%
  count(items, conditions, fewest_per_condition, most_per_condition,
        name = "participants") %>%
  print(width = Inf)

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
print(sample_sizes)

## ---- (1) The five models on item x condition means ----
## As in the main analysis, the models are fit to the mean rating of each
## sentence (one base item at one length in one version) over the trials
## analysed. SLOR is computed as in the main analysis.
item_means <- analysis_trials %>%
  group_by(subset, item, cond, length_level, condition_label, n_words,
           sum_surprisal, dl_distance, weighted_dl_distance,
           weighted_dl_distance_surprisal, grammatical_numeric,
           log_p_sentence, log_p_unigram) %>%
  summarize(response_num = mean(acceptability), .groups = "drop") %>%
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
print(morcela_weights_by_subset)

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

## scam's inner fitting loop converged (one column per subset); the
## smoothing-parameter search is checked below
print(sapply(models_by_subset, sapply, `[[`, "conv"))

## Marginal log-likelihoods, with the random effects integrated out, as in the
## main analysis. `converged` is scam's inner fitting loop ($conv), which is
## TRUE even when scam's search for the smoothing parameters stopped short;
## bfgs_convergence is how that search ended ("Full convergence" when the
## gradient of its smoothness score, GCV/UBRE, became small enough).
model_fits <- tibble()
for (subset_iter in names(models_by_subset)) {
  model_fits <- bind_rows(model_fits, tibble(
    subset = subset_iter,
    fit = names(models_by_subset[[subset_iter]]),
    log_likelihood = vapply(models_by_subset[[subset_iter]], marginal_logLik, numeric(1)),
    converged = sapply(models_by_subset[[subset_iter]], `[[`, "conv"),
    bfgs_convergence = sapply(lapply(models_by_subset[[subset_iter]], `[[`, "bfgs.info"), `[[`, "conv"),
    optimizer = "bfgs"
  ))
}

model_fits %>% count(bfgs_convergence) %>% print()

## A fit whose search did not reach full convergence (it reached its 200-step
## limit, or its last step found no lower score) is refit with scam's two other
## optimizers, nlm and optim, as ../analysis/main_analysis.R checks its fits.
## The fit reported is the one with the lowest smoothness score among the
## original and the refits whose search converged; `optimizer` records which
## (bfgs is scam's default). nlm's two messages of success end in "probably
## solution"; optim's is "Successful completion".
for (row_index in which(model_fits$bfgs_convergence != "Full convergence")) {
  m <- models_by_subset[[model_fits$subset[row_index]]][[model_fits$fit[row_index]]]
  refit_fits <- list(
    bfgs = m,
    nlm = scam(m$formula, data = m$model, optimizer = "nlm"),
    optim = scam(m$formula, data = m$model, optimizer = "optim")
  )
  refits <- tibble(
    subset = model_fits$subset[row_index],
    fit = model_fits$fit[row_index],
    optimizer = names(refit_fits),
    search_convergence = c(m$bfgs.info$conv, refit_fits$nlm$nlm.info$conv,
                           refit_fits$optim$optim.info$conv),
    smoothness_score = sapply(refit_fits, `[[`, "gcv.ubre"),
    log_likelihood = vapply(refit_fits, marginal_logLik, numeric(1))
  )
  print(as.data.frame(refits), digits = 7)

  kept <- refits %>%
    filter(optimizer == "bfgs" |
             str_detect(search_convergence, "probably solution|Successful completion")) %>%
    slice_min(smoothness_score, n = 1, with_ties = FALSE)
  print(paste0("Kept: ", kept$optimizer, "; marginal log-likelihood change ",
               round(kept$log_likelihood - model_fits$log_likelihood[row_index], 3)))

  model_fits$optimizer[row_index] <- kept$optimizer
  model_fits$log_likelihood[row_index] <- kept$log_likelihood
  model_fits$converged[row_index] <- refit_fits[[kept$optimizer]]$conv
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
         likelihood_diff, converged, bfgs_convergence, optimizer,
         n_participants, n_trials, n_cells)

likelihoods %>%
  select(subset, information, model, log_likelihood, likelihood_diff, converged,
         optimizer) %>%
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

print(lapply(lmer_fits, summary))

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

print(lapply(error_free_fits, summary))

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
write_csv(likelihoods, paste0("outputs/s5_", task_name, "_likelihoods.csv"))
write_csv(lmm_dl_by_surprisal, paste0("outputs/s5_", task_name, "_lmm_dl_by_surprisal.csv"))
write_csv(error_free_surprisal_slopes,
          paste0("outputs/s5_", task_name, "_error_free_surprisal_slopes.csv"))
write_csv(morcela_weights_by_subset, paste0("outputs/s5_", task_name, "_morcela_weights.csv"))

## ---- Figure: acceptability ----
## The human-judgment curves of the main text's qualitative-predictions figure,
## with its plotting code (p_baseline and p_gam in analysis/main_analysis.R),
## for all trials of this study: it has no second task. Points are sentence
## means (one item at one length in one condition); lines are GAM smooths of
## those means over surprisal, one per DL distance, with 95% confidence bands.
## Like the main figure, it is drawn at its print size, one PNAS column (8.7 cm)
## wide, in the main figure's 8-pt print theme (theme_print_8pt there).
## Where it differs from the main code: it has no "Human Judgments" title (there
## are no model columns to tell it from) and shows its own legend, at the
## bottom like the main figure's; the legend is centred on the whole plot, as
## the main figure's collected legend is, so that it fits at this width; and it
## leaves out "4+" and labels the "1-2" bucket "1", since these sentences have a
## DL distance of 0, 1 or 3 only.
plot_data_by_item <- analysis_trials %>%
  filter(subset == "All trials") %>%
  group_by(subset, item, length_level, condition_label, sum_surprisal, dl_distance) %>%
  summarize(response = mean(acceptability), .groups = "drop") %>%
  mutate(dl_distance_bucket = bucket_dl(dl_distance))

p_acceptability <- ggplot(
  plot_data_by_item,
  aes(x = sum_surprisal, y = response, color = dl_distance_bucket)
) +
  geom_point(alpha = 0.5, size = 1.3, stroke = 0) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, bs = "tp"),
    se = TRUE,
    linewidth = 0.5,
    show.legend = FALSE
  ) +
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
    labels = c("0" = "0", "1-2" = "1", "3" = "3"),
    limits = DL_BUCKET_LEVELS,
    drop   = FALSE,
    guide  = guide_legend(
      override.aes = list(alpha = 1, linewidth = 0, size = 1.6)
    )
  )

if (interactive()) print(p_acceptability)

ggsave(paste0("figures/s5_", task_name, "_all_trials.png"),
       p_acceptability, width = 3.42, height = 2.6, units = "in", dpi = 300)


print("Done.")


## #############################################################################
## 5.3  Binary error detection + acceptability ratings (constructed sentences)
## #############################################################################
rm(list = setdiff(ls(), shared_objects))

## =============================================================================
## SI Section 5.3: Binary error detection + acceptability ratings (constructed
## sentences)
##
## en_artificial_binary_error_detection. The 243 constructed sentences (27 base
## items x short, medium and long x original, 1 swap and 3 swaps; see
## the shared set-up above), and every participant sees every base item once. On each
## trial the participant first answers "Are there any errors in this sentence?"
## (Yes/No) and then, on the next screen, rates how acceptable the sentence is
## (1-7).
##
## Judged correctly: the answer to the binary error detection question matches
## literal_response, i.e. "No" (no error) for an original and "Yes" for a
## swapped sentence.
##
## The main text's analyses, with sentence surprisal as Information(s) and the
## specifications of the main analysis's controlled English experiment, each on
## all trials and again on the trials judged correctly only:
##   (1) the five models (Binary, Gradient, Communication, SLOR, MORCELA) as
##       shape-constrained GAMs on item x condition means, compared on their
##       marginal log-likelihood;
##   (2) the linear mixed model of acceptability on DL distance x surprisal,
##       over trials;
##   (3) the surprisal slope among the original (error-free) sentences;
##   (4) the length check, in the CSV only: Communication with length in place
##       of surprisal.
## Restricting to the trials judged correctly selects on a response that is
## correlated with the rating, so all trials are the baseline to read the
## correct-trial results against.
##
## Mixed models have the maximal random-effects structure. As in
## ../analysis/main_analysis.R, the linear mixed models of (2) are fit with
## nlminb and the error-free models with bobyqa (see (2)). Singular fits are
## kept. For every mixed model the CSVs record the optimizer, every message
## and warning printed while fitting it (lme4's and lmerTest's) and lme4's
## gradient and Hessian checks.
##
## Figures: figures/s5_binary_error_detection_acceptability_ratings_yes_responses.png
## (proportion of "yes" answers, over all binary error detection judgments) and
## figures/s5_binary_error_detection_acceptability_ratings_correct_trials.png
## (the main text's human-judgment curves, on the trials judged correctly)
## Numbers: outputs/s5_binary_error_detection_acceptability_ratings_*.csv
## The log-likelihoods are drawn together with the other two studies' by
## the figure at the end of this script.
## =============================================================================

task_name <- "binary_error_detection_acceptability_ratings"

## ---- Responses: one row per screen: the binary error detection question (exp_b_comprehension), then the rating (exp_b) ----
responses <- read_csv("../data/en_artificial_binary_error_detection/en_artificial_binary_error_detection_data.csv",
                      show_col_types = FALSE)

cat("\n---- export integrity ----\n")
cat("participants in the file:", n_distinct(responses$participant_id), "\n")
responses %>% distinct(participant_id, native_language) %>% count(native_language) %>% print()

## The study is of native English speakers: drop anyone who answered No (nobody
## did here). A blank answer is kept, as in the other two studies.
non_native <- responses %>%
  distinct(participant_id, native_language) %>%
  filter(native_language %in% "No") %>%
  pull(participant_id)
cat("excluded as non-native:", length(non_native), "participant(s)\n")
responses <- responses %>% filter(!participant_id %in% non_native)

## Participants who also took the exchange identification + acceptability
## ratings study (section 5.4) saw the same sentences twice; they are excluded
## from both studies. The studies of constructed sentences share one set of
## pseudonyms (en_artificial_<n>), so the same id in both is the same person.
took_both <- intersect(
  responses$participant_id,
  read_csv("../data/en_artificial_exchange_identification_acceptability/en_artificial_exchange_identification_acceptability_data.csv",
           show_col_types = FALSE)$participant_id)
cat("excluded for also taking the exchange identification + acceptability ratings study:",
    length(took_both), "participant(s)\n")
responses <- responses %>% filter(!participant_id %in% took_both)

judgments <- responses %>%
  filter(experiment == "exp_b_comprehension") %>%
  transmute(participant_id, item = as.integer(item), cond,
            judgment_order = presentation_order,
            ## Judged correctly: "No" for an original, "Yes" for a swap
            correct = response == literal_response,
            ## Reported an error: answered "Yes"
            reported_error = response == "Yes")
## A "yes" is correct for every sentence except an original
stopifnot(all(judgments$reported_error == (judgments$correct != str_detect(judgments$cond, "^original"))))

ratings <- responses %>%
  filter(experiment == "exp_b") %>%
  transmute(participant_id, item = as.integer(item), cond,
            rating_order = presentation_order,
            acceptability = as.numeric(substr(response, 1, 1)))

responses %>% count(participant_id, name = "rows") %>% count(rows, name = "participants") %>% print()
cat("(54 rows = 27 sentences x 2 questions; Surveyor writes no row for a skipped question)\n")
cat("judgments without a rating:",
    nrow(anti_join(judgments, ratings, by = c("participant_id", "item", "cond"))),
    "| ratings without a judgment:",
    nrow(anti_join(ratings, judgments, by = c("participant_id", "item", "cond"))), "\n")

## Condition assignment: every participant judges and rates each base item at
## most once, in one of its 9 conditions, and those who answered all 27 binary
## error detection questions saw each condition 3 times. The participants by
## how many items they judged, in how many of the 9 conditions, and the fewest
## and most items they judged in one condition:
stopifnot(
  !anyDuplicated(judgments[c("participant_id", "item")]),
  !anyDuplicated(ratings[c("participant_id", "item")])
)
judgments %>%
  count(participant_id, cond) %>%
  group_by(participant_id) %>%
  summarize(items = sum(n), conditions = n(), fewest_per_condition = min(n),
            most_per_condition = max(n)) %>%
  count(items, conditions, fewest_per_condition, most_per_condition,
        name = "participants") %>%
  print(width = Inf)

## Every binary error detection judgment, rated or not, with its sentence: the
## correctness table and figure use all of them, as section 5.4 uses all
## exchange identification trials
judgment_trials <- judgments %>%
  left_join(sentences, by = c("item", "cond"))

## A trial is a sentence with both its binary error detection judgment and its
## rating: the acceptability analyses use these
trials <- judgments %>%
  inner_join(ratings, by = c("participant_id", "item", "cond")) %>%
  left_join(sentences, by = c("item", "cond"))

## The rating is always the screen right after its binary error detection
## question
stopifnot(
  !anyNA(judgment_trials),
  !anyNA(trials),
  all(trials$rating_order == trials$judgment_order + 1)
)
cat("binary error detection judgments:", nrow(judgment_trials),
    "| participants analysed:", n_distinct(trials$participant_id),
    "| trials (judged and rated):", nrow(trials), "\n")

cells <- trials %>% count(item, cond)
cat("item x condition cells seen:", nrow(cells), "of", nrow(sentences),
    "| responses per cell: median", median(cells$n),
    "range", min(cells$n), "-", max(cells$n), "\n")

## ---- Binary error detection: proportion judged correctly ----
## Over all binary error detection judgments, rated or not (as section 5.4
## uses all exchange identification trials), by condition and length, and over
## all lengths
correctness <- bind_rows(
  judgment_trials %>%
    group_by(condition_label, length_level = as.character(length_level)) %>%
    summarize(n_trials = n(), proportion_correct = mean(correct),
              proportion_reported_error = mean(reported_error), .groups = "drop"),
  judgment_trials %>%
    group_by(condition_label) %>%
    summarize(n_trials = n(), proportion_correct = mean(correct),
              proportion_reported_error = mean(reported_error), .groups = "drop") %>%
    mutate(length_level = "all")
) %>%
  mutate(task = task_name, .before = 1) %>%
  arrange(condition_label, factor(length_level, levels = c("short", "medium", "long", "all")))

print(correctness, n = Inf)
write_csv(correctness, paste0("outputs/s5_", task_name, "_correctness.csv"))

## ---- The trials analysed ----
## One long table, with every analysis below run on each subset of it: all
## trials, and the trials judged correctly
analysis_trials <- bind_rows(
  "All trials" = trials,
  "Judged correctly" = trials %>% filter(correct),
  .id = "subset"
)

sample_sizes <- analysis_trials %>%
  group_by(subset) %>%
  summarize(n_participants = n_distinct(participant_id), n_trials = n())
print(sample_sizes)

cat("\n---- mean rating by condition and length ----\n")
analysis_trials %>%
  group_by(subset, condition_label, length_level) %>%
  summarize(mean_rating = round(mean(acceptability), 2), .groups = "drop") %>%
  pivot_wider(names_from = length_level, values_from = mean_rating) %>%
  print()

## ---- (1) The five models on item x condition means ----
## As in the main analysis, the models are fit to the mean rating of each
## sentence (one base item at one length in one version) over the trials
## analysed. SLOR is computed as in the main analysis.
item_means <- analysis_trials %>%
  group_by(subset, item, cond, length_level, condition_label, n_words,
           sum_surprisal, dl_distance, weighted_dl_distance,
           weighted_dl_distance_surprisal, grammatical_numeric,
           log_p_sentence, log_p_unigram) %>%
  summarize(response_num = mean(acceptability), .groups = "drop") %>%
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
print(morcela_weights_by_subset)

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

## scam's inner fitting loop converged (one column per subset); the
## smoothing-parameter search is checked below
print(sapply(models_by_subset, sapply, `[[`, "conv"))

## Marginal log-likelihoods, with the random effects integrated out, as in the
## main analysis. `converged` is scam's inner fitting loop ($conv), which is
## TRUE even when scam's search for the smoothing parameters stopped short;
## bfgs_convergence is how that search ended ("Full convergence" when the
## gradient of its smoothness score, GCV/UBRE, became small enough).
model_fits <- tibble()
for (subset_iter in names(models_by_subset)) {
  model_fits <- bind_rows(model_fits, tibble(
    subset = subset_iter,
    fit = names(models_by_subset[[subset_iter]]),
    log_likelihood = vapply(models_by_subset[[subset_iter]], marginal_logLik, numeric(1)),
    converged = sapply(models_by_subset[[subset_iter]], `[[`, "conv"),
    bfgs_convergence = sapply(lapply(models_by_subset[[subset_iter]], `[[`, "bfgs.info"), `[[`, "conv"),
    optimizer = "bfgs"
  ))
}

model_fits %>% count(bfgs_convergence) %>% print()

## A fit whose search did not reach full convergence (it reached its 200-step
## limit, or its last step found no lower score) is refit with scam's two other
## optimizers, nlm and optim, as ../analysis/main_analysis.R checks its fits.
## The fit reported is the one with the lowest smoothness score among the
## original and the refits whose search converged; `optimizer` records which
## (bfgs is scam's default). nlm's two messages of success end in "probably
## solution"; optim's is "Successful completion".
for (row_index in which(model_fits$bfgs_convergence != "Full convergence")) {
  m <- models_by_subset[[model_fits$subset[row_index]]][[model_fits$fit[row_index]]]
  refit_fits <- list(
    bfgs = m,
    nlm = scam(m$formula, data = m$model, optimizer = "nlm"),
    optim = scam(m$formula, data = m$model, optimizer = "optim")
  )
  refits <- tibble(
    subset = model_fits$subset[row_index],
    fit = model_fits$fit[row_index],
    optimizer = names(refit_fits),
    search_convergence = c(m$bfgs.info$conv, refit_fits$nlm$nlm.info$conv,
                           refit_fits$optim$optim.info$conv),
    smoothness_score = sapply(refit_fits, `[[`, "gcv.ubre"),
    log_likelihood = vapply(refit_fits, marginal_logLik, numeric(1))
  )
  print(as.data.frame(refits), digits = 7)

  kept <- refits %>%
    filter(optimizer == "bfgs" |
             str_detect(search_convergence, "probably solution|Successful completion")) %>%
    slice_min(smoothness_score, n = 1, with_ties = FALSE)
  print(paste0("Kept: ", kept$optimizer, "; marginal log-likelihood change ",
               round(kept$log_likelihood - model_fits$log_likelihood[row_index], 3)))

  model_fits$optimizer[row_index] <- kept$optimizer
  model_fits$log_likelihood[row_index] <- kept$log_likelihood
  model_fits$converged[row_index] <- refit_fits[[kept$optimizer]]$conv
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
         likelihood_diff, converged, bfgs_convergence, optimizer,
         n_participants, n_trials, n_cells)

likelihoods %>%
  select(subset, information, model, log_likelihood, likelihood_diff, converged,
         optimizer) %>%
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

## Both subsets: nlminb (nlminbwrap). With this maximal structure bobyqa stops
## short of the optimum on all trials: lmerTest warns "Model failed to converge
## with 1 negative eigenvalue: -1.4e+02", and lme4's own Hessian has the same
## negative eigenvalue. In lme4's allFit(), nlminbwrap and NLopt's BOBYQA reach
## the same, slightly better optimum (REML criterion 9925.716, against bobyqa's
## 9926.004), with fixed effects identical to 3 decimals, and the nlminbwrap fit
## passes every check. Its fixed effects differ from bobyqa's by at most
## 0.00131 (the surprisal slope). The random effects are the same maximal ones.
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
      control = lmerControl(
        optimizer ="nlminbwrap",
        optCtrl = list(eval.max = 1e5, iter.max = 1e5)
      )
    ),
    type = "message"
  )
}

print(lapply(lmer_fits, summary))

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

print(lapply(error_free_fits, summary))

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
write_csv(likelihoods, paste0("outputs/s5_", task_name, "_likelihoods.csv"))
write_csv(lmm_dl_by_surprisal, paste0("outputs/s5_", task_name, "_lmm_dl_by_surprisal.csv"))
write_csv(error_free_surprisal_slopes,
          paste0("outputs/s5_", task_name, "_error_free_surprisal_slopes.csv"))
write_csv(morcela_weights_by_subset, paste0("outputs/s5_", task_name, "_morcela_weights.csv"))


## ---- Figure: "yes" responses in binary error detection ----
## Over all binary error detection judgments, rated or not, like the
## correctness table. The y axis is the proportion of "yes" answers to "Are
## there any errors in this sentence?": false alarms for the originals, and
## detected errors for the swapped sentences. Points are sentences (one item at one length in one condition). The sentences of one item at one length share
## their x (the surprisal of the intended sentence), so the points are
## jittered sideways, by up to 2 bits, and never up or down, so that no point
## is drawn outside 0-100%. Lines are logistic regressions of the trial-level
## outcome on surprisal, one per condition; their 95% confidence bands are
## computed on the logit scale, so they stay within 0-100%. They ignore that
## trials are grouped by participant and item. The legend lists the conditions
## from top to bottom in the order of the curves, and calls the swaps
## "exchanges", as the text does.
p_yes_responses <- ggplot(mapping = aes(x = sum_surprisal, color = condition_label,
                                        fill = condition_label)) +
  geom_point(data = judgment_trials %>%
               group_by(condition_label, item, length_level, sum_surprisal) %>%
               summarize(value = mean(reported_error), .groups = "drop"),
             aes(y = value), alpha = 0.4,
             position = position_jitter(width = 2, height = 0, seed = 2026)) +
  geom_smooth(data = judgment_trials, aes(y = as.integer(reported_error)),
              method = "glm", method.args = list(family = binomial), formula = y ~ x) +
  scale_color_manual(name = "Condition", values = CONDITION_COLORS,
                     labels = c("Original" = "Original", "1 swap" = "1 exchange",
                                "3 swaps" = "3 exchanges"),
                     aesthetics = c("colour", "fill"),
                     guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous(name = "Sentence surprisal (bits)") +
  scale_y_continuous(name = "Proportion of \"yes\" responses\n(reported an error)",
                     labels = scales::percent, breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1))

if (interactive()) print(p_yes_responses)

ggsave(paste0("figures/s5_", task_name, "_yes_responses.png"),
       p_yes_responses, width = 8, height = 5.5)


## ---- Figure: acceptability on the trials judged correctly ----
## The human-judgment curves of the main text's qualitative-predictions figure,
## with its plotting code (p_baseline and p_gam in analysis/main_analysis.R),
## for the trials judged correctly. Points are sentence means (one item at one
## length in one condition); lines are GAM smooths of those means over
## surprisal, one per DL distance, with 95% confidence bands. Like the main
## figure, it is drawn at its print size, one PNAS column (8.7 cm) wide, in the
## main figure's 8-pt print theme (theme_print_8pt there). Where it differs
## from the main code: it has no "Human Judgments" title (there are no model
## columns to tell it from) and shows its own legend, at the bottom like the
## main figure's; the legend is centred on the whole plot, as the main figure's
## collected legend is, so that it fits at this width; and it leaves out "4+"
## and labels the "1-2" bucket "1", since these sentences have a DL distance of
## 0, 1 or 3 only.
plot_data_by_item <- analysis_trials %>%
  filter(subset == "Judged correctly") %>%
  group_by(subset, item, length_level, condition_label, sum_surprisal, dl_distance) %>%
  summarize(response = mean(acceptability), .groups = "drop") %>%
  mutate(dl_distance_bucket = bucket_dl(dl_distance))

p_acceptability <- ggplot(
  plot_data_by_item,
  aes(x = sum_surprisal, y = response, color = dl_distance_bucket)
) +
  geom_point(alpha = 0.5, size = 1.3, stroke = 0) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, bs = "tp"),
    se = TRUE,
    linewidth = 0.5,
    show.legend = FALSE
  ) +
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
    labels = c("0" = "0", "1-2" = "1", "3" = "3"),
    limits = DL_BUCKET_LEVELS,
    drop   = FALSE,
    guide  = guide_legend(
      override.aes = list(alpha = 1, linewidth = 0, size = 1.6)
    )
  )

if (interactive()) print(p_acceptability)

ggsave(paste0("figures/s5_", task_name, "_correct_trials.png"),
       p_acceptability, width = 3.42, height = 2.6, units = "in", dpi = 300)


print("Done.")


## #############################################################################
## 5.4  Exchange identification + acceptability ratings (constructed sentences)
## #############################################################################
rm(list = setdiff(ls(), shared_objects))

## =============================================================================
## SI Section 5.4: Exchange identification + acceptability ratings
## (constructed sentences)
##
## en_artificial_exchange_identification_acceptability. The 243 constructed
## sentences (27 base items x short, medium and long x original, 1 swap and 3
## swaps; see the shared set-up above), and every participant sees every base item
## once. On each trial the participant first clicks the words they think are
## misplaced (or nothing) and then, on the next screen, rates how acceptable the
## sentence is (1-7).
##
## Identified correctly -- the rule called click_correct in the earlier section
## 4.5 script (git history, commit 2537ec2):
##   - on an original: no word clicked;
##   - on a swapped sentence, with 1 swap or 3: at least one swapped region
##     marked (at least one of its displaced words clicked), and no click away
##     from the swapped regions. Every clicked word must be displaced or be
##     joined to a swapped region by an unbroken run of clicked words, so a
##     click on a word next to a swapped region still counts, but a click on a
##     word cut off from every swapped region by an unclicked, undisplaced word
##     makes the trial incorrect.
## A displaced word is one that differs from the word at the same position in
## the original. Consecutive displaced words form one region, so swaps that
## touch merge into one region. On a 3-swap sentence, then, marking any one
## swapped region is enough: neither every region nor every swap has to be
## marked. (Of the 81 3-swap sentences, 20 have one region, 36 two and 25
## three.) For the stricter reading, the correctness CSV also gives the share of
## trials identified correctly with every swapped region marked
## (proportion_all_regions_marked); on originals and 1-swap sentences, which
## have at most one region, it equals proportion_correct.
##
## The main text's analyses, with sentence surprisal as Information(s) and the
## specifications of the main analysis's controlled English experiment, each on
## all rated trials and again on the rated trials identified correctly only:
##   (1) the five models (Binary, Gradient, Communication, SLOR, MORCELA) as
##       shape-constrained GAMs on item x condition means, compared on their
##       marginal log-likelihood;
##   (2) the linear mixed model of acceptability on DL distance x surprisal,
##       over trials;
##   (3) the surprisal slope among the original (error-free) sentences;
##   (4) the length check, in the CSV only: Communication with length in place
##       of surprisal.
## Restricting to the trials identified correctly selects on a response that is
## correlated with the rating, so all trials are the baseline to read the
## correct-trial results against.
##
## Mixed models have the maximal random-effects structure. As in
## ../analysis/main_analysis.R, the linear mixed models of (2) are fit with
## nlminb and the error-free models with bobyqa (see (2)). Singular fits are kept. For every mixed model the CSVs record the
## optimizer, every message and warning printed while fitting it (lme4's and
## lmerTest's) and lme4's gradient and Hessian checks.
##
## Figures: figures/s5_exchange_identification_acceptability_ratings_correctness.png
## (proportion identified correctly, over all exchange identification trials) and
## figures/s5_exchange_identification_acceptability_ratings_correct_trials.png
## (the main text's human-judgment curves, on the rated trials identified
## correctly)
## Numbers: outputs/s5_exchange_identification_acceptability_ratings_*.csv
## The log-likelihoods are drawn together with the other two studies' by
## the figure at the end of this script.
## =============================================================================

task_name <- "exchange_identification_acceptability_ratings"

## ---- Words of each sentence ----
## One row per word, with its position. A word is displaced (a target) when it
## differs from the word at the same position in the item's original of the
## same length; consecutive targets form a region (one per swapped pair, or one
## per group of swaps that touch).
words <- sentences %>%
  select(item, cond, length_level, condition_label, sentence) %>%
  separate_rows(sentence, sep = " ") %>%
  rename(word = sentence) %>%
  group_by(item, cond) %>%
  mutate(position = row_number()) %>%
  ungroup()

words <- words %>%
  left_join(words %>%
              filter(condition_label == "Original") %>%
              select(item, length_level, position, original_word = word),
            by = c("item", "length_level", "position")) %>%
  mutate(is_target = word != original_word) %>%
  group_by(item, cond) %>%
  mutate(target_run = ifelse(is_target,
                             cumsum(is_target & !lag(is_target, default = FALSE)),
                             NA)) %>%
  ungroup()

stopifnot(!anyNA(words$original_word))

cat("\nswapped regions per sentence (number of sentences):\n")
words %>%
  group_by(condition_label, item, cond) %>%
  summarize(regions = n_distinct(target_run, na.rm = TRUE), .groups = "drop") %>%
  count(condition_label, regions) %>%
  print()

## ---- Responses ----
responses <- read_csv("../data/en_artificial_exchange_identification_acceptability/en_artificial_exchange_identification_acceptability_data.csv",
                      show_col_types = FALSE)

cat("\n---- export integrity ----\n")
cat("participants in the file:", n_distinct(responses$participant_id), "\n")
responses %>% count(participant_id, experiment) %>%
  count(experiment, n, name = "participants") %>% print()
cat("(27 exchange identification rows each = complete: that screen always writes a",
    "row; a skipped rating writes none)\n")
responses %>% distinct(participant_id, native_language) %>% count(native_language) %>% print()

## The study is of native English speakers: drop anyone who answered No (nobody
## did here). A blank answer is kept, as in the other two studies.
non_native <- responses %>%
  distinct(participant_id, native_language) %>%
  filter(native_language %in% "No") %>%
  pull(participant_id)
cat("excluded as non-native:", length(non_native), "participant(s)\n")
responses <- responses %>% filter(!participant_id %in% non_native)

## Participants who also took the binary error detection + acceptability
## ratings study (section 5.3) saw the same sentences twice; they are excluded
## from both studies. The studies of constructed sentences share one set of
## pseudonyms (en_artificial_<n>), so the same id in both is the same person.
took_both <- intersect(
  responses$participant_id,
  read_csv("../data/en_artificial_binary_error_detection/en_artificial_binary_error_detection_data.csv",
           show_col_types = FALSE)$participant_id)
cat("excluded for also taking the binary error detection + acceptability ratings study:",
    length(took_both), "participant(s)\n")
responses <- responses %>% filter(!participant_id %in% took_both)

ratings <- responses %>%
  filter(experiment == "exp_b") %>%
  transmute(participant_id, item = as.integer(item), cond,
            rating_order = presentation_order,
            acceptability = as.numeric(substr(response, 1, 1)))

## One row per exchange identification trial, with its rating when there is one
identification_trials <- responses %>%
  filter(experiment == "annot") %>%
  mutate(trial = row_number()) %>%
  transmute(trial, participant_id, item = as.integer(item), cond,
            identification_order = presentation_order, response) %>%
  left_join(ratings, by = c("participant_id", "item", "cond"))

cat("exchange identification trials:", nrow(identification_trials),
    "| of which rated:", sum(!is.na(identification_trials$acceptability)),
    "| ratings without an exchange identification trial:",
    nrow(anti_join(ratings, identification_trials, by = c("participant_id", "item", "cond"))), "\n")

## Condition assignment: every participant identifies and rates each base item
## at most once, in one of its 9 conditions, and the rating is always the screen
## right after its exchange identification screen. The participants by how many
## items they saw, in how many of the 9 conditions, and the fewest and most
## items they saw in one condition:
stopifnot(
  !anyDuplicated(identification_trials[c("participant_id", "item")]),
  !anyDuplicated(ratings[c("participant_id", "item")]),
  all(is.na(identification_trials$rating_order) == is.na(identification_trials$acceptability)),
  all(identification_trials$rating_order == identification_trials$identification_order + 1,
      na.rm = TRUE)
)
identification_trials %>%
  count(participant_id, cond) %>%
  group_by(participant_id) %>%
  summarize(items = sum(n), conditions = n(), fewest_per_condition = min(n),
            most_per_condition = max(n)) %>%
  count(items, conditions, fewest_per_condition, most_per_condition,
        name = "participants") %>%
  print(width = Inf)

## ---- Scoring each exchange identification trial ----
## One row per clicked word. The response lists the clicked words as
## position:word pairs, with the words percent-escaped for % , and : so that
## they survive this encoding, hence URLdecode.
clicks <- identification_trials %>%
  filter(response != "NONE") %>%
  separate_rows(response, sep = ",") %>%
  separate(response, into = c("position", "clicked_word"), sep = ":", extra = "merge") %>%
  transmute(trial, item, cond, position = as.integer(position),
            clicked_word = map_chr(clicked_word, URLdecode))

## Every click is on the word the participant saw at that position
click_check <- clicks %>%
  left_join(words %>% select(item, cond, position, word), by = c("item", "cond", "position")) %>%
  summarize(clicks = n(), matching_the_stimulus = sum(clicked_word == word, na.rm = TRUE),
            out_of_range = sum(is.na(word)))
print(click_check)
stopifnot(
  click_check$matching_the_stimulus == click_check$clicks,
  click_check$out_of_range == 0
)

## Every word of the sentence each participant saw, marked clicked or not
trial_words <- identification_trials %>%
  select(trial, item, cond) %>%
  inner_join(words %>% select(item, cond, position, is_target, target_run),
             by = c("item", "cond"), relationship = "many-to-many") %>%
  left_join(clicks %>% transmute(trial, position, clicked = TRUE),
            by = c("trial", "position")) %>%
  mutate(clicked = replace_na(clicked, FALSE))

## A stray click: split the clicked-or-target positions of a trial into runs of
## consecutive positions, and flag any run with no target in it. A new run
## starts wherever a position is not the previous position + 1 (this relies on
## + binding tighter than !=)
stray_clicks <- trial_words %>%
  filter(clicked | is_target) %>%
  arrange(trial, position) %>%
  group_by(trial) %>%
  mutate(run = cumsum(position != lag(position, default = -1L) + 1L)) %>%
  group_by(trial, run) %>%
  summarize(stray = !any(is_target), .groups = "drop") %>%
  group_by(trial) %>%
  summarize(any_stray = any(stray), .groups = "drop")

identification_trials <- identification_trials %>%
  left_join(trial_words %>%
              group_by(trial) %>%
              summarize(n_clicked = sum(clicked),
                        n_on_target = sum(clicked & is_target),
                        regions_total = n_distinct(target_run, na.rm = TRUE),
                        regions_marked = n_distinct(target_run[clicked & is_target], na.rm = TRUE),
                        .groups = "drop"),
            by = "trial") %>%
  left_join(stray_clicks, by = "trial") %>%
  left_join(sentences %>% select(-sentence), by = c("item", "cond")) %>%
  mutate(
    flagged = n_clicked > 0,
    ## any_stray is NA when nothing was clicked on an original (there is no
    ## run to check); located is FALSE then, since no target was clicked
    located = n_on_target > 0 & !any_stray,
    ## Identified correctly: nothing clicked on an original, located on a
    ## swapped sentence
    correct = ifelse(condition_label == "Original", !flagged, located),
    ## The stricter reading, reported alongside: correct, with every swapped
    ## region marked (the same as correct on originals and 1-swap sentences)
    all_regions_marked = correct & regions_marked == regions_total
  )

stopifnot(!anyNA(identification_trials %>% select(-acceptability, -rating_order, -any_stray)))

## ---- Exchange identification: proportion identified correctly ----
## Over all exchange identification trials, rated or not, by condition and
## length, and over all lengths. proportion_all_regions_marked is the stricter
## reading: on 3-swap sentences, the share of trials identified correctly that
## also marked every swapped region.
correctness <- bind_rows(
  identification_trials %>%
    group_by(condition_label, length_level = as.character(length_level)) %>%
    summarize(n_trials = n(), proportion_correct = mean(correct),
              proportion_all_regions_marked = mean(all_regions_marked), .groups = "drop"),
  identification_trials %>%
    group_by(condition_label) %>%
    summarize(n_trials = n(), proportion_correct = mean(correct),
              proportion_all_regions_marked = mean(all_regions_marked), .groups = "drop") %>%
    mutate(length_level = "all")
) %>%
  mutate(task = task_name, .before = 1) %>%
  arrange(condition_label, factor(length_level, levels = c("short", "medium", "long", "all")))

print(correctness, n = Inf, width = Inf)
write_csv(correctness, paste0("outputs/s5_", task_name, "_correctness.csv"))

## ---- The trials analysed ----
## The acceptability analyses use the exchange identification trials that have
## a rating. One long table, with every analysis below run on each subset of
## it: all rated trials, and the rated trials identified correctly.
trials <- identification_trials %>%
  filter(!is.na(acceptability))

cat("participants analysed:", n_distinct(trials$participant_id), "| rated trials:", nrow(trials), "\n")

cells <- trials %>% count(item, cond)
cat("item x condition cells seen:", nrow(cells), "of", nrow(sentences),
    "| ratings per cell: median", median(cells$n),
    "range", min(cells$n), "-", max(cells$n), "\n")

analysis_trials <- bind_rows(
  "All trials" = trials,
  "Identified correctly" = trials %>% filter(correct),
  .id = "subset"
)

sample_sizes <- analysis_trials %>%
  group_by(subset) %>%
  summarize(n_participants = n_distinct(participant_id), n_trials = n())
print(sample_sizes)

cat("\n---- mean rating by condition and length ----\n")
analysis_trials %>%
  group_by(subset, condition_label, length_level) %>%
  summarize(mean_rating = round(mean(acceptability), 2), .groups = "drop") %>%
  pivot_wider(names_from = length_level, values_from = mean_rating) %>%
  print()

## ---- (1) The five models on item x condition means ----
## As in the main analysis, the models are fit to the mean rating of each
## sentence (one base item at one length in one version) over the trials
## analysed. SLOR is computed as in the main analysis.
item_means <- analysis_trials %>%
  group_by(subset, item, cond, length_level, condition_label, n_words,
           sum_surprisal, dl_distance, weighted_dl_distance,
           weighted_dl_distance_surprisal, grammatical_numeric,
           log_p_sentence, log_p_unigram) %>%
  summarize(response_num = mean(acceptability), .groups = "drop") %>%
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
print(morcela_weights_by_subset)

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

## scam's inner fitting loop converged (one column per subset); the
## smoothing-parameter search is checked below
print(sapply(models_by_subset, sapply, `[[`, "conv"))

## Marginal log-likelihoods, with the random effects integrated out, as in the
## main analysis. `converged` is scam's inner fitting loop ($conv), which is
## TRUE even when scam's search for the smoothing parameters stopped short;
## bfgs_convergence is how that search ended ("Full convergence" when the
## gradient of its smoothness score, GCV/UBRE, became small enough).
model_fits <- tibble()
for (subset_iter in names(models_by_subset)) {
  model_fits <- bind_rows(model_fits, tibble(
    subset = subset_iter,
    fit = names(models_by_subset[[subset_iter]]),
    log_likelihood = vapply(models_by_subset[[subset_iter]], marginal_logLik, numeric(1)),
    converged = sapply(models_by_subset[[subset_iter]], `[[`, "conv"),
    bfgs_convergence = sapply(lapply(models_by_subset[[subset_iter]], `[[`, "bfgs.info"), `[[`, "conv"),
    optimizer = "bfgs"
  ))
}

model_fits %>% count(bfgs_convergence) %>% print()

## A fit whose search did not reach full convergence (it reached its 200-step
## limit, or its last step found no lower score) is refit with scam's two other
## optimizers, nlm and optim, as ../analysis/main_analysis.R checks its fits.
## The fit reported is the one with the lowest smoothness score among the
## original and the refits whose search converged; `optimizer` records which
## (bfgs is scam's default). nlm's two messages of success end in "probably
## solution"; optim's is "Successful completion".
for (row_index in which(model_fits$bfgs_convergence != "Full convergence")) {
  m <- models_by_subset[[model_fits$subset[row_index]]][[model_fits$fit[row_index]]]
  refit_fits <- list(
    bfgs = m,
    nlm = scam(m$formula, data = m$model, optimizer = "nlm"),
    optim = scam(m$formula, data = m$model, optimizer = "optim")
  )
  refits <- tibble(
    subset = model_fits$subset[row_index],
    fit = model_fits$fit[row_index],
    optimizer = names(refit_fits),
    search_convergence = c(m$bfgs.info$conv, refit_fits$nlm$nlm.info$conv,
                           refit_fits$optim$optim.info$conv),
    smoothness_score = sapply(refit_fits, `[[`, "gcv.ubre"),
    log_likelihood = vapply(refit_fits, marginal_logLik, numeric(1))
  )
  print(as.data.frame(refits), digits = 7)

  kept <- refits %>%
    filter(optimizer == "bfgs" |
             str_detect(search_convergence, "probably solution|Successful completion")) %>%
    slice_min(smoothness_score, n = 1, with_ties = FALSE)
  print(paste0("Kept: ", kept$optimizer, "; marginal log-likelihood change ",
               round(kept$log_likelihood - model_fits$log_likelihood[row_index], 3)))

  model_fits$optimizer[row_index] <- kept$optimizer
  model_fits$log_likelihood[row_index] <- kept$log_likelihood
  model_fits$converged[row_index] <- refit_fits[[kept$optimizer]]$conv
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
         likelihood_diff, converged, bfgs_convergence, optimizer,
         n_participants, n_trials, n_cells)

likelihoods %>%
  select(subset, information, model, log_likelihood, likelihood_diff, converged,
         optimizer) %>%
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

## Both subsets: nlminb (nlminbwrap). With this maximal structure bobyqa stops
## short of the optimum on the trials identified correctly: lmerTest warns
## "Model failed to converge with 1 negative eigenvalue: -9.7e+01", and lme4's
## own Hessian has the same negative eigenvalue. In lme4's allFit(), nlminbwrap
## and NLopt's BOBYQA reach the same, slightly better optimum (REML criterion
## 8633.021, against bobyqa's 8633.331), with fixed effects identical to 3
## decimals, and the nlminbwrap fit passes every check. Its fixed effects differ
## from bobyqa's by at most 0.00103 (the intercept). The random effects are the
## same maximal ones.
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
      control = lmerControl(
        optimizer ="nlminbwrap",
        optCtrl = list(eval.max = 1e5, iter.max = 1e5)
      )
    ),
    type = "message"
  )
}

print(lapply(lmer_fits, summary))

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

print(lapply(error_free_fits, summary))

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
write_csv(likelihoods, paste0("outputs/s5_", task_name, "_likelihoods.csv"))
write_csv(lmm_dl_by_surprisal, paste0("outputs/s5_", task_name, "_lmm_dl_by_surprisal.csv"))
write_csv(error_free_surprisal_slopes,
          paste0("outputs/s5_", task_name, "_error_free_surprisal_slopes.csv"))
write_csv(morcela_weights_by_subset, paste0("outputs/s5_", task_name, "_morcela_weights.csv"))


## ---- Figure: correctness in exchange identification ----
## Over all exchange identification trials, rated or not, like the correctness
## table. Originals are correct when nothing was clicked, swapped sentences when
## they were identified correctly by the rule above. Points are sentences (one
## item at one length in one condition): the proportion identified correctly.
## The sentences of one item at one length share their x (the surprisal of the
## intended sentence), so the points are jittered sideways, by up to 2 bits,
## and never up or down, so that no point is drawn above 100%. Lines are
## logistic regressions of the trial-level outcome on surprisal, one per
## condition; their 95% confidence bands are computed on the logit scale, so
## they stay within 0-100%. They ignore that trials are grouped by participant
## and item. The bands are lighter than ggplot's default, so that where they
## overlap the conditions' colours stay apart. The legend calls the swaps
## "exchanges", as the text does.
p_correctness <- ggplot(mapping = aes(x = sum_surprisal, color = condition_label,
                                      fill = condition_label)) +
  geom_point(data = identification_trials %>%
               group_by(condition_label, item, length_level, sum_surprisal) %>%
               summarize(value = mean(correct), .groups = "drop"),
             aes(y = value), alpha = 0.4,
             position = position_jitter(width = 2, height = 0, seed = 2026)) +
  geom_smooth(data = identification_trials, aes(y = as.integer(correct)),
              method = "glm", method.args = list(family = binomial), formula = y ~ x,
              alpha = 0.2) +
  scale_color_manual(name = "Condition", values = CONDITION_COLORS,
                     labels = c("Original" = "Original", "1 swap" = "1 exchange",
                                "3 swaps" = "3 exchanges"),
                     aesthetics = c("colour", "fill")) +
  scale_x_continuous(name = "Sentence surprisal (bits)") +
  scale_y_continuous(name = "Proportion identified correctly", labels = scales::percent,
                     breaks = seq(0, 1, 0.25)) +
  coord_cartesian(ylim = c(0, 1))

if (interactive()) print(p_correctness)

ggsave(paste0("figures/s5_", task_name, "_correctness.png"),
       p_correctness, width = 8, height = 5.5)


## ---- Figure: acceptability on the trials identified correctly ----
## The human-judgment curves of the main text's qualitative-predictions figure,
## with its plotting code (p_baseline and p_gam in analysis/main_analysis.R),
## for the rated trials identified correctly. Points are sentence means (one
## item at one length in one condition); lines are GAM smooths of those means
## over surprisal, one per DL distance, with 95% confidence bands. Like the main
## figure, it is drawn at its print size, one PNAS column (8.7 cm) wide, in the
## main figure's 8-pt print theme (theme_print_8pt there). Where it differs
## from the main code: it has no "Human Judgments" title (there are no model
## columns to tell it from) and shows its own legend, at the bottom like the
## main figure's; the legend is centred on the whole plot, as the main figure's
## collected legend is, so that it fits at this width; and it leaves out "4+"
## and labels the "1-2" bucket "1", since these sentences have a DL distance of
## 0, 1 or 3 only.
plot_data_by_item <- analysis_trials %>%
  filter(subset == "Identified correctly") %>%
  group_by(subset, item, length_level, condition_label, sum_surprisal, dl_distance) %>%
  summarize(response = mean(acceptability), .groups = "drop") %>%
  mutate(dl_distance_bucket = bucket_dl(dl_distance))

p_acceptability <- ggplot(
  plot_data_by_item,
  aes(x = sum_surprisal, y = response, color = dl_distance_bucket)
) +
  geom_point(alpha = 0.5, size = 1.3, stroke = 0) +
  geom_smooth(
    method  = "gam",
    formula = y ~ s(x, bs = "tp"),
    se = TRUE,
    linewidth = 0.5,
    show.legend = FALSE
  ) +
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
    labels = c("0" = "0", "1-2" = "1", "3" = "3"),
    limits = DL_BUCKET_LEVELS,
    drop   = FALSE,
    guide  = guide_legend(
      override.aes = list(alpha = 1, linewidth = 0, size = 1.6)
    )
  )

if (interactive()) print(p_acceptability)

ggsave(paste0("figures/s5_", task_name, "_correct_trials.png"),
       p_acceptability, width = 3.42, height = 2.6, units = "in", dpi = 300)


print("Done.")


## #############################################################################
## Figure: the five models' log-likelihoods in the three studies (S11)
## #############################################################################
rm(list = setdiff(ls(), shared_objects))

## =============================================================================
## SI Section 5: the five models' log-likelihoods in the three studies of
## constructed sentences
##
## One row per study, drawn like the by-experiment panel of the main text's
## log-likelihood figure (p_language_likelihoods in analysis/main_analysis.R):
## each model's marginal log-likelihood minus that of the best model in the same
## study, with sentence surprisal as the information measure. The binary error
## detection + acceptability ratings study is shown on the trials judged
## correctly, and the exchange identification + acceptability ratings study on
## the trials identified correctly; the acceptability ratings only study has no
## second task, so all its trials are used.
##
## Reads the log-likelihoods that parts 5.2-5.4 above write, so it runs after
## them.
##
## Figure: figures/s5_likelihood_differences.png
## Numbers drawn: outputs/s5_likelihood_differences.csv
## =============================================================================

## ---- The rows of the figure ----
## The trials each study contributes, and its label, from top to bottom. The
## labels are broken over lines, as the main figure breaks its long row label
## ("English\n(controlled)"), with the trials used in parentheses on the last
## line.
figure_rows <- tribble(
  ~task,                                           ~subset,                ~dataset,
  "acceptability_ratings_only",                    "All trials",           "Acceptability ratings only\n(all trials)",
  "binary_error_detection_acceptability_ratings",  "Judged correctly",     "Binary error detection +\nacceptability ratings\n(trials judged correctly)",
  "exchange_identification_acceptability_ratings", "Identified correctly", "Exchange identification +\nacceptability ratings\n(trials identified correctly)"
)

likelihood_files <- paste0("outputs/s5_", figure_rows$task, "_likelihoods.csv")
if (!all(file.exists(likelihood_files))) {
  stop("Run parts 5.2-5.4 first; missing: ",
       paste(likelihood_files[!file.exists(likelihood_files)], collapse = ", "))
}

## Surprisal as the information measure; the length check stays in the CSVs
s5_likelihoods <- bind_rows(lapply(likelihood_files, read_csv, show_col_types = FALSE)) %>%
  filter(information == "Surprisal") %>%
  inner_join(figure_rows, by = c("task", "subset")) %>%
  mutate(
    ## The first study on top: a discrete y axis puts the first level at the
    ## bottom
    dataset = factor(dataset, levels = rev(figure_rows$dataset)),
    model = factor(model, levels = MODEL_LEVELS)
  )

stopifnot(
  nrow(s5_likelihoods) == nrow(figure_rows) * length(MODEL_LEVELS),
  all(s5_likelihoods$converged)
)

## The fits whose smoothing-parameter search did not reach full convergence with
## scam's default optimizer, and the optimizer of the fit reported (see (1) in
## parts 5.2-5.4)
s5_likelihoods %>%
  filter(bfgs_convergence != "Full convergence") %>%
  select(task, subset, model, bfgs_convergence, optimizer) %>%
  print(width = Inf)

s5_likelihoods %>%
  select(task, subset, model, log_likelihood, likelihood_diff) %>%
  mutate(log_likelihood = round(log_likelihood, 1),
         likelihood_diff = round(likelihood_diff, 1)) %>%
  print(n = Inf, width = Inf)

write_csv(
  s5_likelihoods %>%
    arrange(desc(dataset), model) %>%
    select(task, subset, model, log_likelihood, max_likelihood, likelihood_diff,
           n_participants, n_trials, n_cells),
  "outputs/s5_likelihood_differences.csv"
)

## ---- All five sets of trials ----
## The figure shows one set of trials per study; the studies with a second task
## were also analysed on all trials. In each of the five sets: the models from
## best to worst (surprisal as the information measure), and the length check:
## Communication's log-likelihood with length minus that with surprisal
s5_all_likelihoods <- bind_rows(lapply(likelihood_files, read_csv, show_col_types = FALSE))

s5_all_likelihoods %>%
  filter(information == "Surprisal") %>%
  arrange(task, subset, desc(log_likelihood)) %>%
  group_by(task, subset) %>%
  summarize(order = paste(model, collapse = " > "), .groups = "drop") %>%
  print(width = Inf)

s5_all_likelihoods %>%
  filter(model == "Communication") %>%
  select(task, subset, information, log_likelihood) %>%
  pivot_wider(names_from = information, values_from = log_likelihood) %>%
  mutate(length_minus_surprisal = Length - Surprisal) %>%
  print(width = Inf)

## ---- Figure ----
## The by-experiment panel of the main figure, with the legend of that figure's
## top panel, since this figure has no pooled panel to carry it. Five entries
## do not fit on one row at this width, so Binary, Gradient and Communication go
## on the first row and SLOR and MORCELA on the second. The figure is as wide as the
## main one (8 in), so that printed one PNAS column (8.7 cm) wide its smallest
## text (tick labels and legend) comes out at about 6 pt, as in the main
## figure. The long row labels push the panel to the right, and a legend
## centred over the panel would run off the right edge, so the legend is
## centred on the whole plot instead (legend.location = "plot").
p_dataset_likelihoods <- ggplot(s5_likelihoods) +
  geom_point(
    aes(
      y=dataset,
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
    name = "Log-Likelihood Difference"
  ) +
  ## The colour, fill and shape scales share a name, labels and guide, so ggplot
  ## merges them into one legend whose keys show each model's marker.
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE),
    fill = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme(legend.position="top", legend.location = "plot") +
  scale_y_discrete(name="")

if (interactive()) print(p_dataset_likelihoods)

ggsave("figures/s5_likelihood_differences.png",
       plot = p_dataset_likelihoods, width = 8, height = 5.5)


print("Done.")
