# Evidence for a communication-based theory of linguistic acceptability judgments

Data and analysis code for the paper *"Evidence for a communication-based theory of
linguistic acceptability judgments."* The study comprises 9 acceptability-rating
experiments across 8 languages (Danish, English, French, German, Hindi, Korean,
Mandarin, Russian, plus a length-controlled English experiment). The Supporting
Information (SI) adds several English follow-up experiments (SI section 5).

In the main analyses, a sentence's information, Information(s), is its surprisal
under a language model (summed over the sentence, in bits). Five accounts of
acceptability are compared: binary grammaticality, gradient grammaticality, the
communication account, SLOR (Lau, Clark & Lappin 2017) and MORCELA (Tjuatja,
Neubig, Linzen & Hao 2025). SI section 3 repeats the same analyses with
Information(s) measured as sentence length in tokens.

## Requirements

- **R** (developed with R 4.4) with: `tidyverse`, `lme4`, `lmerTest`, `mgcv`,
  `scam`, `patchwork`, `ggnewscale`, `colorspace`; `stringdist` for stimulus preprocessing and
  SI section 5.
  (The Mandarin materials were word-segmented with `jiebaR` beforehand and are
  distributed pre-segmented, so `jiebaR` is not needed to run the pipeline.)
- A **UTF-8 locale.** This is the default in RStudio. If you run the scripts with
  `Rscript` from a terminal, set one explicitly, e.g. `LC_ALL=en_US.UTF-8 Rscript ...`,
  otherwise non-Latin first-language/country values (Hindi, Korean, Cyrillic, `是`)
  will not compare correctly.
- To score the stimuli with the language models (steps 1.3 and 1.4 below): a Python
  environment, given as a Conda file in `data/goldfish_environment.yml` (PyTorch,
  `transformers`, `wordfreq`).

Each script sets its own working directory, so it can be run from RStudio, via
`Rscript`, or `source()`d, without `setwd()`-ing by hand.

**From a clean clone.** The shipped `data/full_data.csv`,
`data/slor_morcela_scores.csv`, `data/slor_morcela_scores_english_settings.csv`
and SI follow-up data folders are all the analyses need; steps 1.1–1.3 below
rebuild the first two, and step 1.4 the third (see also the data dictionary).
From the repository root, run:

```
LC_ALL=en_US.UTF-8 Rscript analysis/main_analysis.R
LC_ALL=en_US.UTF-8 Rscript supporting_information/run_all.R
```

Running them from inside their own folders works too, e.g.
`cd supporting_information && LC_ALL=en_US.UTF-8 Rscript run_all.R`.

The two can run in either order: SI section 4 checks its fits against two
shipped outputs of the main analysis (`analysis/outputs/likelihoods_by_language.csv`
and `morcela_weights.csv`), so after changing `main_analysis.R`, rerun it first.
Each script creates the `models/`, `figures/` and `outputs/` folders it writes to.

## Repository structure

- `data/` — per-language raw experiment + materials data, the preprocessing
  scripts, the combined `full_data.csv`, the language-model scoring script
  `get_slor_morcela_scores.py` and its output `slor_morcela_scores.csv`, the
  English scores under other language and unigram models that SI section 4 uses
  (`slor_morcela_scores_english_settings.csv`) and the scripts that compute them
  (`probability_based_models_scoring/`), the SLOR and MORCELA values of the main
  text's two example sentences (`sample_slor_morcela_scores.py` and its output
  `sample_slor_morcela_scores.csv`), and the
  data of the SI follow-up experiments: on the 27 artificial items,
  `en_artificial_acceptability_ratings/` (5.2, acceptability ratings only),
  `en_artificial_binary_error_detection/` (5.3, binary error detection +
  acceptability ratings) and
  `en_artificial_exchange_identification_acceptability/` (5.4, exchange
  identification + acceptability ratings), and, on the 72 items of the English
  experiment, `en_binary_error_detection/` (5.1, binary error detection +
  acceptability ratings). The folder `en_artificial_exchange_identification/`
  holds no responses, only the materials that the three artificial studies'
  materials are built from (`fillers_exchanges.csv`, the `make_materials.R` that
  reads it, and the Surveyor files it writes). Each of the four follow-up
  folders holds the responses, the Surveyor stimuli (`materials_*.csv`) and survey
  (`demographics_*.csv`) files and, where they exist, the `make_materials.R` that
  wrote them and the debriefing text. See `data/DATA_DICTIONARY.md`.
- `analysis/` — `main_analysis.R` (main-text analyses and figures), `figures/`,
  `outputs/` (the statistics reported in the paper, as CSV files) and `models/`
  (cached fits).
- `supporting_information/` — one script per SI section plus the shared
  `setup.R`, and `descriptives.R`, `check_marginal_loglik.R`, `run_all.R` and
  `multilingual_instructions_and_prompts.txt` (the instructions and rating prompt
  of the main experiments in the 8 languages), with `figures/`, `outputs/` and
  `models/` (cached fits).

## Reproducing the analyses

### 1. Prepare stimuli and data (under `data/`)
1. Run `preprocess_stimuli.R` and `preprocess_en_experimental.R` to compute
   utterance length and edit distance for the materials in all languages. This
   writes a `stimuli_data.csv` into each language folder. (`anonymize_participants.R`
   records how the experimenters replaced the Prolific IDs in the nine
   `experiment_data.csv` files, once, before release. It cannot be re-run: the raw
   IDs are not distributed, on the released files it would renumber the
   already-anonymized IDs, and it loops over every `data/` subfolder, so it stops at
   the first one without an `experiment_data.csv`. The Prolific IDs of the SI
   follow-up studies were replaced before release by code that is not in this
   repository.)
2. Run `concatenate_data_and_stimuli.R` to join the experimental data with the
   stimulus data and concatenate all languages into `full_data.csv` (this also
   adds the harmonized `first_language_is_target` and unified `country_of_birth`
   columns — see the data dictionary).
3. Score the stimuli with the language models. In the Conda environment of
   `goldfish_environment.yml`, run `cd data && python get_slor_morcela_scores.py`.
   It scores every sentence as the participant saw it (`perceived_sentence`:
   lowercased, without punctuation), in every condition, with the Goldfish
   1000mb model of its language, and writes `slor_morcela_scores.csv`: each
   sentence's log2 probability and the summed unigram log2 probability of its
   words. Besides the nine experiments, it scores the 243 sentences of the
   artificial English follow-ups (SI section 5), read from
   `en_artificial_acceptability_ratings/`.

   Apart from SI section 4's other language models, every analysis takes its
   surprisal from this one file. (SI section 4 also scores SLOR and MORCELA under
   other language and unigram models, and refits Communication with the other
   language models' surprisal of the original sentence, from the shipped
   `slor_morcela_scores_english_settings.csv`; see step 4 and the data
   dictionary.) The communication
   account's Information(s) is the surprisal of the item's original sentence,
   while SLOR and MORCELA use each condition's own scores (see the data
   dictionary). The script reads the `stimuli_data.csv` files from step 1, so
   rerun it whenever those change. The scored file is included in the
   repository.

   `python sample_slor_morcela_scores.py`, in the same environment, computes the
   SLOR and MORCELA values the main text gives for its two example sentences
   ("ball" and "demonstrated", in all four versions), writes them to
   `sample_slor_morcela_scores.csv` and prints them rounded. It scores them the
   same way, with MORCELA's English weights read from
   `analysis/outputs/morcela_weights.csv` (row `en`), which
   `analysis/main_analysis.R` writes.

   The first run downloads the Goldfish models from the Hugging Face Hub. Later
   runs use the downloaded copies but still need network access, because the
   script asks the Hub which safetensors revision to load; it fails with
   `HF_HUB_OFFLINE=1`. On import, PyTorch prints a NumPy 1.x/2.x compatibility
   message that looks like an error; it is harmless here (the scores are not
   affected).
4. The English scores of SI section 4 under other language and unigram models,
   `slor_morcela_scores_english_settings.csv`, are written by the scripts in
   `probability_based_models_scoring/`. In the same environment, from that
   folder, run in this order:
   - `python build_corpus.py <path to WikiText-103-raw-v1 wiki.train.raw>`,
     which writes `corpus/`: the first 20 million words of WikiText-103,
     normalized like the stimuli, and their unigram counts;
   - `python train_kn_on_corpus.py`, which writes the Kneser-Ney 2-, 3- and
     4-gram models to `models/` (the 4-gram file alone is about 1 GB, so check
     free memory);
   - `python get_slor_morcela_scores_english_settings.py`, which overwrites
     `../slor_morcela_scores_english_settings.csv`.

   The GPT-2, GPT-2 XL and English Goldfish models are loaded from the local
   Hugging Face cache only, so download them first. `kn_ngram.py` and
   `neural_backends.py` are helpers that the other scripts import. `corpus/` and
   `models/` are git-ignored.

### 2. Main analysis
Run `analysis/main_analysis.R`. Information(s) is sentence surprisal throughout.
The script:

- fits the mixed-effects models of DL distance × surprisal in each experiment
  (`outputs/lmm_dl_by_surprisal.csv`), and reports how strongly surprisal and
  length are correlated (`outputs/surprisal_length_correlations.csv`);
- tests the effect of surprisal on the error-free sentences across experiments
  (`outputs/error_free_surprisal_slopes.csv`, `outputs/error_free_surprisal_ttest.csv`);
- fits the five models (Binary, Gradient, Communication, SLOR and MORCELA) as
  shape-constrained GAMs, to each experiment and pooled over the eight language
  experiments, and compares their marginal log-likelihoods (`outputs/likelihoods_by_language.csv`,
  `outputs/likelihoods_all_experiments.csv`; MORCELA's fitted weights are in
  `outputs/morcela_weights.csv`).

| Figure | Contents |
|---|---|
| `figures/qualitative_predictions.png` | **Main figure: qualitative predictions.** One row per experiment. Columns: the human judgments, then the predictions of the binary, gradient and communication models, over sentence surprisal. Points are item means by DL distance (0, 1–2, 3, 4+); the human-judgment curves are GAM smooths of those points, and the model curves are each model's predictions (±2 SE) at DL distance 0, 1 and 3. Drawn at print size: 7 × 7.5 in (17.8 × 19.1 cm), 8-pt text |
| `figures/observed_vs_predicted.png` | **Main figure: quantitative predictions.** Observed against predicted acceptability, one point per item and condition (DL distance 0, 1 or 3), for all five models (columns: Binary Grammaticality, Gradient Grammaticality, Communication, SLOR, MORCELA): the models pooled over the eight language experiments (top row), then each experiment. Colour shows the DL distance and, within it, sentence surprisal (light = low, dark = high). 7 × 7.6 in. The numbers are in `outputs/observed_vs_predicted.csv`, `outputs/observed_vs_predicted_correlations.csv` (r and RMSE per panel) and `outputs/observed_vs_predicted_errors.csv` (mean error, observed − predicted, per panel at DL 0, 1, 3 and 1 or 3) |
| `figures/likelihood_differences.png` | **Main figure: model comparison.** Each model's marginal log-likelihood relative to the best one, pooled over the eight language experiments (top) and for each experiment (bottom). 8 × 8 in, like the SI's likelihood figures; printed one column wide, its smallest text is about 6 pt. The numbers are in `outputs/likelihoods_all_experiments.csv` and `outputs/likelihoods_by_language.csv` |
| `figures/residuals_by_surprisal.png` | **Main figure: residuals.** Observed minus predicted acceptability against sentence surprisal, for all five models (columns: Binary Grammaticality, Gradient Grammaticality, Communication, SLOR, MORCELA) and every row of the figure above, with a smooth line per DL distance (0, 1, 3). Above zero: the model underpredicts those sentences. 7 × 7.9 in |

The pooled fits, here and in SI sections 1–3, include the eight language
experiments and leave out the controlled English experiment; the figures
label them "All Languages". The controlled English experiment is fit only on its own, where its
three length versions of each item share one item intercept. Names that contain
`all_experiments` (output files, CSV rows, model caches) refer to this
eight-language pool.

**Fitting takes a long time.** `main_analysis.R` and SI sections 1, 2 and 3 each
fit the five models pooled over the eight language experiments, and a fresh
clone has no cached fits (see Model caching), so its first run fits every
model. Later runs load the cached fits.

### 3. Supporting information
`supporting_information/` holds one script per section of the SI document. Each is
self-contained and can be run on its own, so reproducing one section does not mean
running the others. Sections 1–5 and `descriptives.R` source `setup.R` (libraries,
helper functions, the `crit` data frame with its surprisal columns, the by-item
frames and the shared plotting helpers). The section 5 script runs its parts in
order: 5.1, then a shared set-up of the constructed sentences, which builds the
table of the 243 artificial sentences with their scores from
`data/slor_morcela_scores.csv`, then 5.2–5.4 and the figure that combines their
log-likelihoods. Parts 5.3 and 5.4 each also read the other's response file, to
exclude the participants who took both studies. Unless a row says otherwise, Information(s) is sentence surprisal, as
in the main text.

| Script | Contents |
|---|---|
| `setup.R` | Shared setup. Not run directly — sections 1–5, `descriptives.R` and `check_marginal_loglik.R` source it. |
| `section1_number_of_variables.R` | Section 1: controlling for the number of variables. The binary and gradient models, and SLOR and MORCELA, get a surprisal term (monotone decreasing, as in the communication model), so every model has as many predictors as the communication model; all five models are compared. Figures `with_surprisal_term_figure.png` and `likelihood_differences_with_surprisal_term.png`; `outputs/s1_surprisal_term_likelihoods_pooled.csv` and `_by_language.csv` |
| `section2_shuffled_condition.R` | Section 2: the five-model comparison with the shuffled condition included. MORCELA's weights are learned on these data (`outputs/s2_morcela_weights.csv`). Figures `with_shuffle_figure.png` and `likelihood_differences_with_shuffled.png`; `outputs/s2_shuffled_surprisal_likelihoods_pooled.csv` and `_by_language.csv` |
| `section3_length.R` | Section 3: length in tokens instead of surprisal. The main analyses (mixed models, the error-free length test, the five-model comparison) with Information(s) measured as sentence length. Figures `curves_length.png` (the main text's qualitative-predictions figure, over length) and `likelihood_differences_length.png`; `outputs/s3_length_*.csv` |
| `section4_slor_morcela_settings.R` | Section 4: SLOR and MORCELA in the English experiment under 18 settings (6 language models × 3 unigram models), from `data/slor_morcela_scores_english_settings.csv`, against the main text's other three models; Communication is refit under each of the 6 language models, with that model's surprisal of the original sentence. Stops unless the main text's setting reproduces `../analysis/outputs/likelihoods_by_language.csv` and `morcela_weights.csv`. Also prints how the neural language models rank, at each unigram model and at each one's best setting, and the best unigram model for each language model. Figure `figures/s4_slor_morcela_settings.png`; `outputs/s4_slor_morcela_settings_likelihoods.csv` |
| `section5_awareness_of_grammatical_errors.R` | Section 5: probing participants' awareness of grammatical errors, one script in parts, in the order of the SI. |
| ↳ Part 5.1 | Binary error detection + acceptability ratings, on the original materials of the English experiment (its 72 items, original, with 1 swap, with 3 swaps or shuffled): participants first judge whether the sentence contains an error, then rate it. The same analyses as 5.3, on all trials and on the trials judged correctly, with the English experiment's specifications (by-item intercepts in the models, DL distance up to 3; binary error detection is reported for the shuffled sentences too). Starts from `setup.R`, not from the shared set-up of 5.2–5.4. Data `data/en_binary_error_detection/`. Figures `figures/s5_original_binary_error_detection_yes_responses` (proportion of "yes" answers to the binary error detection question) and `_correct_trials` (`.png`); `outputs/s5_original_binary_error_detection_*.csv` |
| ↳ Shared set-up of 5.2–5.4 | It sources `setup.R` and builds the table of the 243 artificial sentences that all three studies show (27 items, each in a short, medium and long version, each original, with 1 swap or with 3 swaps), with their DL distance and their scores from `data/slor_morcela_scores.csv` |
| ↳ Part 5.2 | Acceptability ratings only, on the artificial items. The main text's analyses: the five-model comparison, the DL distance × surprisal mixed model, the error-free surprisal slope and, in the CSV only, Communication with length in place of surprisal. Data `data/en_artificial_acceptability_ratings/`. Figure `figures/s5_acceptability_ratings_only_all_trials.png`; `outputs/s5_acceptability_ratings_only_*.csv` |
| ↳ Part 5.3 | Binary error detection + acceptability ratings, on the artificial items: participants first judge whether the sentence contains an error, then rate it. The same analyses, on all trials and on the trials judged correctly. Data `data/en_artificial_binary_error_detection/`. Figures `figures/s5_binary_error_detection_acceptability_ratings_yes_responses` (proportion of "yes" answers to the binary error detection question) and `_correct_trials` (`.png`); `outputs/s5_binary_error_detection_acceptability_ratings_*.csv` |
| ↳ Part 5.4 | Exchange identification + acceptability ratings, on the artificial items: participants first click the words they think are misplaced, then rate the sentence. The same analyses, on all rated trials and on the trials identified correctly (the rule is in the script's header). Data `data/en_artificial_exchange_identification_acceptability/`. Figures `figures/s5_exchange_identification_acceptability_ratings_correctness` and `_correct_trials` (`.png`); `outputs/s5_exchange_identification_acceptability_ratings_*.csv` |
| ↳ Final figure | The five models' log-likelihoods in the three studies, one row per study (5.2 on all trials, 5.3 on the trials judged correctly, 5.4 on the trials identified correctly): `figures/s5_likelihood_differences.png`, `outputs/s5_likelihood_differences.csv`. Reads the log-likelihoods that 5.2–5.4 write, so run it after them |
| `descriptives.R` | Materials and participant descriptives the paper cites: sentence length per experiment, comprehension accuracy by DL distance and the answers to the first-language question (printed), and the item × condition cells by DL distance (`outputs/dl_distance_counts.txt`). Not part of the SI document. |
| `check_marginal_loglik.R` | Checks `marginal_logLik()`, which the model comparisons use, against `lmer()` for each random-effect term the analyses use (item and language intercepts and slopes). Small fits; stops with an error if a test fails. Not run by `run_all.R` |
| `run_all.R` | Runs six scripts above (sections 1–5 and `descriptives.R`) in document order, each in its own R process and with its own log, `logs/<script name>.log` (it creates `logs/` at the repository root). If a script fails, the others still run (its error is at the end of its log), and at the end it stops with the list of the scripts that failed. It does not run `main_analysis.R`, so it does not draw `likelihood_differences.png`. |

The numbers behind the SI figures are written to `supporting_information/outputs/`
(`s1_*`, `s2_*` (including `s2_morcela_weights.csv`), `s3_length_*`, `s5_*` (the
`s5_*_likelihoods.csv` hold both information measures) and `s4_*` CSVs;
`descriptives.R` writes `dl_distance_counts.txt`). In
`s4_slor_morcela_settings_likelihoods.csv` and the `s5_*_likelihoods.csv`,
`bfgs_convergence` records how scam's search for each fit's smoothing parameters
ended, and `optimizer` which fit the row reports: a fit whose search stopped
short is refit with `nlm` and `optim`, and the fit with the lowest smoothness
score (GCV) is kept.

## Model caching

To save time, fitted models are cached as `.RData` files in each `models/` folder.
These caches are **git-ignored and not distributed**, so a fresh clone refits them
on the first run. To force a refit, delete the contents of the relevant `models/`
folder before running. Do this also after rebuilding `full_data.csv` or
`slor_morcela_scores.csv`: the scripts check only the names of the cached fits,
not the data they were fit to.

## Earlier versions

The submitted, length-based version of the main analysis, the earlier SLOR and
MORCELA comparison and the earlier surprisal pipeline were replaced in the
2026-10-09 restructure and remain in the git history. Their analyses are now SI
section 3 (length) and the main analysis (surprisal, all five models). SI
section 3's `curves_length.png` was the submitted paper's main figure, and the
submitted paper's "All Languages" row pooled the same eight language experiments
as the pooled fits here. The curve figure over surprisal
(`analysis/figures/qualitative_predictions.png`) and the five-model
log-likelihood figure (`analysis/figures/likelihood_differences.png`), formerly
in the SI, are now main-text figures. The SI section 4 of before 2026-11-02
(the English follow-ups of `section4_binary_error_detection.R` and the
artificial-item sections 4.1–4.5) and `diagnostics_and_extras.R` were removed in
the 2026-11-02 rebuild of the SI; their scripts, figures and CSVs remain in the
git history. The parts of SI section 5 were separate scripts until 2026-11-22.

## Data dictionary

See [`data/DATA_DICTIONARY.md`](data/DATA_DICTIONARY.md) for a description of every
column in `full_data.csv`, in `slor_morcela_scores.csv`, in
`slor_morcela_scores_english_settings.csv` and in the data of the SI follow-up
experiments.

## License

MIT — see [`LICENSE`](LICENSE).
