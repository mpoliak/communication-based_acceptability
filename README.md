# Evidence for a communication-based theory of linguistic acceptability judgments

Data and analysis code for the paper *"Evidence for a communication-based theory of
linguistic acceptability judgments."* The study comprises 9 acceptability-rating
experiments across 8 languages (Danish, English, French, German, Hindi, Korean,
Mandarin, Russian, plus a length-controlled English experiment). The Supporting
Information (SI) adds several English follow-up experiments (SI section 4).

In the main analyses, a sentence's information, Information(s), is its surprisal
under a language model (summed over the sentence, in bits). Five accounts of
acceptability are compared: binary grammaticality, gradient grammaticality, the
communication account, SLOR (Lau, Clark & Lappin 2017) and MORCELA (Tjuatja,
Neubig, Linzen & Hao 2025). SI section 3 repeats the same analyses with
Information(s) measured as sentence length in words, which is how the submitted
version of the paper measured it.

## Requirements

- **R** (developed with R 4.4) with: `tidyverse`, `lme4`, `lmerTest`, `mgcv`,
  `scam`, `patchwork`, `ggnewscale`; `stringdist` for stimulus preprocessing and
  SI sections 4.2–4.5.
  (The Mandarin materials were word-segmented with `jiebaR` beforehand and are
  distributed pre-segmented, so `jiebaR` is not needed to run the pipeline.)
- A **UTF-8 locale.** This is the default in RStudio. If you run the scripts with
  `Rscript` from a terminal, set one explicitly, e.g. `LC_ALL=en_US.UTF-8 Rscript ...`,
  otherwise non-Latin first-language/country values (Hindi, Korean, Cyrillic, `是`)
  will not compare correctly.
- To score the stimuli with the language models (step 1.3 below): a Python
  environment, given as a Conda file in `data/goldfish_environment.yml` (PyTorch,
  `transformers`, `wordfreq`).

Each script sets its own working directory, so it can be run from RStudio, via
`Rscript`, or `source()`d, without `setwd()`-ing by hand.

**From a clean clone.** The shipped `data/full_data.csv`,
`data/slor_morcela_scores.csv` and SI follow-up data folders are all the
analyses need; step 1 below only rebuilds the first two. From the repository
root, run:

```
LC_ALL=en_US.UTF-8 Rscript analysis/main_analysis.R
LC_ALL=en_US.UTF-8 Rscript supporting_information/run_all.R
```

The two are independent and can run in either order. Each script creates the
`models/`, `figures/` and `outputs/` folders it writes to.

## Repository structure

- `data/` — per-language raw experiment + materials data, the preprocessing
  scripts, the combined `full_data.csv`, the language-model scoring script
  `get_slor_morcela_scores.py` and its output `slor_morcela_scores.csv`, and the
  data of the SI follow-up experiments: `en_binary_error_detection/`,
  `en_experimental_binary_error_detection/` and `en_exchange_identification/`
  (section 4) and, on the 27 artificial items,
  `en_artificial_exchange_identification/` (4.1, 4.3),
  `en_artificial_binary_error_detection/` (4.2, 4.3),
  `en_artificial_acceptability_ratings/` (4.4) and
  `en_artificial_exchange_identification_acceptability/` (4.5). Each follow-up
  folder holds the responses, the Surveyor stimuli (`materials_*.csv`) and survey
  (`demographics_*.csv`) files and, where they exist, the `make_materials.R` that
  wrote them and the debriefing text. See `data/DATA_DICTIONARY.md`.
- `analysis/` — `main_analysis.R` (main-text analyses and figures), `figures/`,
  `outputs/` (the statistics reported in the paper, as CSV files) and `models/`
  (cached fits).
- `supporting_information/` — one script per SI section plus the shared `setup.R`,
  with `figures/`, `outputs/` and `models/` (cached fits).

## Reproducing the analyses

### 1. Prepare stimuli and data (under `data/`)
1. Run `preprocess_stimuli.R` and `preprocess_en_experimental.R` to compute
   utterance length and edit distance for the materials in all languages. This
   writes a `stimuli_data.csv` into each language folder. (`anonymize_participants.R`
   records how the experimenters replaced the Prolific IDs in the nine
   `experiment_data.csv` files, once, before release. It cannot be re-run: the raw
   IDs are not distributed, on the released files it would renumber the
   already-anonymized IDs, and it loops over every `data/` subfolder, so it stops at
   the first one without an `experiment_data.csv`.)
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
   artificial English follow-ups (SI sections 4.1–4.5), read from
   `en_artificial_acceptability_ratings/`.

   Every analysis takes its surprisal from this one file. The communication
   account's Information(s) is the surprisal of the item's original sentence,
   while SLOR and MORCELA use each condition's own scores (see the data
   dictionary). The script reads the `stimuli_data.csv` files from step 1, so
   rerun it whenever those change. The scored file is included in the
   repository.

   The first run downloads the Goldfish models from the Hugging Face Hub. Later
   runs use the downloaded copies but still need network access, because the
   script asks the Hub which safetensors revision to load; it fails with
   `HF_HUB_OFFLINE=1`. With the models already downloaded, scoring took a few
   minutes on a CPU.

### 2. Main analysis
Run `analysis/main_analysis.R`. Information(s) is sentence surprisal throughout.
The script:

- fits the mixed-effects models of DL distance × surprisal in each experiment
  (`outputs/lmm_dl_by_surprisal.csv`), and reports how strongly surprisal and
  length are correlated (`outputs/surprisal_length_correlations.csv`);
- tests the effect of surprisal on the error-free sentences across experiments
  (`outputs/error_free_surprisal_slopes.csv`, `outputs/error_free_surprisal_ttest.csv`);
- fits the five models (Binary, Gradient, Communication, SLOR and MORCELA) as
  shape-constrained GAMs, to each experiment and pooled over all nine, and
  compares their marginal log-likelihoods (`outputs/likelihoods_by_language.csv`,
  `outputs/likelihoods_all_experiments.csv`; MORCELA's fitted weights are in
  `outputs/morcela_weights.csv`).

| Figure | Contents |
|---|---|
| `figures/observed_vs_predicted.png` / `.pdf` | **Main figure.** Observed against predicted acceptability, one point per item and condition (DL distance 0, 1 or 3): the models pooled over all experiments (top row), then each experiment, with one column per model. The numbers are in `outputs/observed_vs_predicted.csv` and `outputs/observed_vs_predicted_correlations.csv` |
| `figures/likelihood_differences.png` | Model comparison: each model's log-likelihood relative to the best one, pooled over all experiments (top) and for each experiment (bottom) |
| `figures/mod_gradient_dl.png`, `figures/mod_communication_dlweighted_surprisal.png`, `figures/mod_communication_sum_surprisal.png` | The fitted smooths of the pooled gradient and communication models |
| `../supporting_information/figures/curves_surprisal.png` / `.pdf` | An SI figure (with its `_gam`, `_grammatical`, `_distance` and `_weighted_distance` panels): the ratings and the predictions of the binary, gradient and communication models over surprisal |

The pooled fits include all nine experiments, so unlike the "All Languages" row
of the submitted paper they include the controlled English experiment. Its three
length versions of each item share one item intercept.

**Fitting takes a long time.** The pooled five-model fit took about 30 minutes in
the one run that has been timed. `main_analysis.R` and SI `section3_length.R`
each contain one, on top of per-experiment fits and maximal-random-effects mixed
models that have not been timed. With empty model caches, budget several hours
for the whole pipeline. Later runs load the cached fits.

### 3. Supporting information
`supporting_information/` holds one script per section of the SI document. Each is
self-contained and can be run on its own, so reproducing one section does not mean
running the others. Sections 1–4 source `setup.R` (libraries, helper functions, the
`crit` data frame with its surprisal columns, the by-item frames and the shared
plotting helpers). Sections 4.1–4.5 do not source `setup.R`. They read the
`data/en_artificial_*` folders and, for 4.2–4.5, the sentence surprisals in
`data/slor_morcela_scores.csv`. Sections 4.2 and 4.5 each also read the other's
response file, to exclude the participants who took both studies, and 4.3 sources
4.1 and 4.2. Unless a row says otherwise, Information(s) is sentence surprisal, as
in the main text.

| Script | Contents |
|---|---|
| `setup.R` | Shared setup. Not run directly — sections 1–4, `diagnostics_and_extras.R` and `check_marginal_loglik.R` source it. |
| `section1_number_of_variables.R` | Section 1: controlling for the number of variables. The binary and gradient models get a surprisal term, so every model has as many predictors as the communication model |
| `section2_shuffled_condition.R` | Section 2: the model comparison with the shuffled condition included |
| `section3_length.R` | Section 3: length in words instead of surprisal. The main analyses (mixed models, the error-free length test, the five-model comparison) with Information(s) measured as sentence length, as in the submitted paper. Also draws the curve figure that was the submitted paper's main figure, over length (`figures/curves_length.png` / `.pdf`) |
| `../analysis/main_analysis.R` (last block) | The same curve figure over surprisal (`figures/curves_surprisal.png` / `.pdf`). It is drawn by the main analysis, from the main models |
| `section4_binary_error_detection.R` | Section 4: binary error detection — the English experiments in which participants also judged whether each sentence contained an error (`data/en_binary_error_detection/`, `data/en_experimental_binary_error_detection/`) — and the English exchange identification task (`data/en_exchange_identification/`). Information(s) is length here, as in the submitted SI. Figures `binary_error_detection_*.png`, `en_exchange_identification_*.png` and `en_swap_penalty_by_length.png` |
| `section4.1_exchange_identification_artificial.R` | Section 4.1: exchange identification on the 27 artificial items, each in a short, medium and long version: which words participants mark as misplaced. No model comparison. Data `data/en_artificial_exchange_identification/`. Figures `artificial_exchange_identification_{by_length,outcomes*,agreement*}.png` |
| `section4.2_binary_error_detection_artificial.R` | Section 4.2: binary error detection and acceptability on the artificial items. The models are compared with length and with surprisal as Information(s). Data `data/en_artificial_binary_error_detection/`. Figures `artificial_binary_error_detection_*.png` (not `_validated_`); `outputs/artificial_binary_error_detection_likelihoods.csv` |
| `section4.3_binary_error_detection_artificial_validated.R` | Section 4.3: section 4.2 on only the sentences validated in the exchange identification task (it sources 4.1 and 4.2). Length and surprisal. Figures `artificial_binary_error_detection_validated_*.png`; `outputs/artificial_binary_error_detection_validated_likelihoods.csv` |
| `section4.4_acceptability_ratings_artificial.R` | Section 4.4: acceptability alone on the artificial items. Length and surprisal. Data `data/en_artificial_acceptability_ratings/`. Figures `artificial_acceptability_ratings_{by_length,likelihood_differences}.png`; `outputs/artificial_acceptability_ratings_likelihoods.csv` |
| `section4.5_exchange_identification_acceptability_artificial.R` | Section 4.5: exchange identification followed by acceptability on the artificial items. Length and surprisal. Data `data/en_artificial_exchange_identification_acceptability/`. Figures `artificial_exchange_identification_acceptability_*.png` and `artificial_exchange_identification_ratings*.png`; `outputs/artificial_exchange_identification_acceptability_likelihoods.csv` |
| `diagnostics_and_extras.R` | Materials and participant descriptives, parametric/linear models, the maximal-exclusion robustness check, and diagnostic comparison figures. Not part of the SI document. |
| `check_marginal_loglik.R` | Checks `marginal_logLik()`, which the model comparisons use, against `lmer()` for each random-effect term the analyses use (item and language intercepts and slopes). Small fits; stops with an error if a test fails. Not run by `run_all.R` |
| `run_all.R` | Runs the ten section scripts above (sections 1–4.5 and `diagnostics_and_extras.R`) in document order. If a script fails, it reports the error, goes on with the next one, and at the end stops with the list of the scripts that failed. It does not run `main_analysis.R`, so it does not draw `curves_surprisal`. Section 4.3 sources 4.1 (output suppressed) and 4.2 itself, so 4.1 runs twice; it is deterministic, so the second run rewrites identical figures. 4.2 also runs twice, once on all sentences and once, through 4.3, on the validated ones. |

The numbers behind the SI figures are written to `supporting_information/outputs/`
(`s1_*`, `s2_*` and `s3_length_*.csv`, and `artificial_*_likelihoods.csv`, which
hold both information measures).

## Model caching

To save time, fitted models are cached as `.RData` files in each `models/` folder.
These caches are **git-ignored and not distributed**, so a fresh clone refits them
on the first run. To force a refit, delete the contents of the relevant `models/`
folder before running.

## Earlier versions

The submitted, length-based version of the main analysis, the earlier SLOR and
MORCELA comparison and the earlier surprisal pipeline were replaced in the
2026-10-09 restructure and remain in the git history (commit `e0bcc5d` and
earlier). Their analyses are now SI section 3 (length) and the main analysis
(surprisal, all five models). The last edits to the SLOR and MORCELA comparison,
which pooled all nine experiments instead of eight, were never committed, so the
history has only its 8-language version.

## Data dictionary

See [`data/DATA_DICTIONARY.md`](data/DATA_DICTIONARY.md) for a description of every
column in `full_data.csv`, in `slor_morcela_scores.csv` and in the data of the SI
follow-up experiments.

## License

MIT — see [`LICENSE`](LICENSE).
