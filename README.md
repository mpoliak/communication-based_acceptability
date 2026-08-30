# Evidence for a communication-based theory of linguistic acceptability judgments

Data and analysis code for the paper *"Evidence for a communication-based theory of
linguistic acceptability judgments."* The study comprises 9 acceptability-rating
experiments across 8 languages (Danish, English, French, German, Hindi, Korean,
Mandarin, Russian, plus a length-controlled English experiment).

## Requirements

- **R** (developed with R 4.4) with: `tidyverse`, `lme4`, `lmerTest`, `mgcv`,
  `scam`, `patchwork`, `ggnewscale`; `stringdist` for stimulus preprocessing.
  (The Mandarin materials were word-segmented with `jiebaR` beforehand and are
  distributed pre-segmented, so `jiebaR` is not needed to run the pipeline.)
- A **UTF-8 locale.** This is the default in RStudio. If you run the scripts with
  `Rscript` from a terminal, set one explicitly, e.g. `LC_ALL=en_US.UTF-8 Rscript ...`,
  otherwise non-Latin first-language/country values (Hindi, Korean, Cyrillic, `是`)
  will not compare correctly.
- For the surprisal analysis only: a Python environment (see
  `supporting_information/goldfish_surprisals/goldfish_environment.yml`).

Each script sets its own working directory, so it can be run from RStudio, via
`Rscript`, or `source()`d, without `setwd()`-ing by hand.

## Repository structure

- `data/` — per-language raw experiment + materials data, the preprocessing
  scripts, and the combined `full_data.csv` (see `data/DATA_DICTIONARY.md`).
- `analysis/` — `main_analysis.R` (paper analyses and main figures) and `models/`
  (cached fits).
- `supporting_information/` — one script per SI section plus the shared `setup.R`
  they all source, `models/` (cached fits), and `goldfish_surprisals/` (surprisal
  pipeline).
- `slor_morcela/` — comparison against two language-model accounts of
  acceptability, SLOR (Lau, Clark & Lappin 2017) and MORCELA (Tjuatja, Neubig,
  Linzen & Hao 2025). See `slor_morcela/README.md`.

## Reproducing the analyses

### 1. Prepare stimuli and data (under `data/`)
1. Run `preprocess_stimuli.R` and `preprocess_en_experimental.R` to compute
   utterance length and edit distance for the materials in all languages. This
   writes a `stimuli_data.csv` into each language folder. (`anonymize_participants.R`
   was used once by the experimenters to anonymize Prolific IDs; it does not need
   to be re-run.)
2. Run `concatenate_data_and_stimuli.R` to join the experimental data with the
   stimulus data and concatenate all languages into `full_data.csv` (this also
   adds the harmonized `first_language_is_target` and unified `country_of_birth`
   columns — see the data dictionary).

### 2. (Optional) Compute surprisal, for the SI surprisal analysis
3. Run `combine_stimuli_for_surprisal.R`, copy the resulting
   `materials_for_surprisal.csv` into `supporting_information/goldfish_surprisals/`,
   and in that folder run `get_surprisals_goldfish.py` (Conda env in
   `goldfish_environment.yml`). The output `surprisal_goldfish.csv` is already
   included, so this step is only needed to regenerate it.

### 3. Main analysis
Run `analysis/main_analysis.R` to produce the paper's figures and analyses.

### 4. SLOR and MORCELA
Run `slor_morcela/get_slor_morcela_scores.py` once to score the stimuli with the
Goldfish models (it needs the same Conda environment as the surprisal pipeline,
plus `wordfreq`), then `slor_morcela/slor_morcela_analysis.R` to fit the models
and draw the figures.

### 5. Supporting information
`supporting_information/` holds one script per section of the SI document. Each is
self-contained: it sources `setup.R` (libraries, helper functions, the `crit` data
frame, the by-item frames and the shared plotting helpers) and can be run on its
own, so reproducing one section does not mean running the others.

| Script | Contents |
|---|---|
| `setup.R` | Shared setup. Not run directly — every script below sources it. |
| `section1_number_of_variables.R` | Section 1: controlling for the number of variables |
| `section2_shuffled_condition.R` | Section 2: including the shuffled condition |
| `section3_surprisal.R` | Section 3: surprisal instead of length |
| `section4_error_awareness.R` | Section 4: error awareness — the English experiment in which participants also judged whether each sentence contained an error |
| `diagnostics_and_extras.R` | Materials and participant descriptives, parametric/linear models, the maximal-exclusion robustness check, and diagnostic comparison figures. Not part of the SI document. |
| `run_all.R` | Runs all five in document order. |

## Model caching

To save time, fitted models are cached as `.RData` files in each `models/` folder.
These caches are **git-ignored and not distributed**, so a fresh clone refits them
on the first run. To force a refit, delete the contents of the relevant `models/`
folder before running.

## Data dictionary

See [`data/DATA_DICTIONARY.md`](data/DATA_DICTIONARY.md) for a description of every
column in `full_data.csv`.

## License

MIT — see [`LICENSE`](LICENSE).
