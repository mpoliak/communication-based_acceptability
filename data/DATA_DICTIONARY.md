# Data dictionary — `full_data.csv`

`full_data.csv` is the analysis-ready dataset produced by
`concatenate_data_and_stimuli.R`, which joins each language's `experiment_data.csv`
(participants' responses + demographics) with its `stimuli_data.csv` (sentence
properties) and stacks all languages together. One row = one trial for one
participant. The file is UTF-8 encoded, with a byte-order mark; demographic free-text and dropdown
values appear in each experiment's own language/script.

Each participant contributes up to two rows per item: one acceptability-rating trial
(`experiment == "exp_b"`) and one comprehension/attention-check trial
(`experiment == "exp_b_comprehension"`). Demographic columns are repeated on
every row for a participant and are populated only for the experiments that
asked that particular question (so a given demographic column is `NA` for the
other languages — see notes below).

Some participants lack a few rows. Three Hindi participants stopped early
(`hi_26`, `hi_29` and `hi_33`, with 17, 29 and 44 rated items), and `zh_41` has
only 6 rated items.

The language-model scores that the analyses take their surprisal from,
`slor_morcela_scores.csv`, are described after it
([`slor_morcela_scores.csv`](#data-dictionary--slor_morcela_scorescsv)), followed
by the English scores under other language and unigram models that SI section 4
uses ([`slor_morcela_scores_english_settings.csv`](#data-dictionary--slor_morcela_scores_english_settingscsv)). The data
of the SI follow-up experiments (SI section 5) are not in `full_data.csv`. They
are described in the last section of this file
([SI follow-up experiments](#data-dictionary--si-follow-up-experiments-si-section-5)).

## Trial identification
| Column | Description |
|---|---|
| `participant_id` | Anonymized participant identifier (the Prolific ID was anonymized before release). |
| `experiment` | Trial type: `exp_b` = acceptability-rating trial; `exp_b_comprehension` = the comprehension question shown for that item. |
| `item` | Item identifier, prefixed by language (e.g., `en_8`). The English length-controlled experiment (`en_exp`) is prefixed with its folder name instead, `en_experimental_`, and a trailing `s` / `m` / `l` marks the short / medium / long version of the base item (e.g., `en_experimental_62m`). |
| `condition` | Stimulus condition: `original` (unaltered, grammatical), `swap_1` (one adjacent-token exchange), `swap_3` (three exchanges), `random` (fully shuffled). The two **nonsense** conditions — `non_sense` ("incoherent") and `non_sense_random` ("shuffled incoherent") — were run but are not analyzed in the paper. `condition` is the intended manipulation; the realized change is `dl_distance` (62 `swap_1` rating trials have DL distance 0, and `swap_3` trials range from 1 to 4). |
| `presentation_order` | 0-based position of the screen in the participant's session (each comprehension question is the screen after its rating). |
| `timestamp` | Time stamp of the participant's session (the same on all of a participant's rows). Dropped (so `NA`) for `zh` and `en_exp` during processing. |

## Stimulus
| Column | Description |
|---|---|
| `original_sentence` | The uncorrupted source sentence (from the Universal Dependencies corpus), lowercased, with punctuation and symbols removed (see Scored text below). |
| `perceived_sentence` | The (possibly corrupted) sentence actually shown to the participant, lowercased, with punctuation and symbols removed (see Scored text below). |
| `dl_distance` | Token-level Damerau–Levenshtein distance between `perceived_sentence` and `original_sentence` (number of token-level edits). |
| `n_words` | Utterance length: number of space-delimited tokens in the original sentence (Mandarin pre-tokenized with jiebaR). |
| `language` | Experiment: `da`, `de`, `en`, `en_exp` (English length-controlled), `fr`, `hi`, `ko`, `ru`, `zh`, or `NA`. |

> **All five stimulus columns above are `NA` on the nonsense trials.** The
> `non_sense` and `non_sense_random` conditions are dropped during stimulus
> preprocessing, so their trials carry no stimulus data after the join — 14,940
> rows (29.6% of the file) have `NA` for `original_sentence`, `perceived_sentence`,
> `dl_distance`, `n_words`, **and `language`**. These trials are not analyzed in
> the paper. Filter them out (e.g. `filter(!is.na(language))`) before any
> by-language summary, or a `group_by(language)` will silently drop them into an
> `NA` group.

### Provenance: Mandarin word segmentation

Mandarin is not written with spaces, so the Mandarin materials had to be
word-segmented before token-level edit distance and utterance length could be
computed. This was done **once, ahead of time**, with the `jiebaR` package, and
the result is distributed as `data/zh/materials_zh_for_analysis.csv` (the
unsegmented original remains in `data/zh/materials_zh.csv`). Because the
segmented file ships with the repository, `preprocess_stimuli.R` simply reads it
and `jiebaR` is *not* needed to reproduce the analyses.

The segmentation step that produced it was:

```r
library(jiebaR)

engine1 <- worker(bylines = TRUE)
sentences_segmented <- lapply(
  segment(sentences, engine1),
  function(x) x %>% unlist() %>% paste(collapse = " ")
) %>% unlist()
```

This is recorded for provenance. It is not re-run as part of the pipeline, and it
is not guaranteed to reproduce `materials_zh_for_analysis.csv` byte-for-byte,
since `jiebaR`'s bundled dictionaries and defaults may have changed since the
materials were prepared. Treat the shipped segmented file as the authoritative
version.

### Provenance: Russian materials

`data/ru/materials_ru.csv` is the authors' template of the Russian materials, not
the file exported from Surveyor; the two differ only in spaces before
punctuation marks and in two cells that no analysis uses.

## Response
| Column | Description |
|---|---|
| `response` | The participant's response: an acceptability rating from 1 (definitely not acceptable) to 7 (definitely acceptable) on rating trials; a (localized) yes/no answer on comprehension trials. **This is a character column, not numeric.** The two scale endpoints carry a localized label (e.g. `"1 - Definitely not acceptable"`, `"7 - 완전히 허용됨"`) while midpoints are bare digits (`"2"`…`"6"`), and comprehension trials hold text (`"Yes"`, `"Нет"`, `"对"`). Calling `as.numeric(response)` directly therefore silently returns `NA` for both endpoints. Filter to the rating trials and take the leading digit, exactly as the analysis scripts do: `filter(experiment == "exp_b")`, then `as.numeric(substr(response, 1, 1))`. |
| `literal_response` | On comprehension trials, the correct/expected answer (localized yes/no), used to score comprehension accuracy. |

## Demographics (raw, as collected per experiment)
The first-language question and the country-of-birth question were stored under
different column names / answer formats across experiments; the harmonized
columns at the bottom are the recommended ones for analysis.

| Column | Description |
|---|---|
| `first_language` | First language chosen from a localized dropdown (the language name, e.g., `Dansk`, `Deutsch`, `हिंदी`, `한국어`, `Русский`). Populated for `da`, `de`, `fr`, `hi`, `ko`, `ru`. |
| `first_language_text` | Free-text first language where "other" was permitted (`hi` only). |
| `english_first_language` | Yes/No to "is the study's language your first language" for `en`, `en_exp`, `zh` (Mandarin answers are `是` / `不是`). Despite the name, for `zh` it refers to Mandarin. Superseded by `first_language_is_target`. |
| `french_first_language` | `Oui` / `Non` to the same question in the French experiment. |
| `second_language` | Second language chosen from a dropdown (`fr`, `ru`). |
| `second_language_text` | Free-text second language (`da`, `de`, `hi`, `ko`, and a few `ru`). |
| `language_contexts` | Underscore-delimited multi-select of contexts in which the participant uses/used the language (e.g., `with parents_in school_at work_with friends`), for `da`, `de`, `hi`, `ko`, `ru`. |
| `french_contexts` | The French experiment's version of `language_contexts`. |
| `country_birth` | Country of birth as asked in the localized surveys (`da`, `de`, `hi`, `ru`); "Other" appears localized (e.g., `Andet`, `Anderes`). Raw — see `country_of_birth`. |
| `country_birth_text` | Free-text country of birth when "Other" was selected. |

## Harmonized demographic columns (recommended)
Added in `concatenate_data_and_stimuli.R` so that the differently-named/asked
fields above resolve to one canonical column each. The raw columns are retained
for provenance.

| Column | Description |
|---|---|
| `country_of_birth` | Canonical country of birth: `coalesce(country_of_birth, country_birth)`. Holds `country_of_birth` as collected in `en`/`en_exp`/`fr`/`ko`/`zh` (for `zh`, mostly a province or city, e.g. `广东`, `上海`) and the localized `country_birth` for `da`/`de`/`hi`/`ru`. |
| `first_language_is_target` | Logical. `TRUE` if the participant reported the study's language as their first language. Derived as `english_first_language %in% c("Yes","是")` **or** `french_first_language %in% "Oui"` **or** `first_language %in% c("Dansk","Deutsch","हिंदी","한국어","Русский")`. A multi-select answer that includes the target language (e.g. `अंग्रेजी_हिंदी`, English and Hindi) counts as `FALSE`. |

> **Note on locale.** Because the harmonized columns and several analyses compare
> non-Latin string literals (Hindi, Korean, Cyrillic, `是`), the R scripts must be
> run in a UTF-8 locale (the default in RStudio; for `Rscript` set, e.g.,
> `LC_ALL=en_US.UTF-8`). In a `C` locale these comparisons silently fail to match.

---

# Data dictionary — `slor_morcela_scores.csv`

`slor_morcela_scores.csv` holds the language-model scores of every stimulus. It is
written by `get_slor_morcela_scores.py` (run it from `data/` in the Conda
environment of `goldfish_environment.yml`), and it is the source of surprisal in
every analysis except SI section 4 under its other language models, which reads
`slor_morcela_scores_english_settings.csv` (below): the communication account's
Information(s), SLOR and MORCELA are all computed from it.

One row = one sentence version (an item in one condition) under one unigram model.
Every sentence therefore has two rows, one per value of `unigram_model`. The file
has 6,822 rows:

- 6,336 for the nine acceptability experiments: 288 sentences (72 items × 4
  conditions) in each of `da`, `de`, `en`, `fr`, `hi`, `ko`, `ru` and `zh`, and 864
  (216 item versions × 4 conditions) in `en_exp`, each × 2 unigram models;
- 486 for the artificial English follow-ups of SI section 5 (`en_artificial`):
  243 sentences (27 items × 9 conditions) × 2 unigram models.

| Column | Description |
|---|---|
| `language` | Experiment: `da`, `de`, `en`, `en_exp`, `fr`, `hi`, `ko`, `ru`, `zh` (as in `full_data.csv`), or `en_artificial` for the artificial English follow-ups. |
| `item` | Item identifier **without** the language prefix that `full_data.csv` uses: `1`–`72` in the eight languages, `1s`/`1m`/`1l` … `72s`/`72m`/`72l` in `en_exp` (short / medium / long version of the base item), and `1`–`27` in `en_artificial`. Read it as character, since the `en_exp` ids are not numbers. |
| `condition` | For the nine experiments, as in `full_data.csv`: `original`, `swap_1`, `swap_3`, `random`. The nonsense conditions are not scored. For `en_artificial`, the condition string of the materials, `<version>_<length>`: version `original`, `swap-1` or `swap-3` (with a hyphen), and length `short`, `medium` or `long` (e.g. `swap-1_medium`). |
| `unigram_model` | The unigram model behind `log_p_unigram`. `model_initial` (used by every analysis): each word's log2 probability under the same Goldfish model, scored alone at the start of a sentence. It exists in every language and gives every word a probability (the few characters outside the model's vocabulary are scored as `<unk>`). `wordfreq`: the log2 of the word's corpus frequency in the `wordfreq` package, with words it does not know given probability 1e-9. No model is fit on it here; SI section 4 fits wordfreq unigrams in English only. |
| `log_p_sentence` | Log2 probability of the whole sentence under the Goldfish model of its language, summed over its tokens, in bits. A `[CLS]` token is prepended as context, as the Goldfish model cards advise, but is not itself scored. The same in both of a sentence's rows. |
| `log_p_unigram` | Sum over the sentence's space-delimited words of each word's log2 unigram probability, in bits. |

**Sign convention.** Both score columns are log2 probabilities, so they are in
bits and always negative. Surprisal is their negation: a sentence's summed
surprisal is `-log_p_sentence`.

**Language models.** The Goldfish 1000mb model of each language, from the
Hugging Face Hub: `goldfish-models/dan_latn_1000mb` (`da`), `deu_latn_1000mb`
(`de`), `eng_latn_1000mb` (`en`, `en_exp` and `en_artificial`), `fra_latn_1000mb`
(`fr`), `hin_deva_1000mb` (`hi`), `kor_hang_1000mb` (`ko`), `rus_cyrl_1000mb`
(`ru`) and `zho_hans_1000mb` (`zh`).

**Scored text.** The sentence as the participant saw it, in every condition: the
`perceived_sentence` of each experiment's `stimuli_data.csv`. That column is
already normalized by `preprocess_stimuli.R`: punctuation and symbols removed
(`[\p{P}\p{S}]`), leading and trailing spaces trimmed, runs of spaces collapsed,
lowercased. The `en_artificial` sentences have no `stimuli_data.csv`, so the
script reads the critical (`exp_b`) rows of
`en_artificial_acceptability_ratings/materials_en_artificial_acceptability_ratings.csv`
and applies the same normalization itself. The three constructed-sentence
experiments of SI section 5 show the same 243 sentences, so they are scored once. In Mandarin, the jieba-segmented words
define the unigram terms, but the spaces between them are removed before the
sentence is given to the model for `log_p_sentence`.

**Joining to `full_data.csv`.** Keep the `model_initial` rows and rebuild the item
prefix, which for `en_exp` is the folder name:

```r
scores <- read_csv("slor_morcela_scores.csv", col_types = cols(item = col_character())) %>%
  filter(unigram_model == "model_initial") %>%
  mutate(item = paste(ifelse(language == "en_exp", "en_experimental", language),
                      item, sep = "_"))
```

**How the analyses use it.**

- `sum_surprisal`, the communication account's Information(s), is the surprisal
  of the intended sentence: `-log_p_sentence` of the item's `original` row,
  joined to every condition of the item on (`language`, `item`). Like `n_words`,
  it is therefore the same in all of an item's conditions. In `en_exp` the item
  id includes the `s`/`m`/`l` suffix, so each length version of a base item has
  its own value.
- `weighted_dl_distance_surprisal` = `dl_distance / sum_surprisal`, the surprisal
  counterpart of `weighted_dl_distance` (`dl_distance / n_words`).
- SLOR and MORCELA score the sentence that was actually read, so they use each
  condition's own `log_p_sentence` and `log_p_unigram`, joined on (`language`,
  `item`, `condition`). SLOR = (`log_p_sentence` − `log_p_unigram`) / `n_words`;
  MORCELA = (`log_p_sentence` − β·`log_p_unigram` + γ) / `n_words`, with β and γ
  fitted to the ratings (SLOR is β = 1, γ = 0).
- In the `en_artificial` analyses (SI section 5), `sum_surprisal` is
  `-log_p_sentence` of the `original_<length>` row, joined on (`item`, length),
  where length is the suffix of `condition`. The original and swapped versions of
  an item at one length share it; its short, medium and long versions differ.
- The length-controlled follow-ups of SI section 5 (5.2 and 5.3) show sentences
  of `en_exp` and use its rows, as the main analysis does.

---

# Data dictionary — `slor_morcela_scores_english_settings.csv`

The scores of SI section 4 (`supporting_information/section4_slor_morcela_settings.R`),
which asks whether SLOR or MORCELA does better with another language model or
unigram model than the main text's (the original sentences' `log_p_sentence` also
gives the Communication model's surprisal under each language model). The 288 sentences of the English experiment
(`en`: 72 items × 4 conditions) are scored under 18 settings, 6 language models ×
3 unigram models: 5,184 rows, one per sentence per setting. They are computed like
`slor_morcela_scores.csv`: on each condition's own `perceived_sentence` (from
`en/stimuli_data.csv`), in bits (log2), with words split at spaces.

| Column | Description |
|---|---|
| `language` | Always `en`. |
| `item` | `1`–`72`, without the `en_` prefix, as in `slor_morcela_scores.csv`. |
| `condition` | `original`, `swap_1`, `swap_3` or `random`. |
| `language_model` | The model behind `log_p_sentence`: `goldfish` (`goldfish-models/eng_latn_1000mb`, the main text's), `gpt2` and `gpt2-xl` (the Hugging Face GPT-2 models of those names), or `kn2`, `kn3` and `kn4` (interpolated Kneser-Ney 2-, 3- and 4-gram models, discount 0.75, trained on about 20 million words of WikiText-103 normalized like the stimuli). |
| `unigram_model` | The unigram model behind `log_p_unigram`: `model_initial` (each word scored alone at the start of a sentence by the same language model; the main text's), `wordfreq` (as in `slor_morcela_scores.csv`), or `corpus` (the word's relative frequency in the n-gram models' WikiText-103 text, an unseen word counting as half an occurrence). |
| `log_p_sentence` | Log2 probability of the sentence under `language_model`. The neural models condition on a start token that is not itself scored (`[CLS]` for Goldfish, `<\|endoftext\|>` for GPT-2), the n-gram models on `<s>` padding. No end-of-sentence token is scored. |
| `log_p_unigram` | Sum over the sentence's words of each word's log2 unigram probability. |

**Provenance.** Written by
`probability_based_models_scoring/get_slor_morcela_scores_english_settings.py`.
Run from that folder, in the Conda environment of `goldfish_environment.yml`, in
this order:

1. `python build_corpus.py <path to WikiText-103-raw-v1 wiki.train.raw>` writes
   `corpus/train.txt` and `corpus/unigram_counts.csv`: the first 20 million words
   of WikiText-103's training text, normalized like the stimuli, counting only
   sentences of 3 or more words.
2. `python train_kn_on_corpus.py` writes `models/kn2_train.pkl`, `kn3_train.pkl`
   and `kn4_train.pkl` (about 74, 361 and 983 MB; check free memory). The models
   are trained on all of that text except every 200th sentence, which is held
   out for the perplexity check it prints.
3. `python get_slor_morcela_scores_english_settings.py` reads
   `../en/stimuli_data.csv`, `corpus/` and `models/`, and overwrites
   `../slor_morcela_scores_english_settings.csv`.

`kn_ngram.py` (the Kneser-Ney model) and `neural_backends.py` (the Goldfish,
GPT-2 and GPT-2 XL scorers) are helpers the scripts import. The neural models
are loaded from the local Hugging Face cache only, so `goldfish-models/eng_latn_1000mb`,
`gpt2` and `gpt2-xl` must be downloaded first. `corpus/` and `models/` are
git-ignored. The shipped file was written by these scripts before they were
moved into the repository, with other column names; the scripts here differ
from that run only in their paths, comments and output format, and the scores
are the same. The Goldfish rows with
`model_initial` and `wordfreq` equal the `en` rows of `slor_morcela_scores.csv`
to about 1e-13 bits, which section 4 checks.

---

# Data dictionary — SI follow-up experiments (SI section 5)

Each folder below holds one Surveyor study of English speakers, run after the nine
acceptability experiments. None of them is part of `full_data.csv`; the SI
section 5 script reads six of them directly (see Used by). The seventh,
`en_artificial_exchange_identification/`, holds no responses, only the materials
that the three constructed-sentence studies are built from.

| Folder | Responses | Design | Used by |
|---|---|---|---|
| `en_binary_error_detection/` | `data_en_binary_error_detection.csv`, 80 participants (`en_awareness_<n>`) | The 72 items of the English experiment (`en`) in `original`, `swap_1`, `swap_3` and `random`. On each trial a binary error detection question (`exp_b_comprehension`), then a rating (`exp_b`) | SI 5.1, binary error detection + acceptability ratings on the original materials (`supporting_information/section5_awareness_of_grammatical_errors.R`, part 5.1) |
| `en_experimental_binary_error_detection/` | `data_en_exp_binary_error_detection.csv`, 40 participants (`en_exp_awareness_<nn>`) | The 864 sentences of the English (controlled) experiment (`en_exp`, materials in `en_experimental/`): 72 base items, each in a short, medium and long version, each `original`, `swap-1`, `swap-3` or `random`: 12 conditions (`<version>_<length>`). On each trial a binary error detection question (`exp_b_comprehension`), then a rating (`exp_b`) | SI 5.2, binary error detection + acceptability ratings on the length-controlled materials |
| `en_highlighted_errors/` | `data_en_highlighted_errors.csv`, 40 participants (`en_highlighted_<nn>`) | The 648 sentences of the English (controlled) experiment that are not shuffled (9 conditions), with every displaced word in bold red (`<b style="color:red;">`…`</b>` in `sentence`). Rating only (`exp_b`), with the prompt "Please rate how acceptable the sentence is (note that errors consisting of exchanged words are shown in red)" | SI 5.3, acceptability ratings with the errors highlighted |
| `en_artificial_exchange_identification/` | None (materials only) | 27 constructed items, each in a short, medium and long version, each `original`, `swap-1` or `swap-3`: 9 conditions (`<version>_<length>`), 243 sentences | The source of the materials of the three studies below |
| `en_artificial_binary_error_detection/` | `en_artificial_binary_error_detection_data.csv`, 100 participants | The same 243 sentences. Binary error detection, then rating | SI 5.5, binary error detection + acceptability ratings (5.6 reads its responses, to drop the people who took both) |
| `en_artificial_acceptability_ratings/` | `en_artificial_acceptability_ratings_data.csv`, 101 participants | The same 243 sentences. Rating only | SI 5.4, acceptability ratings only, and the source of the `en_artificial` sentences in `slor_morcela_scores.csv` |
| `en_artificial_exchange_identification_acceptability/` | `en_artificial_exchange_identification_acceptability_data.csv`, 100 participants | The same 243 sentences. Exchange identification, then rating | SI 5.6, exchange identification + acceptability ratings (5.5 reads its responses likewise) |

Three people took both `en_artificial_binary_error_detection` and
`en_artificial_exchange_identification_acceptability`. Sections 5.5 and 5.6 match
the two response files on `participant_id` and drop them from both.

## Response files

One row per trial (per screen). The follow-up question of a trial (the rating after
the binary error detection question or after exchange identification) is the
screen right after it, so its `presentation_order` is one higher.

| Column | Description |
|---|---|
| `participant_id` | Pseudonym that replaced the Prolific ID before release. The `en_artificial_*` studies share one numbering (`en_artificial_<n>`), so a person who took two of them has the same id in both; `en_binary_error_detection` is numbered on its own (`en_awareness_1` to `en_awareness_80`), as are `en_experimental_binary_error_detection` (`en_exp_awareness_01` to `en_exp_awareness_40`) and `en_highlighted_errors` (`en_highlighted_01` to `en_highlighted_40`). |
| `experiment` | Trial type: `exp_b` = acceptability rating (1–7, the same response format as in `full_data.csv`); `exp_b_comprehension` = here the binary error detection question, "Are there any errors in this sentence?" (`Yes`/`No`); `annot` = exchange identification, "Click on all the words that are misplaced". |
| `item` | Item number (`1`–`72`, or `1`–`27` in the artificial studies), without a language prefix. In the studies with length versions, it is the base item, and `cond` gives the length. |
| `cond` | Condition, as in the table above. |
| `presentation_order` | 0-based position of the screen in the participant's session. |
| `literal_response` | On the binary error detection question, the correct answer: `No` for originals, `Yes` for every altered sentence. `NA` on rating and click rows. |
| `response` | On rating rows, the rating (character, as in `full_data.csv`: take the leading digit). On binary error detection rows, `Yes` or `No`. On click rows, the clicked words as comma-separated `position:word` pairs, with 1-based word positions and the words percent-encoded for `%`, `,` and `:` (decode with `URLdecode`), or `NONE` when nothing was clicked. |
| `timestamp` | Time stamp of the participant's session (the same on all of a participant's rows). |
| `native_language` | `Yes`/`No` to "Is English your first language?". Sections 5.1–5.6 drop the participants who answered `No` (blank answers are kept). |
| `other_languages` | Free text: any languages other than English the participant speaks. |
| `country_of_birth` | Country of birth, chosen from a list. |

## Other files

- `materials_*.csv` — the Surveyor stimuli file (`experiment`, `item`,
  `condition`, `sentence`, `prompt`, `question_type`, `options`,
  `literal_response`, `config`). On two-screen trials the follow-up row is linked to
  the first screen through `config` (`coordinate`) and has an empty `sentence`.
- `demographics_*.csv` — the Surveyor survey-properties file: consent, the
  participant-ID field, the language and country questions, and the instructions
  (`name`, `type`, `content`, `options`, `config`).
- `debriefing_*.txt` — the debriefing shown at the end of the study.
- `make_materials.R` — writes that folder's `materials_*.csv` and
  `demographics_*.csv`. The materials are shipped, and the analyses do not
  regenerate them. `en_binary_error_detection/` has no `make_materials.R` and no
  debriefing file, `en_experimental_binary_error_detection/` and
  `en_highlighted_errors/` have no debriefing file, and
  `en_artificial_exchange_identification/` has no responses and no debriefing
  file. The two `make_materials.R` of the length-controlled studies read the
  sentences of `en_experimental/materials_en_experimental.csv`.
- `en_artificial_exchange_identification/fillers_exchanges.csv` — the table of
  the 243 artificial sentences, which that folder's `make_materials.R` reads. The
  authors constructed 27 base sentences in a short, medium and long version (the
  sentences were drafted with an AI assistant, then all checked and adjusted by
  the authors); the exchange positions were chosen by helper scripts (not kept)
  under the rules below and checked by hand. Columns: `item_id`, `length`,
  `version` (`original`, `1-exchange`, `3-exchange`), `structure` (a short label
  of the item's syntactic structure, e.g. `simple transitive, SVO + PP`, the same
  in all 9 rows of an item), `n_words`, `sentence`, `swap_slots` and
  `shared_1exch_across_lengths` (`True` if the 1-exchange swaps the same two
  words in the short, medium and long version of the item). `swap_slots` lists
  the exchanges, applied in order to the original: each value `i` exchanges
  words `i + 1` and `i + 2` (1-based). The exchange rules: each exchange swaps
  two adjacent words and must make the sentence ungrammatical on its own, so it
  moves a function word, and never moves an `-ly` adverb, *not*, a floating
  adverb (*now*, *still*, *just*, *only*, *even*, *ever*, *so*, *almost*,
  *nearly*, *always*), a subject–auxiliary pair (which would make a question) or
  a complementizer leftward past its verb (which would make a well-formed noun
  phrase); the only exception is the short version of item 10, which has too few
  such positions. The 3-exchange version contains the 1-exchange, applied first,
  and no word moves more than 2 positions.
  The other three artificial `make_materials.R` scripts read the materials file of
  `en_artificial_exchange_identification/`, so build that one first.
