# Data dictionary — `full_data.csv`

`full_data.csv` is the analysis-ready dataset produced by
`concatenate_data_and_stimuli.R`, which joins each language's `experiment_data.csv`
(participants' responses + demographics) with its `stimuli_data.csv` (sentence
properties) and stacks all languages together. One row = one trial for one
participant. The file is UTF-8 encoded; demographic free-text and dropdown
values appear in each experiment's own language/script.

Each participant contributes two rows per item: one acceptability-rating trial
(`experiment == "exp_b"`) and one comprehension/attention-check trial
(`experiment == "exp_b_comprehension"`). Demographic columns are repeated on
every row for a participant and are populated only for the experiments that
asked that particular question (so a given demographic column is `NA` for the
other languages — see notes below).

## Trial identification
| Column | Description |
|---|---|
| `participant_id` | Anonymized participant identifier (the Prolific ID was anonymized before release). |
| `experiment` | Trial type: `exp_b` = acceptability-rating trial; `exp_b_comprehension` = the comprehension question shown for that item. |
| `item` | Item identifier, prefixed by language (e.g., `en_8`). For the English length-controlled experiment (`en_exp`), a trailing `s` / `m` / `l` marks the short / medium / long version of the base item. |
| `condition` | Stimulus condition: `original` (unaltered, grammatical), `swap_1` (one adjacent-token exchange), `swap_3` (three exchanges), `random` (fully shuffled). The two **nonsense** conditions — `non_sense` ("incoherent") and `non_sense_random` ("shuffled incoherent") — were run but are not analyzed in the paper. |
| `presentation_order` | Position (order) in which the item was presented to the participant. |
| `timestamp` | Trial timestamp. Dropped (so `NA`) for `zh` and `en_exp` during processing. |

## Stimulus
| Column | Description |
|---|---|
| `original_sentence` | The uncorrupted source sentence (from the Universal Dependencies corpus). |
| `perceived_sentence` | The (possibly corrupted) sentence actually shown to the participant. |
| `dl_distance` | Token-level Damerau–Levenshtein distance between `perceived_sentence` and `original_sentence` (number of token-level edits). |
| `n_words` | Utterance length: number of space-delimited tokens in the original sentence (Mandarin pre-tokenized with jiebaR). |
| `language` | Experiment: `da`, `de`, `en`, `en_exp` (English length-controlled), `fr`, `hi`, `ko`, `ru`, `zh`, or `NA`. |

> **All five stimulus columns above are `NA` on the nonsense trials.** The
> `non_sense` and `non_sense_random` conditions are dropped during stimulus
> preprocessing, so their trials carry no stimulus data after the join — 15,036
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
| `first_language_text` | Free-text first language where "other" was permitted (mostly `hi`). |
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
| `country_of_birth` | Canonical country of birth: `coalesce(country_of_birth, country_birth)`. Holds `country_of_birth` as collected in `en`/`en_exp`/`fr`/`ko`/`zh` and the localized `country_birth` for `da`/`de`/`hi`/`ru`. |
| `first_language_is_target` | Logical. `TRUE` if the participant reported the study's language as their first language. Derived as `english_first_language %in% c("Yes","是")` **or** `french_first_language == "Oui"` **or** `first_language %in% c("Dansk","Deutsch","हिंदी","한국어","Русский")`. |

> **Note on locale.** Because the harmonized columns and several analyses compare
> non-Latin string literals (Hindi, Korean, Cyrillic, `是`), the R scripts must be
> run in a UTF-8 locale (the default in RStudio; for `Rscript` set, e.g.,
> `LC_ALL=en_US.UTF-8`). In a `C` locale these comparisons silently fail to match.
