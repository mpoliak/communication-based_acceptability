"""
Score the English experiment's stimuli under 6 language models x 3 unigram
models: the 18 settings of SLOR and MORCELA in SI section 4
(supporting_information/section4_slor_morcela_settings.R).

Language models (log_p_sentence):
    goldfish  goldfish-models/eng_latn_1000mb  (the main text's model)
    gpt2      gpt2
    gpt2-xl   gpt2-xl
    kn2/3/4   interpolated Kneser-Ney 2-, 3- and 4-gram models, discount 0.75,
              trained on about 20 million words of WikiText-103 (kn_ngram.py)

Unigram models (log_p_unigram):
    wordfreq       corpus word frequencies from the `wordfreq` package
    corpus         the word's relative frequency in the same WikiText-103 text
                   that trains the n-gram models (an unseen word counts as half
                   an occurrence)
    model_initial  log2 p(word | start of sentence) from the language model
                   itself, each word scored on its own (the main text's)

All logs are base 2. A word is a whitespace-delimited token of the normalized
sentence (perceived_sentence: lowercased, punctuation and symbols removed), so
the length matches n_words in ../en/stimuli_data.csv; the script checks this.

Reads ../en/stimuli_data.csv, corpus/unigram_counts.csv (from build_corpus.py)
and models/kn2_train.pkl, kn3_train.pkl, kn4_train.pkl (from
train_kn_on_corpus.py). Writes ../slor_morcela_scores_english_settings.csv
(one row per sentence x language model x unigram model).

Run order, from this folder, in the goldfish conda environment
(../goldfish_environment.yml):
    python build_corpus.py <path to WikiText-103-raw-v1 wiki.train.raw>
    python train_kn_on_corpus.py
    python get_slor_morcela_scores_english_settings.py
The three neural models are loaded from the local Hugging Face cache only
(local_files_only in neural_backends.py), so download goldfish-models/
eng_latn_1000mb, gpt2 and gpt2-xl first. The 4-gram model alone is a 983 MB
file: check free memory before running. The corpus/ and models/ folders are
not in the repository (.gitignore).
"""
import csv
import math
import os
import warnings

# torch 2.2.2 was built against NumPy 1.x, so importing it next to NumPy 2 warns
# "Failed to initialize NumPy". Nothing here passes tensors to NumPy, so only
# that warning is hidden; every other warning is shown.
warnings.filterwarnings("ignore", message="Failed to initialize NumPy")

import pandas as pd
from wordfreq import word_frequency

# The helpers live next to this script
from neural_backends import load_backend
from kn_ngram import BOS, KNLanguageModel

# Paths below are relative to this folder, data/probability_based_models_scoring/
os.chdir(os.path.dirname(os.path.abspath(__file__)))

NEURAL = ["goldfish", "gpt2", "gpt2-xl"]
KN_ORDERS = [2, 3, 4]
UNIGRAMS = ["wordfreq", "corpus", "model_initial"]

# p_u floor for words wordfreq has never seen (the number of such words is
# printed below).
FLOOR_P = 1e-9


def load_corpus_unigram():
    path = os.path.join("corpus", "unigram_counts.csv")
    counts, total = {}, 0
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            c = int(row["count"])
            counts[row["word"]] = c
            total += c
    # Unseen words get half a count, the standard mild floor.
    return counts, total


def main():
    stim = pd.read_csv(os.path.join("..", "en", "stimuli_data.csv"))
    stim = stim[stim["perceived_sentence"].notna()].reset_index(drop=True)
    sentences = [str(s) for s in stim["perceived_sentence"]]
    print(f"{len(sentences)} English sentences", flush=True)

    uni_counts, uni_total = load_corpus_unigram()
    print(f"corpus unigram: {len(uni_counts)} types, {uni_total} tokens", flush=True)

    # Cache the two LM-independent unigram tables over the stimulus vocabulary.
    vocab = sorted({w for s in sentences for w in s.split()})
    floors = {"wordfreq": 0, "corpus": 0}

    log_pu_wordfreq = {}
    for w in vocab:
        f = word_frequency(w, "en")
        if f <= 0.0:
            f, floors["wordfreq"] = FLOOR_P, floors["wordfreq"] + 1
        log_pu_wordfreq[w] = math.log2(f)

    log_pu_corpus = {}
    for w in vocab:
        c = uni_counts.get(w, 0)
        if c == 0:
            floors["corpus"] += 1
            p = 0.5 / uni_total
        else:
            p = c / uni_total
        log_pu_corpus[w] = math.log2(p)

    print(f"stimulus vocab {len(vocab)} types; "
          f"unseen by wordfreq {floors['wordfreq']}, by corpus {floors['corpus']}",
          flush=True)

    rows = []

    def emit(lm_name, per_sentence_pm, init_log2p):
        """per_sentence_pm: list of per-word log2 p_m lists, aligned to `sentences`."""
        for i, per_word in enumerate(per_sentence_pm):
            words = sentences[i].split()
            n = len(words)
            assert len(per_word) == n, (lm_name, i, len(per_word), n)
            log_pm = math.fsum(per_word)
            pu_tables = {
                "wordfreq": math.fsum(log_pu_wordfreq[w] for w in words),
                "corpus": math.fsum(log_pu_corpus[w] for w in words),
                "model_initial": math.fsum(init_log2p(w) for w in words),
            }
            # The length must match the repository's n_words
            assert n == stim["n_words"][i], (lm_name, i, n, stim["n_words"][i])
            for uname in UNIGRAMS:
                rows.append({
                    "language": "en",
                    "item": stim["item"][i], "condition": stim["condition"][i],
                    "language_model": lm_name, "unigram_model": uname,
                    "log_p_sentence": log_pm, "log_p_unigram": pu_tables[uname],
                })

    # ---- Neural backends -----------------------------------------------------
    for name in NEURAL:
        print(f"scoring with {name} ...", flush=True)
        be = load_backend(name)
        scores = [be.score_sentence(s) for s in sentences]
        cache = {}

        def init_log2p(w, _be=be, _c=cache):
            if w not in _c:
                _c[w] = _be.initial_word_log2p(w)
            return _c[w]

        emit(name, scores, init_log2p)
        del be
        print(f"  done {name}", flush=True)

    # ---- Kneser-Ney n-grams --------------------------------------------------
    for order in KN_ORDERS:
        name = f"kn{order}"
        path = os.path.join("models", f"kn{order}_train.pkl")
        print(f"scoring with {name} (loading {os.path.basename(path)}) ...", flush=True)
        lm = KNLanguageModel.load(path)
        scores = [lm.score_sentence(s) for s in sentences]
        cache = {}

        def init_log2p(w, _lm=lm, _o=order, _c=cache):
            if w not in _c:
                _c[w] = math.log2(_lm.prob(w, [BOS] * (_o - 1)))
            return _c[w]

        emit(name, scores, init_log2p)
        del lm
        print(f"  done {name}", flush=True)

    out = pd.DataFrame(rows)
    dest = os.path.join("..", "slor_morcela_scores_english_settings.csv")
    out.to_csv(dest, index=False)
    print(f"\nwrote data/slor_morcela_scores_english_settings.csv ({len(out)} rows = "
          f"{out['language_model'].nunique()} LMs x {out['unigram_model'].nunique()} unigrams x "
          f"{len(sentences)} sentences)")


if __name__ == "__main__":
    main()
