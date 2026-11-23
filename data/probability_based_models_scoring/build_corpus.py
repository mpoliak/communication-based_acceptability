#!/usr/bin/env python
"""
Build the training text of the Kneser-Ney n-gram models and the "corpus"
unigram model of SI section 4 (step 1 of 3; see
get_slor_morcela_scores_english_settings.py for the run order).

Source: WikiText-103-raw-v1 (wiki.train.raw), the raw release of English
Wikipedia text. Pass its path as the only argument:
    python build_corpus.py <path to wiki.train.raw>

Each line is split into sentences, and each sentence is normalized like the
stimuli (lowercase, punctuation and symbols removed with the same \\p{P}\\p{S}
rule, whitespace collapsed); sentences shorter than 3 words are dropped, and
the text stops at 20 million words. A "word" is a whitespace-delimited token.

Writes corpus/train.txt (one sentence per line, about 118 MB) and
corpus/unigram_counts.csv (word, count), in this folder. Both are git-ignored.
"""
import csv
import os
import sys
from collections import Counter

import regex

RAW = os.path.abspath(sys.argv[1])

# Paths below are relative to this folder, data/probability_based_models_scoring/
os.chdir(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = "corpus"
TRAIN_PATH = os.path.join(OUT_DIR, "train.txt")
UNIGRAM_PATH = os.path.join(OUT_DIR, "unigram_counts.csv")

TOKEN_CAP = 20_000_000
MIN_WORDS = 3

# --- compiled patterns -------------------------------------------------------
# Wikipedia heading lines, e.g. " = = Gameplay = = "
HEADING = regex.compile(r"^\s*=.*=\s*$")
# WikiText's escaped intra-word punctuation: "role @-@ playing" -> "role-playing",
# "1 @.@ 5" -> "1.5", "12 @,@ 000" -> "12,000". Undoing these BEFORE punctuation
# stripping is what makes "role-playing" collapse to the single token
# "roleplaying", exactly as it would in a normalized stimulus sentence.
DETOK = regex.compile(r"\s@([-.,])@\s")
# WikiText also splits a space before clitics: "game 's", "don 't", "they 're".
# Left alone, the shared punctuation rule turns those into the stray tokens
# "game s" / "don t", while a normalized stimulus yields "games" / "dont".
# Rejoin the closed clitic set (only these, so quoted phrases are untouched).
CLITIC = regex.compile(r"([\p{L}\p{N}])\s*['’]\s*(s|t|re|ve|ll|d|m)\b")
# Simple sentence split on sentence-final punctuation followed by whitespace.
SENT_SPLIT = regex.compile(r"(?<=[.!?])\s+")
# The study's punctuation/symbol stripping rule.
PUNCT = regex.compile(r"[\p{P}\p{S}]+")

os.makedirs(OUT_DIR, exist_ok=True)

n_sentences = 0
n_tokens = 0
counts = Counter()
lines_read = 0
hit_cap = False

with open(RAW, "r", encoding="utf-8") as fin, \
     open(TRAIN_PATH, "w", encoding="utf-8") as fout:

    for line in fin:
        if hit_cap:
            break
        lines_read += 1

        stripped = line.strip()
        if not stripped:                 # blank line
            continue
        if HEADING.match(line):          # "= = Heading = ="
            continue

        # undo WikiText's @-@ / @.@ / @,@ escapes
        line = DETOK.sub(r"\1", line)

        for raw_sent in SENT_SPLIT.split(line):
            s = CLITIC.sub(r"\1'\2", raw_sent.lower())
            s = PUNCT.sub("", s)
            words = s.split()
            if len(words) < MIN_WORDS:
                continue

            if n_tokens + len(words) > TOKEN_CAP:
                hit_cap = True
                break

            fout.write(" ".join(words) + "\n")
            n_sentences += 1
            n_tokens += len(words)
            counts.update(words)

# --- unigram counts ----------------------------------------------------------
with open(UNIGRAM_PATH, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f)
    w.writerow(["word", "count"])
    for word, c in counts.most_common():
        w.writerow([word, c])

vocab_size = len(counts)

print(f"source file       : {RAW}")
print(f"lines read        : {lines_read:,}")
print(f"total sentences   : {n_sentences:,}")
print(f"total word tokens : {n_tokens:,}")
print(f"vocabulary size   : {vocab_size:,}")
print()
print("15 most frequent words:")
print(f"{'rank':>4}  {'word':<12} {'count':>10}")
for i, (word, c) in enumerate(counts.most_common(15), 1):
    print(f"{i:>4}  {word:<12} {c:>10,}")

print()
print(f"train    -> {TRAIN_PATH}")
print(f"unigrams -> {UNIGRAM_PATH}")
