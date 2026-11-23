"""
Train the interpolated Kneser-Ney language models (orders 2, 3 and 4, discount
0.75; kn_ngram.py) of SI section 4 on corpus/train.txt (step 2 of 3; see
get_slor_morcela_scores_english_settings.py for the run order):
    python train_kn_on_corpus.py

Every 200th sentence is held out (spread over the whole corpus rather than one
contiguous topic region); the models are trained on the rest. For each order
the script prints the held-out perplexity and checks that the probabilities
sum to 1 in a few contexts. It does this twice: keeping every word type
(min_count = 1, the models used) and, for comparison only, folding singletons
into <unk> (min_count = 2, not saved).

Writes models/kn2_train.pkl, kn3_train.pkl and kn4_train.pkl in this folder
(about 74 MB, 361 MB and 983 MB; git-ignored). Check free memory first.
"""

from __future__ import annotations

import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from kn_ngram import BOS, KNLanguageModel  # noqa: E402

# Paths below are relative to this folder, data/probability_based_models_scoring/
os.chdir(os.path.dirname(os.path.abspath(__file__)))
CORPUS = os.path.join("corpus", "train.txt")
MODEL_DIR = "models"
HELDOUT_EVERY = 200


def load_split():
    train, held = [], []
    with open(CORPUS, "r", encoding="utf-8") as fh:
        for i, ln in enumerate(fh):
            if not ln.split():
                continue
            (held if i % HELDOUT_EVERY == 0 else train).append(ln.rstrip("\n"))
    return train, held


def norm_contexts(order, held):
    """A few real contexts (plus the sentence start and an OOV one) to sum over."""
    ctxs = [[BOS] * (order - 1)]
    for s in held[:2000]:
        t = s.split()
        if len(t) >= order + 2:
            ctxs.append(t[1 : order])
            if len(ctxs) >= 5:
                break
    ctxs.append(["zzzqqqfloop"] * (order - 1))
    return ctxs


def main():
    os.makedirs(MODEL_DIR, exist_ok=True)
    t0 = time.time()
    train, held = load_split()
    print(f"train: {len(train)} sentences | held-out: {len(held)} sentences "
          f"({sum(len(s.split()) for s in held)} words) | load {time.time()-t0:.1f}s",
          flush=True)

    for min_count in (1, 2):
        print(f"\n{'='*70}\nmin_count = {min_count} "
              f"({'keep every observed type' if min_count == 1 else 'fold singletons into <unk>'})"
              f"\n{'='*70}", flush=True)
        for order in (2, 3, 4):
            t = time.time()
            lm = KNLanguageModel(order=order).fit(train, min_count=min_count)
            t_train = time.time() - t

            t = time.time()
            ppl, nw = lm.perplexity(held)
            t_eval = time.time() - t

            print(f"\norder {order}: V={lm.V} tokens={lm.n_train_tokens} "
                  f"| train {t_train:.1f}s, eval {t_eval:.1f}s", flush=True)
            print(f"  held-out perplexity = {ppl:.3f}  ({nw} words, base 2, no EOS)",
                  flush=True)

            worst = 0.0
            for ctx in norm_contexts(order, held):
                s = lm.distribution_sum(ctx)
                worst = max(worst, abs(s - 1.0))
                print(f"    sum_w P(w|{' '.join(ctx)!r}) = {s:.15f}", flush=True)
            print(f"  max |sum - 1| = {worst:.3e}", flush=True)

            if min_count == 1:
                p = os.path.join(MODEL_DIR, f"kn{order}_train.pkl")
                try:
                    lm.save(p)
                    print(f"  saved -> {p} ({os.path.getsize(p)/1e6:.0f} MB)", flush=True)
                except Exception as e:  # noqa: BLE001
                    print(f"  save skipped: {e}", flush=True)
            del lm

    print(f"\ntotal {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
