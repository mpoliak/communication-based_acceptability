"""
Interpolated Kneser-Ney n-gram language model (pure Python), for the n-gram
settings of SI section 4. Imported by train_kn_on_corpus.py (which trains and
saves the models) and get_slor_morcela_scores_english_settings.py (which scores
the stimuli with them). It writes nothing itself, except a model file when run
from the command line with --save (see CLI below).

Reproduces the model class used by Lau, Clark & Lappin: "an N-gram model with
Kneser-Ney interpolation", tested at orders 2, 3 and 4.

SHARED SCORING CONVENTIONS (must match the neural scorers exactly)
-----------------------------------------------------------------
* Sentences arrive already normalized: lowercase, punctuation stripped, single
  spaces. A "word" is a whitespace-delimited token; |x| = number of words.
* ``score_sentence`` returns a list of PER-WORD **log2** probabilities of length
  exactly |x|. Never natural log.
* Left context is padded with ``<s>`` as needed.
* NO end-of-sentence token is scored -- the neural scorers do not score one, so
  this model must not either. ``</s>`` is never added to the training data.
* Out-of-vocabulary words map to ``<unk>``, which is always a member of the
  vocabulary.

MODEL SPECIFICATION (standard interpolated KN, single fixed discount D)
----------------------------------------------------------------------
Highest order n (raw counts):

    P_n(w | h) = max(c(h,w) - D, 0) / c(h)
                 + (D * N1+(h,*) / c(h)) * P_{n-1}(w | h')

Lower orders k < n (CONTINUATION counts):

    P_k(w | h) = max(N1+(*,h,w) - D, 0) / N1+(*,h,*)
                 + (D * N1+(h,*) / N1+(*,h,*)) * P_{k-1}(w | h')

Unigram base:

    P_1(w) = max(N1+(*,w) - D, 0) / N1+(*,*)
             + (D * V_context / N1+(*,*)) * (1 / V)

where N1+(*,w) = number of distinct word types preceding w, N1+(*,*) = number of
distinct bigram types, V_context = number of word types with a nonzero
continuation count, and V = vocabulary size (including ``<unk>``, excluding the
padding symbol ``<s>``, which is never predicted).

D = 0.75.

Two details that are easy to get wrong and that this implementation is careful
about, because they are what make the distribution normalize:

1. ``c(h)`` is the number of times h occurs as the *prefix of an observed
   n-gram*, i.e. sum_w c(h,w) -- NOT the standalone corpus frequency of h.
   With ``<s>`` padding and no ``</s>`` those two numbers differ.
2. ``N1+(h,*)`` in the interpolation weight is the number of distinct w with a
   nonzero count *in the table used at that level* (raw counts at order n,
   continuation counts below). Mixing the two breaks normalization.

Under these definitions sum_w P_k(w | h) = 1 exactly for every level and every
context; see ``test_kn_ngram.py``.

Usage
-----
    from kn_ngram import KNLanguageModel, train_kn

    lm = train_kn("corpus/train.txt", order=3)
    lm.score_sentence("the dog barked")   # -> [log2 p1, log2 p2, log2 p3]

CLI
---
    python kn_ngram.py --train CORPUS --order 3 --save models/kn3.pkl
    python kn_ngram.py --train CORPUS --order 2 --heldout-frac 0.05
    python kn_ngram.py --load models/kn3.pkl --score-file sents.txt
"""

from __future__ import annotations

import argparse
import math
import os
import pickle
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Sequence, Tuple

BOS = "<s>"
UNK = "<unk>"
DEFAULT_D = 0.75

_LOG2 = math.log(2.0)


class KNLanguageModel:
    """Interpolated Kneser-Ney n-gram LM with a fixed discount D."""

    def __init__(self, order: int, D: float = DEFAULT_D):
        if order < 1:
            raise ValueError("order must be >= 1")
        if not (0.0 < D < 1.0):
            raise ValueError("D must be in (0, 1)")
        self.order = int(order)
        self.D = float(D)

        # vocabulary
        self.word2id: Dict[str, int] = {}
        self.id2word: List[str] = []
        self.bos_id: int = -1
        self.unk_id: int = -1
        self.V: int = 0  # size of the predictable vocabulary (incl <unk>, excl <s>)

        # tables[k] : dict mapping k-gram (tuple of ids) -> count
        #   k == order : raw counts;  k < order : continuation counts
        self.tables: List[Dict[Tuple[int, ...], int]] = []
        # ctx_total[k][h]    = sum_w tables[k][h + (w,)]
        # ctx_distinct[k][h] = |{w : tables[k][h + (w,)] > 0}|
        self.ctx_total: List[Dict[Tuple[int, ...], int]] = []
        self.ctx_distinct: List[Dict[Tuple[int, ...], int]] = []

        self.n_train_sentences = 0
        self.n_train_tokens = 0
        self.min_count = 1

    # ------------------------------------------------------------------ vocab

    def _add_word(self, w: str) -> int:
        i = self.word2id.get(w)
        if i is None:
            i = len(self.id2word)
            self.word2id[w] = i
            self.id2word.append(w)
        return i

    def _encode_train(self, tokens: Sequence[str]) -> List[int]:
        out = []
        for w in tokens:
            if w == BOS:  # a literal "<s>" in the data must not alias the pad symbol
                w = UNK
            out.append(self.word2id.get(w, self.unk_id))
        return out

    def encode(self, tokens: Sequence[str]) -> List[int]:
        """Map surface words to ids, sending OOV (and literal <s>) to <unk>."""
        unk = self.unk_id
        w2i = self.word2id
        return [unk if w == BOS else w2i.get(w, unk) for w in tokens]

    # --------------------------------------------------------------- training

    def fit(self, sentences: Iterable[str], min_count: int = 1) -> "KNLanguageModel":
        """Train from an iterable of whitespace-tokenizable sentence strings.

        ``min_count``: training words occurring fewer than this many times are
        folded into <unk>.  Default 1 keeps every observed type (in which case
        <unk> receives only the smoothing floor mass).
        """
        self.min_count = int(min_count)
        n = self.order

        tokenized: List[List[str]] = []
        freq: Counter = Counter()
        for line in sentences:
            toks = line.split()
            if not toks:
                continue
            tokenized.append(toks)
            freq.update(toks)

        # ---- vocabulary: <s> first (id 0, never predicted), then <unk>, then types
        self.word2id, self.id2word = {}, []
        self.bos_id = self._add_word(BOS)
        self.unk_id = self._add_word(UNK)
        for w, c in sorted(freq.items()):
            if w in (BOS, UNK):
                continue
            if c >= self.min_count:
                self._add_word(w)
        # V counts every id that can be predicted: everything except <s>
        self.V = len(self.id2word) - 1

        self.n_train_sentences = len(tokenized)
        self.n_train_tokens = sum(len(t) for t in tokenized)

        # ---- highest-order raw counts
        tables: List[Dict[Tuple[int, ...], int]] = [dict() for _ in range(n + 1)]
        top: Dict[Tuple[int, ...], int] = defaultdict(int)
        bos = self.bos_id
        pad = (bos,) * (n - 1)
        for toks in tokenized:
            ids = pad + tuple(self._encode_train(toks))
            for i in range(n - 1, len(ids)):
                top[ids[i - n + 1 : i + 1]] += 1
        tables[n] = dict(top)

        # ---- lower-order CONTINUATION counts, derived from the type set above
        for k in range(n - 1, 0, -1):
            cont: Dict[Tuple[int, ...], int] = defaultdict(int)
            for gram in tables[k + 1]:
                cont[gram[1:]] += 1
            tables[k] = dict(cont)

        # ---- context totals and distinct-continuation counts per level
        ctx_total: List[Dict[Tuple[int, ...], int]] = [dict() for _ in range(n + 1)]
        ctx_distinct: List[Dict[Tuple[int, ...], int]] = [dict() for _ in range(n + 1)]
        for k in range(1, n + 1):
            tot: Dict[Tuple[int, ...], int] = defaultdict(int)
            dis: Dict[Tuple[int, ...], int] = defaultdict(int)
            for gram, c in tables[k].items():
                h = gram[:-1]
                tot[h] += c
                dis[h] += 1
            ctx_total[k] = dict(tot)
            ctx_distinct[k] = dict(dis)

        self.tables = tables
        self.ctx_total = ctx_total
        self.ctx_distinct = ctx_distinct
        return self

    @classmethod
    def train_from_file(
        cls, path: str, order: int, D: float = DEFAULT_D, min_count: int = 1
    ) -> "KNLanguageModel":
        with open(path, "r", encoding="utf-8") as fh:
            return cls(order=order, D=D).fit(fh, min_count=min_count)

    # --------------------------------------------------------------- scoring

    def _p(self, level: int, ctx: Tuple[int, ...], w: int) -> float:
        """P_level(w | ctx), ctx being a tuple of (level-1) ids."""
        tot = self.ctx_total[level].get(ctx, 0)
        if tot == 0:
            # unseen context: back off wholly to the next lower order, which
            # keeps the distribution normalized
            if level == 1:
                return 1.0 / self.V
            return self._p(level - 1, ctx[1:], w)
        c = self.tables[level].get(ctx + (w,), 0)
        nd = self.ctx_distinct[level].get(ctx, 0)
        lam = self.D * nd / tot
        lower = (1.0 / self.V) if level == 1 else self._p(level - 1, ctx[1:], w)
        return max(c - self.D, 0.0) / tot + lam * lower

    def prob(self, word: str, context: Sequence[str]) -> float:
        """Linear-domain P(word | context); context is a sequence of words."""
        n = self.order
        ctx_ids = tuple(([self.bos_id] * (n - 1) + self.encode(list(context)))[-(n - 1) :]) if n > 1 else ()
        w = self.encode([word])[0]
        return self._p(n, ctx_ids, w)

    def score_sentence(self, sentence: str) -> List[float]:
        """Per-word log2 probabilities; the returned list has length |x|."""
        toks = sentence.split()
        if not toks:
            return []
        n = self.order
        ids = [self.bos_id] * (n - 1) + self.encode(toks)
        out: List[float] = []
        for i in range(n - 1, len(ids)):
            p = self._p(n, tuple(ids[i - n + 1 : i]), ids[i])
            out.append(math.log(p) / _LOG2)
        assert len(out) == len(toks)
        return out

    def logprob(self, sentence: str) -> float:
        """Total log2 probability of the sentence."""
        return float(sum(self.score_sentence(sentence)))

    def perplexity(self, sentences: Iterable[str]) -> Tuple[float, int]:
        """(perplexity, n_words) over the given sentences, base-2, no EOS."""
        total, n_words = 0.0, 0
        for s in sentences:
            sc = self.score_sentence(s)
            total += sum(sc)
            n_words += len(sc)
        if n_words == 0:
            return float("nan"), 0
        return 2.0 ** (-total / n_words), n_words

    # ------------------------------------------------------- diagnostics / io

    def predictable_ids(self) -> List[int]:
        """Every id the model can predict (all of V; excludes <s>)."""
        return [i for i in range(len(self.id2word)) if i != self.bos_id]

    def distribution_sum(self, context: Sequence[str]) -> float:
        """Sum of P(w | context) over the whole vocabulary. Must be 1.0."""
        n = self.order
        ctx_ids = tuple(([self.bos_id] * (n - 1) + self.encode(list(context)))[-(n - 1) :]) if n > 1 else ()
        return math.fsum(self._p(n, ctx_ids, w) for w in self.predictable_ids())

    def save(self, path: str) -> None:
        os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
        with open(path, "wb") as fh:
            pickle.dump(self, fh, protocol=pickle.HIGHEST_PROTOCOL)

    @staticmethod
    def load(path: str) -> "KNLanguageModel":
        with open(path, "rb") as fh:
            return pickle.load(fh)

    def __repr__(self) -> str:
        return (
            f"KNLanguageModel(order={self.order}, D={self.D}, V={self.V}, "
            f"sents={self.n_train_sentences}, tokens={self.n_train_tokens})"
        )


def train_kn(
    corpus_path: str, order: int, D: float = DEFAULT_D, min_count: int = 1
) -> KNLanguageModel:
    """Convenience entry point: train an order-`order` KN LM from a text file."""
    return KNLanguageModel.train_from_file(corpus_path, order=order, D=D, min_count=min_count)


# --------------------------------------------------------------------- CLI

def _main() -> None:
    ap = argparse.ArgumentParser(description="Interpolated Kneser-Ney n-gram LM")
    ap.add_argument("--train", help="path to a training corpus (one sentence per line)")
    ap.add_argument("--order", type=int, default=3, help="n-gram order (2, 3 or 4)")
    ap.add_argument("--discount", type=float, default=DEFAULT_D)
    ap.add_argument("--min-count", type=int, default=1)
    ap.add_argument("--heldout-frac", type=float, default=0.0,
                    help="hold out this fraction of lines (from the end) for perplexity")
    ap.add_argument("--save", help="pickle the trained model here")
    ap.add_argument("--load", help="load a pickled model instead of training")
    ap.add_argument("--score-file", help="score each line of this file")
    args = ap.parse_args()

    if args.load:
        lm = KNLanguageModel.load(args.load)
    elif args.train:
        with open(args.train, "r", encoding="utf-8") as fh:
            lines = [ln.rstrip("\n") for ln in fh]
        lines = [ln for ln in lines if ln.split()]
        if args.heldout_frac > 0:
            k = max(1, int(len(lines) * args.heldout_frac))
            train_lines, held = lines[:-k], lines[-k:]
        else:
            train_lines, held = lines, []
        lm = KNLanguageModel(order=args.order, D=args.discount).fit(
            train_lines, min_count=args.min_count
        )
        print(lm)
        if held:
            ppl, nw = lm.perplexity(held)
            print(f"held-out: {len(held)} sentences, {nw} words, perplexity = {ppl:.3f}")
    else:
        ap.error("one of --train or --load is required")

    if args.save:
        lm.save(args.save)
        print(f"saved -> {args.save}")

    if args.score_file:
        with open(args.score_file, "r", encoding="utf-8") as fh:
            for ln in fh:
                s = ln.strip()
                if not s:
                    continue
                sc = lm.score_sentence(s)
                print(f"{sum(sc):.6f}\t{len(sc)}\t{s}")


if __name__ == "__main__":
    _main()
