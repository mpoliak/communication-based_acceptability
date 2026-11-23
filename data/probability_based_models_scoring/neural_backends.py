"""
Multi-backend neural sentence scorer for the SLOR / Lau-Clark-Lappin measures:
the Goldfish, GPT-2 and GPT-2 XL settings of SI section 4. Imported by
get_slor_morcela_scores_english_settings.py; it writes nothing. Run directly, it
only runs the self-validation below and prints the result.

Exposes
    load_backend(name) -> backend
with

    backend.name                       -> str
    backend.score_sentence(sentence)   -> list[float], length |x|
                                          per-word log2 p(word_i | preceding words)
    backend.initial_word_log2p(word)   -> float
                                          log2 p(word | BOS), for the sentence-initial
                                          "unigram" variant p_u

Backends
    "goldfish"  goldfish-models/eng_latn_1000mb  (GPT2LMHeadModel + ALBERT
                SentencePiece tokenizer; BOS marker "[CLS]")
    "gpt2"      gpt2                             (BOS marker "<|endoftext|>")
    "gpt2-xl"   gpt2-xl                          (BOS marker "<|endoftext|>")

SHARED CONVENTIONS (must match the n-gram scorers exactly)
    * Sentences arrive already normalized: lowercase, no punctuation, single spaces.
    * A "word" is a whitespace-delimited token; |x| = len(sentence.split()).
      Bare split() -- whitespace runs, no empty words -- is what the n-gram
      scorers use, so the two families always agree on |x|.
    * Output is ALWAYS log2, never natural log.
    * The BOS marker is prepended as context only and is never scored.
    * No end-of-sentence token is scored (the n-gram side must not score one either).
    * Every one of the |x| words receives exactly one score; subword log2 probs are
      summed WITHIN each whitespace word.

WORD ALIGNMENT -- the one thing that has to be exactly right
    Both tokenizer families here return `offset_mapping` spans that INCLUDE the
    leading space of a word:

        goldfish/SentencePiece   "they have"  -> "_they"(0,4)  "_have"(4,9)
        gpt2/byte-level BPE      "they have"  -> "they"(0,4)   "Ġhave"(4,9)

    so a subword's start offset can point at the SPACE in front of its word rather
    than at the word's first letter. Assigning by raw start offset would therefore
    push every non-initial word's first subword onto the PREVIOUS word.

    The rule used here is explicit rather than incidental: for each subword we take
    the first NON-SPACE character position inside its span, and assign the subword
    to the whitespace word that owns that character (via a precomputed
    char-index -> word-index map). A subword whose span is entirely whitespace
    (possible in principle for byte-level BPE) is attached to the word that follows
    it. This is exact for both leading-space and no-leading-space conventions, and
    `alignment_report()` prints the resulting grouping so it can be inspected.

Loading note
    The goldfish checkpoint on the `main` revision ships only `pytorch_model.bin`,
    and transformers >= 4.56 refuses `torch.load` on torch < 2.6 (CVE-2025-32434).
    The local env is torch 2.2.2, so `_load_model` falls back to instantiating the
    model from its config and loading the state dict with `torch.load` directly.
    The fallback asserts that no real parameter is missing; the only unexpected keys
    accepted are the legacy `attn.bias` / `attn.masked_bias` causal-mask buffers,
    which modern transformers builds as non-persistent buffers.

Run this file directly to execute the self-validation:
    python neural_backends.py      (in the goldfish conda environment)
"""

import math
import os
import warnings

# torch 2.2.2 was built against NumPy 1.x, so importing it next to NumPy 2 warns
# "Failed to initialize NumPy". Nothing here passes tensors to NumPy, so only
# that warning is hidden; every other warning is shown.
warnings.filterwarnings("ignore", message="Failed to initialize NumPy")
os.environ.setdefault("TRANSFORMERS_VERBOSITY", "error")
os.environ.setdefault("TOKENIZERS_PARALLELISM", "false")

import torch
from transformers import AutoConfig, AutoModelForCausalLM, AutoTokenizer
from transformers.utils import cached_file

LN2 = math.log(2.0)

# name -> (hf model id, BOS marker token string)
BACKEND_SPECS = {
    "goldfish": ("goldfish-models/eng_latn_1000mb", "[CLS]"),
    "gpt2": ("gpt2", None),        # None -> use tokenizer.bos_token
    "gpt2-xl": ("gpt2-xl", None),
}

# Buffers that old GPT-2 checkpoints carry but current transformers rebuilds itself.
_IGNORABLE_UNEXPECTED = ("attn.bias", "attn.masked_bias")


# --------------------------------------------------------------------------- #
# model loading
# --------------------------------------------------------------------------- #
def _load_model(model_id):
    """from_pretrained, with a torch.load fallback for .bin-only checkpoints."""
    try:
        return AutoModelForCausalLM.from_pretrained(
            model_id, local_files_only=True, low_cpu_mem_usage=True
        ).eval()
    except Exception:
        pass

    config = AutoConfig.from_pretrained(model_id, local_files_only=True)
    model = AutoModelForCausalLM.from_config(config)
    path = cached_file(model_id, "pytorch_model.bin", local_files_only=True)
    state = torch.load(path, map_location="cpu", weights_only=True)
    result = model.load_state_dict(state, strict=False)
    if result.missing_keys:
        raise RuntimeError(f"{model_id}: missing weights {result.missing_keys}")
    stray = [k for k in result.unexpected_keys
             if not any(k.endswith(s) for s in _IGNORABLE_UNEXPECTED)]
    if stray:
        raise RuntimeError(f"{model_id}: unexpected weights {stray}")
    model.tie_weights()
    return model.eval()


# --------------------------------------------------------------------------- #
# word / character bookkeeping
# --------------------------------------------------------------------------- #
def word_spans(sentence):
    """(start, end) character span of each whitespace-delimited word, in order.

    Deliberately matches ``str.split()`` semantics: ANY run of whitespace is a
    separator, and leading/trailing whitespace yields no empty word. This must
    agree with the n-gram scorers, which all tokenize with bare
    ``sentence.split()`` (kn_ngram.py, kenlm_ngram.py) -- if the two families
    disagreed about |x| the SLOR denominators would silently desync.

    Splitting on the literal " " instead would be wrong in two silent ways:
    a trailing space produces a phantom empty word (|x| + 1 scores), and a tab
    or newline separator fails to split at all (|x| - 1 scores).
    """
    spans = []
    i, n = 0, len(sentence)
    while i < n:
        while i < n and sentence[i].isspace():
            i += 1
        if i >= n:
            break
        start = i
        while i < n and not sentence[i].isspace():
            i += 1
        spans.append((start, i))
    return spans


def _char_to_word(sentence, spans):
    """char index -> word index, or -1 for characters that belong to no word."""
    owner = [-1] * len(sentence)
    for wi, (s, e) in enumerate(spans):
        for i in range(s, e):
            owner[i] = wi
    return owner


def _owning_word(span, owner, n_words):
    """Word index a subword belongs to, ignoring any leading whitespace it covers.

    Returns None for a span that covers no character at all (empty span), and for a
    whitespace-only span returns the word that starts after it (a byte-level BPE
    token made only of spaces belongs with the word it precedes).
    """
    s, e = span
    if s >= e:
        return None
    for i in range(s, min(e, len(owner))):
        if owner[i] >= 0:
            return owner[i]
    for i in range(min(e, len(owner)), len(owner)):   # whitespace-only span
        if owner[i] >= 0:
            return owner[i]
    return n_words - 1


# --------------------------------------------------------------------------- #
# backend
# --------------------------------------------------------------------------- #
class HFCausalBackend:
    """Autoregressive HF causal LM wrapped in the shared scoring convention."""

    def __init__(self, name, model_id, bos_marker=None):
        self.name = name
        self.model_id = model_id
        self.tokenizer = AutoTokenizer.from_pretrained(
            model_id, use_fast=True, local_files_only=True
        )
        if not self.tokenizer.is_fast:
            raise RuntimeError(f"{model_id}: need a fast tokenizer for offsets")
        self.model = _load_model(model_id)

        if bos_marker is None:
            bos_marker = self.tokenizer.bos_token
        if bos_marker is None:
            raise RuntimeError(f"{model_id}: no BOS marker available")
        self.bos_marker = bos_marker
        bos_id = self.tokenizer.convert_tokens_to_ids([bos_marker])[0]
        # NB: for GPT-2, bos_token == eos_token == unk_token == "<|endoftext|>", so
        # "is it the unk id" is not a usable check; round-trip the id instead.
        if bos_id is None or self.tokenizer.convert_ids_to_tokens(bos_id) != bos_marker:
            raise RuntimeError(f"{model_id}: BOS marker {bos_marker!r} is not a token")
        self.bos_id = int(bos_id)
        self._initial_cache = {}

    # -- low level ---------------------------------------------------------- #
    def _encode(self, text):
        enc = self.tokenizer(text, add_special_tokens=False,
                             return_offsets_mapping=True)
        return list(enc["input_ids"]), [tuple(o) for o in enc["offset_mapping"]]

    def _token_log2_probs(self, text):
        """log2 p(token | BOS + preceding) for every subword of `text`.

        The BOS marker is prepended as a token id, so the returned offsets index
        `text` itself. Returns (ids, tokens, log2probs, offsets), all length T,
        one entry per subword of `text` -- the BOS position is context only.
        """
        ids, offsets = self._encode(text)
        if not ids:
            return [], [], [], []
        full = torch.tensor([[self.bos_id] + ids])
        with torch.no_grad():
            logits = self.model(input_ids=full).logits
        log_probs = torch.log_softmax(logits.float(), dim=-1)
        targets = full[:, 1:]
        picked = log_probs[:, :-1, :].gather(2, targets.unsqueeze(-1)).squeeze(-1)[0]
        lps = [float(v) / LN2 for v in picked]           # nats -> log2
        toks = self.tokenizer.convert_ids_to_tokens(ids)
        return ids, toks, lps, offsets

    # -- public API --------------------------------------------------------- #
    def score_sentence(self, sentence):
        """Per-word log2 p(word_i | preceding words); length == number of words."""
        return [g["log2p"] for g in self.alignment_report(sentence)]

    def alignment_report(self, sentence):
        """Per-word records: {word, index, subwords, offsets, log2p}.

        Same computation as score_sentence, but keeps the subword grouping so the
        alignment can be printed and checked.
        """
        spans = word_spans(sentence)
        owner = _char_to_word(sentence, spans)
        _, toks, lps, offsets = self._token_log2_probs(sentence)

        groups = [{"index": i, "word": sentence[s:e], "subwords": [],
                   "offsets": [], "log2p": 0.0}
                  for i, (s, e) in enumerate(spans)]
        for tok, lp, off in zip(toks, lps, offsets):
            wi = _owning_word(off, owner, len(spans))
            if wi is None:                                # zero-width token
                continue
            groups[wi]["subwords"].append(tok)
            groups[wi]["offsets"].append(off)
            groups[wi]["log2p"] += lp

        empty = [g["word"] for g in groups if not g["subwords"]]
        if empty:
            raise RuntimeError(
                f"{self.name}: no subword assigned to word(s) {empty} in {sentence!r}")
        return groups

    def initial_word_log2p(self, word):
        """log2 p(word | BOS) -- the sentence-initial 'unigram' variant of p_u."""
        if word not in self._initial_cache:
            _, _, lps, offsets = self._token_log2_probs(word)
            self._initial_cache[word] = sum(
                lp for lp, (s, e) in zip(lps, offsets) if s != e)
        return self._initial_cache[word]

    def __repr__(self):
        return f"<{type(self).__name__} {self.name} ({self.model_id})>"


_LOADED = {}


def load_backend(name):
    """Return the (cached) backend named `name`."""
    if name not in BACKEND_SPECS:
        raise KeyError(f"unknown backend {name!r}; have {sorted(BACKEND_SPECS)}")
    if name not in _LOADED:
        model_id, bos = BACKEND_SPECS[name]
        _LOADED[name] = HFCausalBackend(name, model_id, bos)
    return _LOADED[name]


def available_backends():
    return list(BACKEND_SPECS)


# --------------------------------------------------------------------------- #
# validation
# --------------------------------------------------------------------------- #
def _reference_goldfish_scores(backend, sentence):
    """Re-implementation of the word alignment of an earlier Goldfish scorer
    (not in this repository).

    Kept here so the new code can be shown to reproduce the already-validated
    Goldfish numbers exactly.
    """
    cls = "[CLS]"
    prefixed = cls + " " + sentence
    n_prefix = len(cls) + 1
    tok = backend.tokenizer
    enc = tok(prefixed, add_special_tokens=False, return_offsets_mapping=True)
    ids = torch.tensor([enc["input_ids"]])
    with torch.no_grad():
        logits = backend.model(input_ids=ids).logits
    log_probs = torch.log_softmax(logits.float(), dim=-1)
    tok_lp = log_probs[:, :-1, :].gather(
        2, ids[:, 1:].unsqueeze(-1)).squeeze(-1)[0]
    lps = [float(v) / LN2 for v in tok_lp]
    offsets = [tuple(o) for o in enc["offset_mapping"][1:]]

    spans, cursor = [], n_prefix
    for w in prefixed[n_prefix:].split(" "):
        start = prefixed.index(w, cursor)
        spans.append((start, start + len(w)))
        cursor = start + len(w)

    per_word = [0.0] * len(spans)
    wi = 0
    for lp, (s, e) in zip(lps, offsets):
        if s == e:
            continue
        while wi < len(spans) - 1 and s >= spans[wi][1]:
            wi += 1
        per_word[wi] += lp
    return per_word


def _finite(values):
    return all(math.isfinite(v) for v in values)


def _validate():
    import json
    import random

    sentence = "they have a warm local atmosphere and taste great"
    words = sentence.split(" ")
    n = len(words)

    # The required sentence happens to be one subword per word on all three
    # models, which does not exercise the within-word summation at all. This
    # probe deliberately uses words that shatter into many subwords.
    probe = ("the lugubrious antidisestablishmentarianism unquestionably "
             "bewildered several onlookers")
    probe_words = probe.split(" ")

    scrambled_words = list(words)
    random.Random(20260825).shuffle(scrambled_words)
    scrambled = " ".join(scrambled_words)

    print("=" * 78)
    print("NEURAL BACKEND VALIDATION")
    print("=" * 78)
    print(f"sentence : {sentence!r}  (|x| = {n})")
    print(f"scrambled: {scrambled!r}  (|x| = {len(scrambled_words)})")
    print()

    ok = True
    working, counts, totals, notes = [], {}, {}, []

    for name in available_backends():
        print("-" * 78)
        print(f"BACKEND {name}")
        print("-" * 78)
        try:
            b = load_backend(name)
        except Exception as exc:
            ok = False
            notes.append(f"{name}: failed to load ({type(exc).__name__}: {exc})")
            print(f"  LOAD FAILED: {type(exc).__name__}: {exc}")
            continue

        print(f"  model_id   : {b.model_id}")
        print(f"  tokenizer  : {type(b.tokenizer).__name__}")
        print(f"  BOS marker : {b.bos_marker!r}  (id {b.bos_id}, context only, never scored)")

        report = b.alignment_report(sentence)
        per_word = [g["log2p"] for g in report]

        # raw subword stream, to show what the offsets actually look like
        _, toks, lps, offs = b._token_log2_probs(sentence)
        print(f"  raw subwords ({len(toks)}): "
              + " ".join(f"{t}{o}" for t, o in zip(toks, offs)))
        print("  subword grouping (offsets include the leading space -> "
              "assigned by first non-space char):")
        print(f"    {'#':>2}  {'word':<12} {'log2 p':>10}   subwords (offsets)")
        for g in report:
            sub = " ".join(f"{t}{o}" for t, o in zip(g["subwords"], g["offsets"]))
            print(f"    {g['index']:>2}  {g['word']:<12} {g['log2p']:>10.4f}   {sub}")

        assert len(per_word) == n, f"{name}: got {len(per_word)} scores, expected {n}"
        assert _finite(per_word), f"{name}: non-finite per-word score"
        # every subword accounted for exactly once
        n_assigned = sum(len(g["subwords"]) for g in report)
        assert n_assigned == len([o for o in offs if o[0] != o[1]]), \
            f"{name}: {n_assigned} subwords grouped but {len(offs)} produced"
        assert abs(sum(per_word) - sum(lps)) < 1e-6, \
            f"{name}: word sums do not equal the subword sum"
        print(f"  len(score_sentence) == {len(per_word)}  [OK, expected {n}]")
        print(f"  every subword assigned exactly once: {n_assigned}/{len(offs)}  [OK]")
        print(f"  sum(per-word) == sum(per-subword): "
              f"{sum(per_word):.6f} == {sum(lps):.6f}  [OK]")

        # goldfish must reproduce the already-validated reference implementation
        if name == "goldfish":
            ref = _reference_goldfish_scores(b, sentence)
            delta = max(abs(a - c) for a, c in zip(per_word, ref))
            assert delta < 1e-9, f"goldfish diverges from reference by {delta}"
            print(f"  matches the earlier reference implementation "
                  f"(max |diff| = {delta:.2e})  [OK]")
            notes.append("goldfish reproduces the earlier scorer's "
                         "per-word scores exactly "
                         f"(max abs diff {delta:.1e})")

        # multi-subword probe: proves subwords are summed within a word and that
        # the leading-space offsets still land on the right word
        p_report = b.alignment_report(probe)
        _, p_toks, p_lps, p_offs = b._token_log2_probs(probe)
        print(f"  multi-subword probe ({len(p_toks)} subwords -> "
              f"{len(probe_words)} words):")
        for g in p_report:
            sub = " ".join(f"{t}{o}" for t, o in zip(g["subwords"], g["offsets"]))
            print(f"    {g['index']:>2}  {g['word']:<32} {g['log2p']:>10.4f}   {sub}")
        p_scores = b.score_sentence(probe)
        assert len(p_scores) == len(probe_words), f"{name}: probe length mismatch"
        assert _finite(p_scores), f"{name}: non-finite probe score"
        assert abs(sum(p_scores) - sum(p_lps)) < 1e-6, f"{name}: probe sum mismatch"
        assert len(p_toks) > len(probe_words), \
            f"{name}: probe did not produce multi-subword words"
        assert all(g["subwords"] for g in p_report)
        print(f"    -> {len(p_toks)} subwords grouped into {len(p_scores)} words, "
              f"sums match ({sum(p_scores):.6f})  [OK]")

        scr = b.score_sentence(scrambled)
        assert len(scr) == n and _finite(scr), f"{name}: bad scrambled scores"
        t_orig, t_scr = sum(per_word), sum(scr)
        totals[name] = [t_orig, t_scr]
        counts[name] = len(per_word)

        init = b.initial_word_log2p(words[0])
        assert math.isfinite(init)
        # Mathematically these are identical (causal LM: a prefix's log-probs
        # cannot depend on what follows), but they are computed in float32 at
        # two different sequence lengths, so multi-subword first words differ by
        # ~1e-5 bits of kernel noise. 1e-9 only passes for single-subword words.
        assert abs(init - per_word[0]) < 1e-3, \
            f"{name}: initial_word_log2p disagrees with the first scored word"

        print(f"  initial_word_log2p({words[0]!r}) = {init:.4f} "
              f"(== first word of score_sentence)  [OK]")
        print(f"  total log2 p  original  = {t_orig:10.4f}")
        print(f"  total log2 p  scrambled = {t_scr:10.4f}")
        if t_scr < t_orig:
            print(f"  original > scrambled by {t_orig - t_scr:.4f} bits  [OK]")
        else:
            ok = False
            print("  *** FAIL: scrambled scored >= original ***")
            notes.append(f"{name}: scrambled did not score lower")
        working.append(name)
        print()

    print("=" * 78)
    print("SUMMARY")
    print("=" * 78)
    for name in working:
        o, s = totals[name]
        print(f"  {name:<10} |x|={counts[name]}  orig={o:10.4f}  "
              f"scrambled={s:10.4f}  diff={o - s:9.4f}")

    result = {
        "ok": bool(ok and len(working) == len(BACKEND_SPECS)),
        "script_path": os.path.abspath(__file__),
        "backends_working": working,
        "per_word_counts": counts,
        "original_vs_scrambled": totals,
        "notes": notes,
    }
    print()
    print("JSON " + json.dumps(result))
    return result


if __name__ == "__main__":
    _validate()
