# Compute the language-model quantities every analysis in the paper uses, for the
# stimuli of all nine experiments and of the artificial English follow-ups: the
# log-probability of each sentence, and the unigram log-probability of the words
# in it.
#
# Both are computed on the PERCEIVED sentence -- the string the participant
# actually saw -- for every condition, not just the original. SLOR and MORCELA
# are properties of the string on the screen, so they use each condition's own
# row. The communication account's Information(s) is the surprisal of the
# intended sentence, so the analyses take it as -log_p_sentence of the item's
# "original" row.
#
# Run this from the data/ folder with the goldfish conda environment
# (goldfish_environment.yml, next to this script):
#   cd data && python get_slor_morcela_scores.py
# It writes data/slor_morcela_scores.csv

import math
import os
import warnings

warnings.filterwarnings("ignore")

import pandas as pd
import regex
import torch
from transformers import AutoTokenizer, AutoModelForCausalLM
from wordfreq import word_frequency

# Paths below are relative to data/, where this script lives
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# One Goldfish monolingual model per language, all trained on 1GB of text, so the
# scores are as comparable across languages as possible. The length-controlled
# English experiment (en_experimental) and the artificial English follow-ups
# (en_artificial) use the English model.
GOLDFISH_MODELS = {
    "da":              "goldfish-models/dan_latn_1000mb",   # Danish (Latin)
    "de":              "goldfish-models/deu_latn_1000mb",   # German (Latin)
    "en":              "goldfish-models/eng_latn_1000mb",   # English (Latin)
    "fr":              "goldfish-models/fra_latn_1000mb",   # French (Latin)
    "hi":              "goldfish-models/hin_deva_1000mb",   # Hindi (Devanagari)
    "ko":              "goldfish-models/kor_hang_1000mb",   # Korean (Hangul)
    "ru":              "goldfish-models/rus_cyrl_1000mb",   # Russian (Cyrillic)
    "zh":              "goldfish-models/zho_hans_1000mb",   # Chinese (Han Simpl.)
    "en_experimental": "goldfish-models/eng_latn_1000mb",   # English (Latin)
    "en_artificial":   "goldfish-models/eng_latn_1000mb",   # English (Latin)
}

# The language code wordfreq uses for each experiment.
WORDFREQ_LANGUAGES = {
    "da": "da", "de": "de", "en": "en", "fr": "fr", "hi": "hi",
    "ko": "ko", "ru": "ru", "zh": "zh", "en_experimental": "en",
    "en_artificial": "en",
}

# The four artificial English experiments (SI 4.1-4.5) all show the same 243
# sentences (27 items x 9 conditions), so we score them once, from the
# acceptability-ratings materials.
EN_ARTIFICIAL_MATERIALS = os.path.join(
    "en_artificial_acceptability_ratings",
    "materials_en_artificial_acceptability_ratings.csv",
)

# Probability given to a word wordfreq has never seen. Danish and Korean need it
# most (about 9% and 10% of their word tokens); English, Mandarin and the
# controlled English experiment need it for none.
UNSEEN_WORD_PROBABILITY = 1e-9

device = "cuda" if torch.cuda.is_available() else "cpu"

_loaded = {}


def get_tokenizer_model(experiment):
    if experiment not in _loaded:
        model_id = GOLDFISH_MODELS[experiment]
        tok = AutoTokenizer.from_pretrained(model_id, use_fast=True)
        # Force safetensors to avoid .bin safety gate; keep memory low on CPU
        mdl = AutoModelForCausalLM.from_pretrained(
            model_id,
            use_safetensors=True,
            low_cpu_mem_usage=True
        ).to(device).eval()
        _loaded[experiment] = (tok, mdl)
    return _loaded[experiment]


def sentence_log_probability(experiment, sentence):
    """Log2 probability of a sentence, summed over all of its tokens.

    The Goldfish model cards advise prepending [CLS], so we do; that token is
    context only and is never itself scored.

    Mandarin is the one special case. Its stimuli are jieba-segmented, because
    that segmentation is what defines n_words and dl_distance. But the Goldfish
    Chinese model is a SentencePiece model trained on running text, so feeding it
    the spaces would make it score every space as a token. We therefore strip the
    spaces before the model sees the sentence. Word segmentation still decides the
    length and the unigram probabilities -- only the model input is unsegmented.
    """
    words = sentence.split()
    text = "".join(words) if experiment == "zh" else " ".join(words)

    tokenizer, model = get_tokenizer_model(experiment)
    cls_id = tokenizer.convert_tokens_to_ids("[CLS]")
    token_ids = tokenizer(text, add_special_tokens=False)["input_ids"]
    input_ids = torch.tensor([[cls_id] + token_ids]).to(device)

    with torch.no_grad():
        logits = model(input_ids=input_ids).logits

    # log p(x_t | x_<t) for each next token (the [CLS] itself is never scored)
    log_probs = torch.log_softmax(logits.float(), dim=-1)
    next_ids = input_ids[:, 1:]
    ctx_log_probs = log_probs[:, :-1, :]
    tok_log_probs = ctx_log_probs.gather(2, next_ids.unsqueeze(-1)).squeeze(-1)

    # Sum in full precision: adding up the tensor first would do it in the
    # model's 32-bit arithmetic and lose about 1e-4 of a bit per sentence.
    return math.fsum(float(lp) for lp in tok_log_probs[0]) / math.log(2.0)


def unigram_log_probability_from_model(experiment, word, cache):
    """Log2 p(word | start of sentence), read off the model itself.

    This is the unigram model we use by default. It is available in every
    language, it is matched to the model's own vocabulary so no word is ever
    unknown, and it needs no outside corpus. It is a sentence-INITIAL
    distribution rather than a true unigram frequency, which is why we also
    compute the wordfreq version below and check that nothing depends on the
    choice.
    """
    if word not in cache:
        cache[word] = sentence_log_probability(experiment, word)
    return cache[word]


def unigram_log_probability_from_wordfreq(experiment, word):
    """Log2 of the word's corpus frequency, from the wordfreq package."""
    frequency = word_frequency(word, WORDFREQ_LANGUAGES[experiment])
    if frequency <= 0.0:
        frequency = UNSEEN_WORD_PROBABILITY
    return math.log2(frequency)


def normalize_sentence(sentence):
    """The same normalization data/preprocess_stimuli.R applies to build
    perceived_sentence: drop punctuation and symbols, trim, collapse repeated
    whitespace, lowercase. The `regex` module is used because Python's `re` has
    no \\p{P} / \\p{S} Unicode classes."""
    sentence = regex.sub(r"[\p{P}\p{S}]+", "", sentence)
    sentence = regex.sub(r"^\s+|\s+$", "", sentence)
    sentence = regex.sub(r"\s+", " ", sentence)
    return sentence.lower()


def load_stimuli(experiment):
    """One row per sentence, with language, item, condition, perceived_sentence."""
    if experiment == "en_artificial":
        # The artificial English materials have no stimuli_data.csv, so we build
        # perceived_sentence here from the critical (exp_b) trials.
        materials = pd.read_csv(EN_ARTIFICIAL_MATERIALS, dtype={"item": str})
        materials = materials[
            (materials["experiment"] == "exp_b") & materials["sentence"].notna()
        ].reset_index(drop=True)
        return pd.DataFrame({
            "language": "en_artificial",
            "item": materials["item"],
            "condition": materials["condition"],
            "perceived_sentence": materials["sentence"].map(normalize_sentence),
        })
    # Each experiment's stimuli live in its own folder, as written by
    # data/preprocess_stimuli.R and data/preprocess_en_experimental.R.
    stimuli = pd.read_csv(os.path.join(experiment, "stimuli_data.csv"))
    return stimuli[stimuli["perceived_sentence"].notna()].reset_index(drop=True)


rows = []

for experiment in GOLDFISH_MODELS:
    stimuli = load_stimuli(experiment)
    print("Scoring " + experiment + ": " + str(len(stimuli)) + " sentences")

    initial_word_cache = {}

    for i, stimulus in stimuli.iterrows():
        sentence = str(stimulus["perceived_sentence"])
        words = sentence.split()

        log_p_sentence = sentence_log_probability(experiment, sentence)
        log_p_unigram_model = sum(
            unigram_log_probability_from_model(experiment, w, initial_word_cache)
            for w in words
        )
        log_p_unigram_wordfreq = sum(
            unigram_log_probability_from_wordfreq(experiment, w) for w in words
        )

        for unigram_name, log_p_unigram in (
            ("model_initial", log_p_unigram_model),
            ("wordfreq", log_p_unigram_wordfreq),
        ):
            rows.append({
                "language": stimulus["language"],
                "item": stimulus["item"],
                "condition": stimulus["condition"],
                "unigram_model": unigram_name,
                "log_p_sentence": log_p_sentence,
                "log_p_unigram": log_p_unigram,
            })

scores = pd.DataFrame(rows)
scores.to_csv("slor_morcela_scores.csv", index=False)
print("Done. Wrote " + str(len(scores)) + " rows to data/slor_morcela_scores.csv")
