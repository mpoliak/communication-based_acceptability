import math
import torch
from transformers import AutoTokenizer, AutoModelForCausalLM
import pandas as pd


# One Goldfish monolingual model per language, all trained on 1GB of text so that
# surprisal is as comparable across languages as possible. The length-controlled
# English experiment (en_experimental) uses the same model as English.
#   Model cards (details + advice to prepend [CLS]):
#   https://huggingface.co/goldfish-models/eng_latn_1000mb
GOLDFISH_MODELS = {
    "en": "goldfish-models/eng_latn_1000mb",   # English (Latin)
    "da": "goldfish-models/dan_latn_1000mb",   # Danish (Latin)
    "de": "goldfish-models/deu_latn_1000mb",   # German (Latin)
    "fr": "goldfish-models/fra_latn_1000mb",   # French (Latin)
    "hi": "goldfish-models/hin_deva_1000mb",   # Hindi (Devanagari)
    "ko": "goldfish-models/kor_hang_1000mb",   # Korean (Hangul)
    "ru": "goldfish-models/rus_cyrl_1000mb",   # Russian (Cyrillic)
    "zh": "goldfish-models/zho_hans_1000mb",   # Chinese (Han Simpl.)
    "en_experimental": "goldfish-models/eng_latn_1000mb",   # English (Latin)
}

device = "cuda" if torch.cuda.is_available() else "cpu"

## Read materials
materials = pd.read_csv("materials_for_surprisal.csv")  # columns: language, item, sentence

## Add [CLS] to all sentences (preserve existing cols; create if missing)
if "sentence_prefixed" not in materials.columns:
    materials["sentence_prefixed"] = "[CLS]" + materials["sentence"]
else:
    # normalize to ensure a single space after [CLS]
    materials["sentence_prefixed"] = (
        materials["sentence_prefixed"]
        .astype(str)
        .apply(lambda s: "[CLS] " + s.lstrip().removeprefix("[CLS]").lstrip())
    )

## Create/overwrite a column to store sum surprisal per sentence
materials["sum_surprisal"] = 0.0


# --------------------------- minimal multi-sentence scorer ---------------------------

# cache per-language tokenizer/model
_loaded = {}

def get_tokenizer_model(lang_code: str):
    if lang_code not in _loaded:
        model_id = GOLDFISH_MODELS[lang_code]
        tok = AutoTokenizer.from_pretrained(model_id, use_fast=True)
        # Force safetensors to avoid .bin safety gate; keep memory low on CPU
        mdl = AutoModelForCausalLM.from_pretrained(
            model_id,
            use_safetensors=True,
            low_cpu_mem_usage=True
        ).to(device).eval()
        _loaded[lang_code] = (tok, mdl)
    return _loaded[lang_code]


for idx, r in materials.iterrows():
    lang = str(r["language"]).strip()
    s = str(r["sentence_prefixed"])

    # minimal normalization (idempotent)
    sp = s.lstrip()
    if not sp.startswith("[CLS]"):
        sp = "[CLS] " + sp
    elif not sp.startswith("[CLS] "):
        sp = sp.replace("[CLS]", "[CLS] ", 1)

    tokenizer, model = get_tokenizer_model(lang)

    # tokenize and score
    enc = tokenizer(sp, return_tensors="pt", add_special_tokens=False)
    input_ids = enc["input_ids"].to(device)

    with torch.no_grad():
        logits = model(input_ids=input_ids).logits  # [1, seq, vocab]

    # log p(x_t | x_<t) for each next token (skip the very first token)
    log_probs = torch.log_softmax(logits, dim=-1)
    next_ids = input_ids[:, 1:]
    ctx_log_probs = log_probs[:, :-1, :]
    tok_log_probs = ctx_log_probs.gather(2, next_ids.unsqueeze(-1)).squeeze(-1)  # [1, T-1]

    # total sentence surprisal in bits
    total_bits = float(sum((-lp) / math.log(2.0) for lp in tok_log_probs[0]))
    materials.at[idx, "sum_surprisal"] = total_bits

# Save the updated dataframe with the new column
materials = materials.drop(columns=["sentence_prefixed"])
materials.to_csv("surprisal_goldfish.csv", index=False)
print("Done. Wrote sentence-level sums to surprisal_goldfish.csv")
