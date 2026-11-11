# Compute SLOR and MORCELA for the example sentences in the main text: the "ball"
# and "demonstrated" examples, each in its error-free, 1-exchange, 3-exchange and
# shuffled version. These are illustrations written for the paper, not all of
# them stimuli, so they are scored here rather than read from
# slor_morcela_scores.csv.
#
# The steps are those of get_slor_morcela_scores.py and analysis/main_analysis.R:
#   - each sentence is normalized as perceived_sentence is (no punctuation,
#     lowercased) and scored by the English Goldfish model with [CLS] prepended,
#     in log2 probabilities;
#   - the unigram term is the model-initial one: the sum over the sentence's
#     words of log2 p(word | start of sentence);
#   - SLOR = (log p_model - log p_unigram) / length, in words;
#   - MORCELA = (log p_model - beta * log p_unigram + gamma) / length, with beta
#     and gamma learned from the English acceptability ratings exactly as in
#     main_analysis.R: regress the item-by-condition mean ratings (DL distance
#     up to 3, after the comprehension exclusions) on log p_model / length,
#     log p_unigram / length and 1 / length.
#
# Run this from the data/ folder with the goldfish conda environment
# (goldfish_environment.yml, next to this script):
#   cd data && python sample_slor_morcela_scores.py
# It only prints; it writes no files.

import math
import os
import warnings

warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import regex
import torch
from transformers import AutoTokenizer, AutoModelForCausalLM

# Paths below are relative to data/, where this script lives
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# The examples as they appear in the main text
examples = pd.DataFrame([
    ("ball", "a", "A ball flying in the air can hurt."),
    ("ball", "b", "A ball flying in air the can hurt."),
    ("ball", "c", "A ball flying air in can the hurt."),
    ("ball", "d", "In flying can ball a hurt air the."),
    ("demonstrated", "a", "You will already have demonstrated a pretty long successful record of writing stuff for us before you will be in the room on a daily basis."),
    ("demonstrated", "b", "You will already have demonstrated a pretty long successful record of writing stuff for us before you will be in room the on a daily basis."),
    ("demonstrated", "c", "You will already have a demonstrated pretty long successful record of writing stuff for before us you will be in room the on a daily basis."),
    ("demonstrated", "d", "Stuff basis of will daily you for in have writing room a us be long will pretty before on the already record you demonstrated a successful."),
], columns=["example", "version", "sentence"])

# ---- MORCELA's English weights, as in main_analysis.R ----
dat = pd.read_csv("full_data.csv", low_memory=False)

# Exclude participants with a low comprehension score
comprehension = dat[dat["experiment"].str.contains("comprehension")]
correct = comprehension["literal_response"] == comprehension["response"]
accuracy = correct.groupby(comprehension["participant_id"]).mean()
exclusions = accuracy[accuracy < .8].index

english = dat[
    ~dat["participant_id"].isin(exclusions)
    & (dat["experiment"] == "exp_b")
    & ~dat["condition"].str.contains("non_sense")
    & (dat["language"] == "en")
].copy()
english["response"] = english["response"].str[0].astype(float)

english_by_items = (
    english[english["dl_distance"] <= 3]
    .groupby(["item", "condition", "n_words"], as_index=False)["response"]
    .mean()
)

scores = pd.read_csv("slor_morcela_scores.csv", dtype={"item": str})
scores = scores[(scores["language"] == "en") & (scores["unigram_model"] == "model_initial")].copy()
scores["item"] = "en_" + scores["item"]
english_by_items = english_by_items.merge(
    scores[["item", "condition", "log_p_sentence", "log_p_unigram"]],
    on=["item", "condition"], how="left", validate="one_to_one",
)
assert not english_by_items["log_p_sentence"].isna().any()

predictors = np.column_stack([
    np.ones(len(english_by_items)),
    english_by_items["log_p_sentence"] / english_by_items["n_words"],
    english_by_items["log_p_unigram"] / english_by_items["n_words"],
    1 / english_by_items["n_words"],
])
coefficients = np.linalg.lstsq(predictors, english_by_items["response"].to_numpy(), rcond=None)[0]
beta = -coefficients[2] / coefficients[1]
gamma = coefficients[3] / coefficients[1]
print("English MORCELA weights: beta = %.7f, gamma = %.7f" % (beta, gamma))

# ---- Score the examples with the English Goldfish model ----
tokenizer = AutoTokenizer.from_pretrained("goldfish-models/eng_latn_1000mb", use_fast=True)
model = AutoModelForCausalLM.from_pretrained(
    "goldfish-models/eng_latn_1000mb",
    use_safetensors=True,
    low_cpu_mem_usage=True
).eval()
cls_id = tokenizer.convert_tokens_to_ids("[CLS]")

# The same normalization as perceived_sentence
examples["words"] = [
    regex.sub(r"\s+", " ", regex.sub(r"^\s+|\s+$", "", regex.sub(r"[\p{P}\p{S}]+", "", s))).lower().split()
    for s in examples["sentence"]
]
examples["n_words"] = examples["words"].map(len)

# Every string we need a log probability for: each whole sentence, and each word
# on its own (the model-initial unigram term)
strings = list(dict.fromkeys(
    [" ".join(w) for w in examples["words"]] + [x for w in examples["words"] for x in w]
))
log_p = {}
for text in strings:
    token_ids = tokenizer(text, add_special_tokens=False)["input_ids"]
    input_ids = torch.tensor([[cls_id] + token_ids])
    with torch.no_grad():
        logits = model(input_ids=input_ids).logits
    log_probs = torch.log_softmax(logits.float(), dim=-1)
    tok_log_probs = log_probs[:, :-1, :].gather(2, input_ids[:, 1:].unsqueeze(-1)).squeeze(-1)
    log_p[text] = math.fsum(float(lp) for lp in tok_log_probs[0]) / math.log(2.0)

examples["log_p_sentence"] = [log_p[" ".join(w)] for w in examples["words"]]
examples["log_p_unigram"] = [math.fsum(log_p[x] for x in w) for w in examples["words"]]
examples["slor"] = (examples["log_p_sentence"] - examples["log_p_unigram"]) / examples["n_words"]
examples["morcela"] = (
    examples["log_p_sentence"] - beta * examples["log_p_unigram"] + gamma
) / examples["n_words"]

print(examples[["example", "version", "n_words", "slor", "morcela"]].round(2).to_string(index=False))
