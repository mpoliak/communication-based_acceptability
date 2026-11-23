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
#   - MORCELA = (log p_model - beta * log p_unigram + gamma) / length, with the
#     English beta and gamma that analysis/main_analysis.R learned from the
#     English acceptability ratings, read from
#     analysis/outputs/morcela_weights.csv (row en), so run main_analysis.R first.
#
# Run this from the data/ folder with the goldfish conda environment
# (goldfish_environment.yml, next to this script):
#   cd data && python sample_slor_morcela_scores.py
# It writes data/sample_slor_morcela_scores.csv (one row per example sentence,
# full precision) and prints the same values rounded to 2 decimals.

import math
import os
import warnings

# torch 2.2.2 was built against NumPy 1.x, so importing it next to NumPy 2 warns
# "Failed to initialize NumPy". Nothing here passes tensors to NumPy, so only
# that warning is hidden; every other warning is shown.
warnings.filterwarnings("ignore", message="Failed to initialize NumPy")

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

# ---- MORCELA's English weights, from main_analysis.R ----
weights = pd.read_csv(os.path.join("..", "analysis", "outputs", "morcela_weights.csv"))
beta = weights.loc[weights["experiment"] == "en", "beta"].item()
gamma = weights.loc[weights["experiment"] == "en", "gamma"].item()
print("English MORCELA weights: beta = %.7f, gamma = %.7f" % (beta, gamma))

# ---- Score the examples with the English Goldfish model ----
tokenizer = AutoTokenizer.from_pretrained("goldfish-models/eng_latn_1000mb", use_fast=True)
model = AutoModelForCausalLM.from_pretrained(
    "goldfish-models/eng_latn_1000mb",
    use_safetensors=True,
    low_cpu_mem_usage=True
).eval()
cls_id = tokenizer.convert_tokens_to_ids("[CLS]")

# The same normalization as perceived_sentence (normalize_sentence() in
# get_slor_morcela_scores.py), one step at a time:
# drop punctuation and symbols
normalized = [regex.sub(r"[\p{P}\p{S}]+", "", s) for s in examples["sentence"]]
# remove leading and trailing spaces
normalized = [regex.sub(r"^\s+|\s+$", "", s) for s in normalized]
# collapse runs of whitespace into one space, and lowercase
normalized = [regex.sub(r"\s+", " ", s).lower() for s in normalized]
examples["words"] = [s.split() for s in normalized]
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

examples = examples.drop(columns="words")
examples.to_csv("sample_slor_morcela_scores.csv", index=False)
print(examples[["example", "version", "n_words", "slor", "morcela"]].round(2).to_string(index=False))
print("Wrote data/sample_slor_morcela_scores.csv")
