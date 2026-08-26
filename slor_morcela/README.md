# SLOR and MORCELA

Two language-model accounts of acceptability, fit to the same data and compared
against the three models in the paper.

* **SLOR** (Lau, Clark & Lappin 2017) is a sentence's log probability under a
  language model, corrected for how frequent its words are and divided by how
  long it is.
* **MORCELA** (Tjuatja, Neubig, Linzen & Hao 2025) is the same thing with the two
  corrections allowed to be weaker or stronger than SLOR assumes, by learning
  their weights from the ratings.

Both are asked here to do what the communication-based model does: predict how
acceptable a sentence is.

## Running it

1. `get_slor_morcela_scores.py` scores every stimulus with that language's
   Goldfish model and writes `../data/slor_morcela_scores.csv`. Use the Conda
   environment in `supporting_information/goldfish_surprisals/goldfish_environment.yml`,
   plus the `wordfreq` package. This step takes a few minutes and needs to be run
   only once.
2. `slor_morcela_analysis.R` fits the models and draws the figures. Run it in a
   UTF-8 locale, as with the other scripts here.

Fitted models are cached in `models/`; delete that folder to refit.

## What it produces

| File | Contents |
|---|---|
| `figures/main_figure_slor_morcela.png` / `.pdf` | The main figure with SLOR and MORCELA added next to the human ratings |
| `figures/likelihood_differences_slor_morcela.png` | Model comparison for each experiment |
| `figures/likelihoods_by_language.csv` | The numbers behind that figure |

## Two things worth knowing about the figure

**The scores are computed on the sentence the participant saw**, corrupted words
and all. The surprisals in `supporting_information/goldfish_surprisals/` cannot
be reused, because they score only the uncorrupted sentence, which would give
every condition of an item the same value.

**The SLOR and MORCELA panels are drawn differently from the other three.** The
three models from the paper are functions of length and number of errors, so we
can ask them what they predict for any combination of the two, and that is what
their panels show. SLOR and MORCELA are functions of a sentence's probability:
many sentences share a length and an error count but have quite different
scores, so there is no single sentence to ask about. Their panels instead show
what those models predict for the sentences that actually exist, averaged at
each length and then smoothed the same way the human ratings are.
