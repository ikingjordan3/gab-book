---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Risk stratification, not diagnosis

```{admonition} Module at a glance
**Lecture:** Polygenic risk scores 1 &nbsp;·&nbsp; **Prerequisite:** {doc}`prs-defined`

**After this module you should be able to:** explain why a PRS stratifies risk rather than diagnoses disease, and distinguish relative risk, absolute risk and AUC.
```

## The PRS distribution

A PRS adds up many small effects, so in a population the scores are approximately **normally distributed** — the
polygenic model of {doc}`../complex-traits-gwas/polygenic-model`. Scores are usually standardized to a mean of 0
and a standard deviation (SD) of 1. Most people sit near the middle; a small fraction in each tail have markedly
lower or higher genetic risk.

```{figure} images/prs_distribution.png
:alt: Normal distribution of polygenic risk scores from low risk (left) to high risk (right); the top 1%, more than 2.33 SD above the mean, is shaded orange.
:width: 100%

PRS standardized to mean 0 and SD 1. The top 1% of the population (orange) lies more than 2.33 SD above the mean.
```

## Risk rises steeply in the top percentiles

Khera et al. ([2018](https://doi.org/10.1038/s41588-018-0183-z)) built genome-wide polygenic scores for five
common diseases in UK Biobank. For type 2 diabetes, the prevalence of disease rose from about 0.3% in the bottom
percentile of the score to about 6% in the top percentile, against an average of 2.4%. They found that 3.5% of
the population had at least three times the average risk of type 2 diabetes.

The figure below reproduces this pattern with the liability-threshold model: disease occurs when liability (the
PRS plus everything else) exceeds a threshold, with the overall prevalence and the prevalence in the top 1% set
to Khera's values.

```{figure} images/prs_risk_gradient.png
:alt: Type 2 diabetes prevalence by PRS percentile, rising from about 0.7% in the bottom percentile to 6% in the top 1%, against a population average of 2.4%.
:width: 520px

Type 2 diabetes prevalence by PRS percentile. Exact calculation under the liability-threshold model, calibrated to
Khera et al. (2018), Fig. 3b: prevalence 2.4%, top 1% ≈ 6%. To fit these values the PRS must explain about 3% of
the variance in liability.
```

You can reproduce the curve yourself:

```{code-cell} ipython3
import numpy as np
from scipy.stats import norm
from scipy.integrate import quad

K, r2 = 0.024, 0.027                 # T2D prevalence; share of liability variance explained by the PRS
T = norm.isf(K)                      # liability threshold

def prevalence(lo, hi):
    """Prevalence among people whose PRS percentile lies between lo and hi."""
    f = lambda u: norm.sf((T - np.sqrt(r2) * norm.ppf(u)) / np.sqrt(1 - r2))
    return quad(f, lo, hi)[0] / (hi - lo)

top1 = prevalence(0.99, 1.0)
rest = prevalence(0.0, 0.99)
print(f"top 1%: {top1:.1%}   everyone else: {rest:.1%}   ratio: {top1 / rest:.1f}x")
```

## Why a PRS is not a diagnostic test

People in the top 1% of the type 2 diabetes PRS have elevated prevalence, but only about **6%** of them have
the disease. If the PRS were used as a diagnostic test, about **94% of the people it flagged would be false
positives**.

Two ways of describing risk explain the apparent paradox:

- **Relative risk** compares groups. People in the top 1% have about 2.5 times the risk of everyone else — a
  large relative difference.
- **Absolute risk** is the chance that a given person develops the disease. For the top 1% it is still only
  about 6%, because type 2 diabetes in this cohort is uncommon to begin with (2.4%).

A PRS shifts the probability of disease; it does not determine it. Its proper use is **risk stratification**,
a preventive measure:

- **Early intervention**, such as lifestyle programmes for people at high genetic risk.
- **Enhanced screening**, such as earlier or more frequent HbA1c testing.

## How well does a score separate cases from controls? AUC

The usual summary of a score's predictive performance is the **AUC** (area under the receiver operating
characteristic curve). It is the probability that a randomly chosen case has a higher score than a randomly
chosen control:

- AUC = 0.5: no better than chance.
- AUC = 1: perfect separation.

Khera's type 2 diabetes score had an AUC of 0.72 in UK Biobank, and that figure comes from a model that also
included age, sex, genotyping array and four ancestry principal components. Like the risk gradient, the AUC shows
a real but modest separation: useful for stratifying a population, not for diagnosing an individual.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/risk-stratification.json")
```
