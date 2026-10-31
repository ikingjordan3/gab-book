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

# Hard and soft thresholding

```{admonition} Module at a glance
**Lecture:** Polygenic risk scores 2 &nbsp;·&nbsp; **Prerequisite:** {doc}`shrinkage-winners-curse`

**After this module you should be able to:** contrast P-value thresholding (hard) with shrinking all effects (soft).
```

Three facts set the problem:

- Variant effect sizes are estimated with uncertainty (the winner's curse).
- Not all variants influence the trait of interest.
- Using unadjusted effect size estimates in a PRS could give a poorly estimated score with high variance.

There are two broad shrinkage strategies:

1. **Hard thresholding:** use a P-value threshold as the inclusion criterion for SNPs in the score.
2. **Soft thresholding:** shrink the effect estimates of **all** variants, with standard or tailored statistical
   techniques.

## Effect sizes by rank

The plots in this module show a simulated GWAS of 2,000 variants, 5% of them causal. Each variant's estimated
effect (as a z-score) is plotted against its rank, with positive and negative estimates ranked separately from the
largest down. A few large estimates sit at the left; most are small and near zero.

```{figure} images/shrink_original.png
:alt: GWAS effect size estimates (z) of 2,000 variants by rank on a log scale; a few large positive and negative estimates at rank 1, falling toward zero by rank 1000.
:width: 420px

Variant effects estimated with uncertainty. Simulated: z-scores are the exact quantiles of a mixture of 95% null and
5% causal variants. Rank is on a log scale so that the strongest hits are visible.
```

## Hard thresholding

**P-value thresholding** keeps only variants with P below the threshold, at their full GWAS effect size. The
effect sizes of all other variants shrink to 0. This is "hard" thresholding: each variant is either in, unchanged,
or out.

```{figure} images/shrink_hard_strict.png
:alt: Original (grey) and adjusted (orange) effect sizes by rank at P < 5 × 10⁻⁸; only 8 of 2,000 variants keep their effect and all others are set to zero.
:width: 420px

Hard thresholding at P < 5 × 10⁻⁸: 8 of 2,000 variants are kept.
```

As the P-value threshold increases, more and more variants are included:

```{figure} images/shrink_hard_lenient.png
:alt: Original and adjusted effect sizes by rank at P < 0.05; 148 of 2,000 variants keep their effect.
:width: 420px

Hard thresholding at P < 0.05: 148 of 2,000 variants are kept.
```

The best threshold is not known in advance, so scores are usually built at several thresholds
(for example 5 × 10⁻⁸, 10⁻⁵, 10⁻³, 0.05, 0.1, 0.5 and 1) and the best is chosen in validation data
({doc}`prs-workflow`).

```{code-cell} ipython3
import numpy as np
from scipy.stats import norm

rng = np.random.default_rng(7)
n = 2000
causal = rng.uniform(size=n) < 0.05
z = rng.normal(0, 1, n) + np.where(causal, rng.normal(0, 3, n), 0)   # 95% null, 5% causal
p = 2 * norm.sf(np.abs(z))

for threshold in [5e-8, 1e-5, 1e-3, 0.05, 0.5]:
    kept = p < threshold
    print(f"P < {threshold:<6g}: {kept.sum():5d} SNPs kept, {causal[kept].mean():5.0%} of them causal")
```

Raising the threshold brings in more of the causal variants, but an ever larger share of null ones — the
trade-off that the validation data resolve.

## Soft thresholding

More sophisticated techniques shrink the effects of **all** variants. This is "soft" thresholding: no variant is
simply in or out; each estimate is pulled toward zero by an amount that depends on how well supported it is.
Examples are LASSO and Bayesian methods such as LDpred ({doc}`shrinkage-methods`). Soft thresholding includes many
more variants in the PRS — millions, in the scores of Khera et al. ({doc}`prs-workflow`).

```{figure} images/shrink_soft.png
:alt: Three panels of Bayesian shrinkage of the same effect sizes for assumed fractions of causal variants p = 1, 0.05 and 0.005. With p = 1 every effect is shrunk by the same proportion; with smaller p, small effects are shrunk almost to zero while the largest are kept.
:width: 100%

Soft thresholding of the same 2,000 variants: the posterior mean effect under a Bayesian (point-normal) prior, as
in LDpred, for three assumed fractions of causal variants p.
```

The panels show the two behaviours of Bayesian shrinkage:

- With **p = 1** (the infinitesimal model: every variant has some effect), every estimate is shrunk by the same
  proportion, like ridge regression.
- With **small p** (few variants with larger effects), small, noisy estimates are shrunk almost to zero while large,
  well-supported ones are left nearly untouched — a smooth version of hard thresholding.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/hard-soft-thresholding.json")
```
