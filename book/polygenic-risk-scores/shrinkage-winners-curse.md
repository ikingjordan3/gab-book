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

# Shrinkage and the winner's curse

```{admonition} Module at a glance
**Lecture:** Polygenic risk scores 2 &nbsp;·&nbsp; **Prerequisite:** {doc}`prs-workflow`

**After this module you should be able to:** explain why GWAS effect sizes are over-estimated and why they must be shrunk before use in a PRS.
```

## What a PRS needs

Recall the PRS from {doc}`prs-defined`:

$$
\widehat{PRS}_i = \sum_{j=1}^{M} X_{ij}\,\hat{\beta}_j
$$

To calculate it we need:

- **Base data:** a list of trait-associated variants ($X_j$) and their effect sizes ($\hat\beta_j$).
- **Target data:** individual genomes (and phenotypes) for which the PRS will be calculated.
- A way to go from GWAS summary statistics to the list of variants and effect sizes used in the score. This
  step — the **PRS calculation** stage of the workflow — is called **shrinkage**.

How do we get there from a GWAS? A GWAS gives an effect size and P-value for millions of SNPs. Using only the
genome-wide significant peaks throws away most of the polygenic signal; using every SNP with its raw estimate
adds a great deal of noise. Shrinkage is the answer to this dilemma.

## Why effect sizes need shrinking

- GWAS are performed on **finite samples** drawn from a larger target population.
- So each SNP effect size estimate includes both the **true effect** and **stochastic variation** (sampling noise).
- The top associations may show the **winner's curse**: their estimates are too large, and they do not generalize
  well to other samples or populations.

This is the problem of **overfitting**: the estimates fit the discovery sample better than they fit anyone else.
It can be addressed by **shrinkage** of the effect size estimates toward zero.

## The winner's curse

The idea comes from auctions. When many bidders estimate the value of the same item, their bids scatter around its
true value. The winner is the bidder who overestimated the most, so the winning bid is more than the item is worth.
The actual value of the item is best approximated by the **mean** of all the bids, not by the winning bid.

Now think of each bid as an **independent GWAS** of the same SNP, each on a different sample from the larger
population. The estimates scatter around the true effect, and the actual effect size is best estimated as the mean
of the estimates across many independent GWAS. But a SNP only enters the base data if its estimate passed the
significance threshold — that is, if it landed in the upper tail. The base-data effect size is the winning bid.

```{figure} images/winners_curse.png
:alt: Sampling distribution of GWAS effect size estimates for one SNP, centred on the actual effect size; the base-data GWAS estimate lies in the shaded upper tail that reaches p < 5 × 10⁻⁸.
:width: 560px

The winner's curse for one SNP. Illustrative: exact sampling distribution of β̂ with true β = 0.15 and SE = 0.035.
The shaded tail contains the estimates that reach P < 5 × 10⁻⁸.
```

You can see the size of the bias with a quick simulation of many GWAS of the same SNP:

```{code-cell} ipython3
import numpy as np
from scipy.stats import norm

rng = np.random.default_rng(1)
true_beta, se = 0.15, 0.035
estimates = rng.normal(true_beta, se, size=100_000)      # 100,000 independent GWAS of one SNP

z_threshold = norm.isf(5e-8 / 2)                          # |z| needed for P < 5e-8
significant = estimates / se > z_threshold

print(f"true effect:                         {true_beta:.3f}")
print(f"mean of all estimates:               {estimates.mean():.3f}")
print(f"GWAS reaching P < 5e-8:              {significant.mean():.1%}")
print(f"mean estimate among those 'winners': {estimates[significant].mean():.3f}")
```

Averaged over all GWAS, the estimate is unbiased. Among the GWAS where the SNP was significant, the average
estimate is about 40% larger than the true effect (0.21 instead of 0.15). The weaker the true effect relative to its standard error, the
worse the inflation — which is exactly the situation for most SNPs in a polygenic score.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/shrinkage-winners-curse.json")
```
