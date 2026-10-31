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

# Shrinkage methods

```{admonition} Module at a glance
**Lecture:** Polygenic risk scores 2 &nbsp;·&nbsp; **Prerequisite:** {doc}`hard-soft-thresholding`

**After this module you should be able to:** name the main shrinkage methods (LASSO, ridge, LDpred) and explain why GWAS power matters more than the method.
```

## PRS methods that shrink all SNPs

Methods for soft thresholding use common statistical shrinkage techniques. The details are beyond the scope of
this course, but the main ideas are easy to state:

| Method | How it shrinks effect sizes |
|---|---|
| **LASSO regression** | Reduces small effect sizes to exactly zero (and shrinks the rest). |
| **Ridge regression** | Shrinks every effect size, but does not reduce any of them to zero. |
| **Bayesian methods** | Shrink effect sizes by specifying a **prior distribution** for them and using the posterior estimate. LDpred is an example. |

## Effect size shrinkage with LDpred

**LDpred** (Vilhjálmsson et al., [2015](https://doi.org/10.1016/j.ajhg.2015.09.001)) is a Bayesian approach to
effect size shrinkage:

- It starts with the effect sizes from the original GWAS and assumes a prior distribution (Gaussian) for the true
  effect sizes.
- A parameter $p$ sets the **proportion of causal variants**, that is, what proportion of variants should have
  non-zero effect sizes.
  - **Infinitesimal model:** many variants with small effects (high $p$, up to 1).
  - **Non-infinitesimal model:** few variants with larger effects (low $p$).
- An **LD matrix** gives the correlations between variants, estimated from a reference panel, so that correlated
  SNPs do not count the same signal twice.
- From $p$ and the LD matrix it estimates the **posterior mean** effect size of every variant.

The panels in {doc}`hard-soft-thresholding` show what $p$ does: $p = 1$ shrinks every effect by the same
proportion, while small $p$ keeps the largest effects and shrinks the rest almost to zero.

### LDpred in practice: Khera et al. (2018)

For four of the five diseases studied by Khera et al. ({doc}`prs-workflow`), the best score was built with LDpred
and included 6.6–6.9 million variants. Their tuning parameter **ρ** is LDpred's $p$:

| Disease | Variants in score | LDpred ρ (proportion of causal variants) |
|---|---|---|
| Coronary artery disease | 6,630,150 | 0.001 |
| Atrial fibrillation | 6,730,541 | 0.003 |
| **Type 2 diabetes** | **6,917,436** | **0.01** |
| Inflammatory bowel disease | 6,907,112 | 0.1 |

The type 2 diabetes score (PGS000014 in the PGS Catalog, {doc}`base-target-data`) assumed that 1% of variants are
causal.

## Considerations for shrinkage methods

- The optimal shrinkage parameters are unknown, just as the best P-value threshold is unknown in clumping and
  thresholding. So PRS should be calculated across a **range of parameter values** and the best chosen in
  validation data.
- Shrinkage methods include many more variants, so they **should be superior**…
- …but they **may not perform much better** than the classic PRS method, clumping and thresholding
  ({doc}`clumping-thresholding`).

## Data matter more than methods

The **power of the original GWAS** used for the effect size estimates is much more important than the PRS
calculation method:

- **Garbage in, garbage out** — always true in bioinformatics.
- **Data >> methods.** A larger, better-powered and more diverse base GWAS, such as the type 2 diabetes
  meta-analysis of Suzuki et al. (2024) in {doc}`base-target-data`, improves a score far more than a cleverer method.

The practical recommendation is to try **multiple GWAS variant and effect size sets**, along with different methods
and parameters, and compare them in validation data. The PGS Catalog is a good source of already-validated variant
and weight sets: its type 2 diabetes page alone lists about 250 scores built with LDpred, PRS-CS
([Ge et al., 2019](https://doi.org/10.1038/s41467-019-09718-5)), clumping and thresholding and other methods.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/shrinkage-methods.json")
```
