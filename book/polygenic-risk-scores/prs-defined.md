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

# Polygenic risk scores defined

```{admonition} Module at a glance
**Lecture:** Polygenic risk scores 1 &nbsp;·&nbsp; **Prerequisite:** {doc}`gwas-to-prs`

**After this module you should be able to:** define a PRS, write it as a weighted sum of effect alleles, and distinguish the effect allele from the risk allele.
```

A **polygenic risk score (PRS)** is a single value estimate of an individual's genetic "risk" for a trait.

- "Risk" here refers to the risk of disease, such as type 2 diabetes.
- The same idea can be used for the genetic predisposition for **any** trait, including quantitative traits
  such as height or BMI. The more general names are **polygenic score (PGS)** and
  **genome-wide polygenic score (GPS)**.

## A weighted sum of effect alleles

A PRS is calculated as the **weighted sum** of trait-associated effect alleles in the genome:

$$
\widehat{PRS}_i = \sum_{j=1}^{M} X_{ij}\,\hat{\beta}_j
$$

| Symbol | Meaning |
|---|---|
| $\widehat{PRS}_i$ | PRS for individual $i$ |
| $\sum_{j=1}^{M}$ | sum across the $M$ SNPs in the score |
| $X_{ij}$ | number of effect alleles (0, 1 or 2) that individual $i$ carries at SNP $j$ |
| $\hat{\beta}_j$ | GWAS-estimated effect size of SNP $j$ (its weight) |

For each person, count the effect alleles at each SNP, multiply by that SNP's effect size, and add up over all
$M$ SNPs. The hats (^) mark estimates: $\hat\beta_j$ is estimated in a GWAS, so the PRS is itself an estimate.
This is the same formula as $\text{PRS} = \sum_j \beta_j G_j$ in {doc}`../complex-traits-gwas/heritability-to-gwas`,
written with the notation used in the lectures.

For a binary trait such as type 2 diabetes the weight is usually the **log odds ratio**, $\hat\beta = \ln(\text{OR})$.
Logs turn the multiplicative effects of odds ratios into effects that can be added up.

## Variants, effect alleles and weights

The variants, effect alleles and weights all come from **GWAS summary statistics**. The weights can be betas (β)
or odds ratios (OR). Two terms describe which allele is counted, and they are easy to confuse:

- **Effect allele:** the allele that is coded for the association. If the association model regresses the trait
  on the number of alternate alleles {0, 1, 2}, then the alternate allele is the effect allele. Its effect can
  increase or decrease risk.
- **Risk allele:** the allele that is associated with an increase in the trait — an increase in the risk of
  disease, or in the value of a quantitative trait. Its OR is always greater than 1.

At a SNP where the alternate allele is protective, the effect allele and the risk allele are different alleles,
and the OR of one is the inverse of the OR of the other. {doc}`base-target-data` works through real type 2 diabetes
examples.

Building a score then takes two more steps:

1. Variants are selected and the effect-allele weights are refined (the subject of the second lecture,
   starting with {doc}`shrinkage-winners-curse`).
2. The refined variants are applied to individual-level genotype (and phenotype) data — the **target data**.

## A first example

Here is a three-SNP score for one person, using the log odds ratios as weights:

```{code-cell} ipython3
import math

odds_ratios = [1.37, 0.89, 0.88]   # OR for the counted (effect) allele at each SNP
genotypes   = [2, 0, 1]            # copies of the effect allele carried by this person

weights = [math.log(o) for o in odds_ratios]
prs = sum(x * b for x, b in zip(genotypes, weights))
print("weights:", [round(b, 3) for b in weights])
print(f"PRS = {prs:.3f}")
```

A real genome-wide score does exactly the same thing at thousands to millions of SNPs.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/prs-defined.json")
```
