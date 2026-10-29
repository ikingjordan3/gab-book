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

# The case-control design

```{admonition} Module at a glance
**Lecture:** Genome-wide association studies &nbsp;·&nbsp; **Prerequisite:** {doc}`heritability-to-gwas`

**After this module you should be able to:** describe the case-control GWAS design, count alleles in cases and controls, and identify the effect allele at a SNP.
```

Most GWAS of disease use the simplest possible comparison: people who have the disease against people who do
not. This module describes that design and the first step of the analysis, counting alleles.

## The design

```{figure} images/case_control_design.png
:alt: Case-control GWAS design. DNA from cases (patients) and controls (non-patients) is genotyped on arrays that measure hundreds of thousands of SNPs per sample, giving genotype calls CC, CT or TT at each SNP. Allele frequencies are then compared between cases and controls to find disease-associated SNPs.
:width: 100%

The case-control GWAS design.
```

1. **Recruit cases and controls.** Cases have the disease, for example type 2 diabetes diagnosed by a
   physician or by blood glucose. Controls do not. Both groups should come from the same population, so that
   the only systematic difference between them is the disease; {doc}`covariates` explains what happens when
   they are not.
2. **Genotype everyone.** DNA from each person is genotyped on a **SNP array** that measures hundreds of
   thousands of SNPs per sample (more are then filled in statistically, by imputation). At each SNP, each
   person gets a genotype call such as CC, CT or TT.
3. **Compare allele frequencies.** At every SNP, test whether the allele frequencies differ between cases and
   controls.

The underlying question is the same at every SNP:

- **Are there genetic variants (alleles or genotypes) that are found at significantly higher frequency in
  cases than in controls?**
- If yes, such variants are said to be **associated with** the trait (disease) of interest.
- The **effect allele** is the allele that is associated with the trait. When it is more common in cases, it
  is also called the **risk allele**.

"Associated with" is deliberately careful wording. An associated SNP may cause the disease, or it may simply
be in linkage disequilibrium with a nearby causal variant, which is why GWAS identify **loci** (regions)
rather than single causal variants ({doc}`significance-and-plots`).

## Counting alleles

For each genotyped SNP, the analysis counts the number of each allele in cases and in controls. Each person
carries two alleles, so 1,000 people contribute 2,000 alleles.

```{figure} images/allele_frequencies.png
:alt: 1,000 cases and 1,000 controls at one SNP. In cases, C is 62% of alleles and T is 38%. In controls, C is 49% and T is 51%.
:width: 520px

Allele frequencies at one SNP in 1,000 cases and 1,000 controls.
```

For the SNP in the figure, start from the genotype counts and convert them to allele counts:

| | CC | CT | TT | C alleles | T alleles | freq(C) |
|---|---|---|---|---|---|---|
| **Cases** (n = 1,000) | 384 | 472 | 144 | 2 × 384 + 472 = **1,240** | 472 + 2 × 144 = **760** | 0.62 |
| **Controls** (n = 1,000) | 240 | 500 | 260 | 2 × 240 + 500 = **980** | 500 + 2 × 260 = **1,020** | 0.49 |

C makes up 62% of alleles in cases but only 49% in controls. **For this SNP, C is the effect allele.** Whether
a difference of this size is large enough to matter, and unlikely enough to be chance, is the subject of the
next module, {doc}`odds-ratios`.

The same counting, in code:

```{code-cell} ipython3
genotypes = {                 # counts of CC, CT, TT
    "cases":    {"CC": 384, "CT": 472, "TT": 144},
    "controls": {"CC": 240, "CT": 500, "TT": 260},
}

for group, g in genotypes.items():
    n_C = 2 * g["CC"] + g["CT"]          # each CC person carries two C alleles
    n_T = 2 * g["TT"] + g["CT"]
    print(f"{group:8s}: C = {n_C:,}, T = {n_T:,}, freq(C) = {n_C / (n_C + n_T):.2f}")
```

A GWAS repeats exactly this calculation, followed by a statistical test, at every one of hundreds of
thousands to millions of SNPs.

## What makes a good case-control study

- **Large samples.** Effects of common variants are small, so tens to hundreds of thousands of people are
  needed to detect them at genome-wide significance.
- **Careful phenotyping.** Misclassifying cases as controls (or the reverse) dilutes the difference between
  the groups.
- **Matched populations.** Cases and controls should have the same ancestry and come from the same source
  population, or allele-frequency differences unrelated to disease will look like associations.
- **Quality control.** Poorly genotyped SNPs or samples can create false differences, especially if cases and
  controls were genotyped in different batches.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/case-control-design.json")
```
