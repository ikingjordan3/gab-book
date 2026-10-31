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

# From GWAS to polygenic risk scores

```{admonition} Module at a glance
**Lecture:** Polygenic risk scores 1 &nbsp;·&nbsp; **Prerequisite:** {doc}`../complex-traits-gwas/significance-and-plots`

**After this module you should be able to:** explain how GWAS summary statistics become the input for a polygenic risk score.
```

A genome-wide association study tests millions of SNPs, one at a time, for association with a trait. The
Manhattan plot below is the simulated type 2 diabetes GWAS from {doc}`../complex-traits-gwas/significance-and-plots`.
Each point is a SNP, and each tower of points that rises above the Bonferroni-adjusted threshold of
P = 5 × 10⁻⁸ is a **locus**: a region of the genome associated with type 2 diabetes.

```{figure} ../complex-traits-gwas/images/manhattan.png
:alt: Simulated type 2 diabetes Manhattan plot with peaks above the 5 × 10⁻⁸ line at TCF7L2 (tallest), CDKAL1, FTO, IGF2BP2, SLC30A8, PPARG and JAZF1.
:width: 100%

Illustrative, simulated type 2 diabetes GWAS. Peaks are placed at well-known type 2 diabetes loci; the dashed
line is the Bonferroni-adjusted threshold, P = 5 × 10⁻⁸.
```

Each of these loci changes risk only slightly. Type 2 diabetes is **polygenic**: its heritability is spread
over hundreds of loci, and over many more variants that do not reach genome-wide significance. A polygenic
risk score combines all of these small effects into a single number for each person.

## Summary statistics

A GWAS reports its results as **summary statistics**: a summary of every association tested between a variant
(SNP) and the trait. There is one row per SNP, and each row has two kinds of information.

| Variant (SNP) location information | GWAS analysis results (statistical coefficients) |
|---|---|
| Chromosome | Effect size (β or OR) |
| Position | Standard error |
| rsID | P-value |
| Effect allele, reference allele, alternate allele | |

Summary statistics are the **base data** for a PRS. They contain everything needed to build a score — which
variants to use, which allele counts, and how much weight each allele gets — but no information about any
individual person. That is why large GWAS consortia can share them publicly, and why one GWAS can be reused to
build many scores.

The rest of this chapter follows the path from these summary statistics to a score for each person:
{doc}`prs-defined` defines the score, {doc}`prs-workflow` and {doc}`base-target-data` describe the data, and the
modules of the second lecture describe how the GWAS effect sizes are refined into PRS weights.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/gwas-to-prs.json")
```
