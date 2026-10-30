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

# From heritability to GWAS

```{admonition} Module at a glance
**Lecture:** Both lectures (end of *Genetics of complex, common disease*; start of *Genome-wide association studies*)
&nbsp;·&nbsp; **Prerequisite:** {doc}`heritability`

**After this module you should be able to:** describe what a GWAS looks for, connect heritability to the specific loci GWAS find, describe what GWAS reveal about genetic architecture, and explain how GWAS results feed SNP-heritability and polygenic risk scores.
```

The modules so far established that complex traits are heritable and that their genetic basis is many variants
of small effect. This module bridges the chapter's two lectures. It introduces the tool that finds those
variants, the **genome-wide association study (GWAS)**, and previews what GWAS results are used for. The rest of
the chapter works through how a GWAS is done, from study design to reading the results.

## How do we study complex traits that are heritable?

**Genome-wide association studies (GWAS)** seek to identify the single nucleotide polymorphisms (SNPs) that are
common in the human genome and contribute to the variance in observed phenotypes. Very simply, this means finding
**SNPs that vary in frequency between cases and controls**.

At each of hundreds of thousands to millions of SNPs, a GWAS asks whether one allele is more common in people
with the disease (cases) than in people without it (controls). For each SNP it reports:

- an **effect size** ($\beta$, or an odds ratio for a disease), which says how much each copy of the allele
  changes the trait or the risk; and
- a **P-value**, which says how surprising an association that strong would be if the SNP had no effect.

{doc}`case-control-design` and {doc}`odds-ratios` show how these numbers are computed for one SNP.

## Heritability and GWAS

Heritability and GWAS answer two different questions about the same genetic variation:

- **Heritability** quantifies the **overall** genetic contribution to trait variance. It says *how much*
  genes matter, without saying which genes.
- **GWAS** discover and quantify the **specific** genetic factors, loci and variants, that drive heritability.
  The results can be used
  - **individually**, as variant–trait associations, each with its own effect size and P-value; and
  - **collectively**, to estimate SNP-heritability and to build polygenic risk scores (PRS).

You can think of the variance of a trait as a set of nested pieces:

```{figure} images/variance_explained.png
:alt: Four nested bars. The top bar is all phenotypic variance. Below it, the variance splits into a genetic part (heritability, from twins) and an environmental part. Below the genetic part, common SNPs capture SNP-heritability, and the rest is not tagged by SNP arrays. Below that, only part of the SNP-heritability comes from genome-wide significant GWAS loci; the rest comes from SNPs whose effects are below the significance threshold.
:width: 100%

How trait variance divides up, from all differences between people down to the loci a GWAS can name.
Bar widths are illustrative.
```

Each step down the figure is a smaller slice. Heritability is part of the phenotypic variance.
**SNP-heritability** ($h^2_\text{SNP}$) is the part of heritability captured by the common variants included in a
study; it can be estimated from GWAS data, as {doc}`snp-heritability` explains. The **genome-wide significant
loci** are the part of SNP-heritability whose effects are large enough, in a given study, to pass the significance
threshold. As GWAS grow larger, more small effects pass the threshold and the bottom bar grows toward the one
above it. In general,

$$
h^2_\text{SNP} \;\le\; h^2 \;\le\; H^2
$$

## GWAS and genetic architecture

GWAS are used to study the **genetic architecture** of heritable traits: how their genetic basis is organized.
For each trait, the trait-associated variants tell us:

- **how many** loci are involved (the **polygenicity** of the trait);
- **where** they are in the genome, which points to genes and biological pathways; and
- **how much** of the variance they explain, individually and together.

Large GWAS have shown that most common traits and diseases are highly polygenic. Type 2 diabetes, the running
example in this chapter, has more than 600 associated loci in a multi-ancestry GWAS of 2.5 million people
(Suzuki et al., [2024](https://doi.org/10.1038/s41586-024-07019-6)), and for adult height a GWAS of about 5 million
people identified over 12,000 independent associated SNPs (Yengo et al.,
[2022](https://doi.org/10.1038/s41586-022-05275-y)). Each of these variants changes the trait only slightly;
it is their combined effect that matters.

## From GWAS to polygenic risk

A GWAS tests every SNP and reports the result for each one. The complete table of results is called the GWAS
**summary statistics**: one row per SNP giving its position, its alleles, its effect size and standard error,
and its P-value ({doc}`association-models`, {doc}`significance-and-plots`). Summary statistics are usually shared
publicly, and they make GWAS results reusable.

In particular, they are the raw material for **polygenic risk**, which quantifies a person's overall genetic
liability for a trait or disease, in the sense of the liability-threshold model ({doc}`polygenic-model`). A
**polygenic risk score (PRS)** is the weighted sum of the risk alleles a person carries,

$$
\text{PRS} = \sum_j \beta_j \, G_j
$$

where $\beta_j$ is the GWAS effect size of SNP $j$ and $G_j$ is the number of risk alleles (0, 1 or 2) the person
carries at SNP $j$. It is the genetic score of {doc}`polygenic-model`, with each allele weighted by its effect.
Like everything else in common-disease genetics, a PRS is **probabilistic**: it shifts a person's chance of
disease without saying who will become ill.

Polygenic risk scores are the subject of the next chapter, {doc}`../polygenic-risk-scores/index`, which covers how
they are built, evaluated and used, and their limitations.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/heritability-to-gwas.json")
```
