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
**Lecture:** Genome-wide association studies &nbsp;·&nbsp; **Prerequisite:** {doc}`gwas-to-prs`

**After this module you should be able to:** connect heritability to the specific loci found by GWAS, describe what GWAS reveal about the genetic architecture of a trait, and explain what GWAS summary statistics are.
```

This module opens the second lecture of the chapter. The first lecture ended with the idea of a GWAS
({doc}`gwas-to-prs`); the rest of the chapter works through how a GWAS is actually done, from study design
to reading the results.

## Heritability and GWAS

Heritability and GWAS answer two different questions about the same genetic variation:

- **Heritability** quantifies the **overall** genetic contribution to trait variance. It says *how much*
  genes matter, without saying which genes.
- **GWAS** discover and quantify the **specific** genetic factors, loci and variants, that drive heritability.
  The results can be used
  - **individually**, as variant–trait associations, each with its own effect size and P-value; and
  - **collectively**, to build polygenic risk scores (PRS) and to estimate SNP-heritability.

You can think of the variance of a trait as a set of nested pieces:

```{figure} images/variance_explained.png
:alt: Four nested bars. The top bar is all phenotypic variance. Below it, the variance splits into a genetic part (heritability, from twins) and an environmental part. Below the genetic part, common SNPs capture SNP-heritability, and the rest is not tagged by SNP arrays. Below that, only part of the SNP-heritability comes from genome-wide significant GWAS loci; the rest comes from SNPs whose effects are below the significance threshold.
:width: 100%

How trait variance divides up, from all differences between people down to the loci a GWAS can name.
Bar widths are illustrative.
```

Each step down the figure is a smaller slice. Heritability is part of the phenotypic variance; SNP-heritability
is the part of heritability captured by common SNPs; and the **genome-wide significant loci** are the part of
SNP-heritability whose effects are large enough, in a given study, to pass the significance threshold. As GWAS
grow larger, more small effects pass the threshold and the bottom bar grows toward the one above it.

## GWAS and genetic architecture

GWAS are used to study the **genetic architecture** of heritable traits: how their genetic basis is organized.
For each trait, the trait-associated variants tell us:

- **how many** loci are involved (the **polygenicity** of the trait);
- **where** they are in the genome, which points to genes and biological pathways; and
- **how much** of the variance they explain, individually and together.

Large GWAS have shown that most common traits and diseases are highly polygenic. Type 2 diabetes, the running
example in this lecture, has more than 600 associated loci in a multi-ancestry GWAS of 2.5 million people
(Suzuki et al., [2024](https://doi.org/10.1038/s41586-024-07019-6)), and for adult height a GWAS of about 5 million
people identified over 12,000 independent associated SNPs (Yengo et al.,
[2022](https://doi.org/10.1038/s41586-022-05275-y)). Each of these variants changes the trait only slightly;
it is their combined effect that matters.

## From GWAS to polygenic risk: summary statistics

A GWAS tests every SNP and reports the result for each one. The complete table of results is called the GWAS
**summary statistics**: a summary of the statistical association tests for all loci genome-wide, with one row
per SNP giving its position, its alleles, its effect size and standard error, and its P-value.
{doc}`significance-and-plots` shows what such a table looks like and how it is drawn as a Manhattan plot.

Summary statistics are what make GWAS results reusable. They are usually shared publicly, and they are the raw
material for polygenic risk scores:

- a PRS provides a single **point estimate of the genetic risk** (liability) for a trait or disease that a GWAS
  has studied; and
- it can be used to make predictions about the trait from a **personal genome sequence** or genotype data, by
  weighting each allele the person carries by its effect size from the summary statistics.

{doc}`../polygenic-risk-scores/index` builds on this in detail.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/heritability-to-gwas.json")
```
