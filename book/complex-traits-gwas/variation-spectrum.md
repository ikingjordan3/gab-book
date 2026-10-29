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

# Genetic variation in disease

```{admonition} Module at a glance
**Lecture:** Genetics of complex, common disease &nbsp;·&nbsp; **Prerequisite:** {doc}`../monogenic-disease/index`

**After this module you should be able to:** contrast rare, monogenic disease with common, polygenic disease in terms of allele frequency, effect size and penetrance, and explain why genetic risk for common disease is probabilistic.
```

Suppose you have re-sequenced a study participant or a patient. You have accessed the data, checked the
quality of the reads, mapped them to the reference genome and called variants. The last step is to
**interpret** those variants: what do they mean for this person's health?

For a rare, monogenic disease the answer is often a **diagnosis**: one or two variants in a single gene
explain the condition. For a common disease such as type 2 diabetes, the
answer is a **risk estimate**. Risk comes from the combined effect of many variants across the genome,
each with a small effect, so the tools change too: genome-wide association studies (GWAS) find the
variants, and polygenic risk scores (PRS) add them up for one person.

## Three kinds of genetic variation in health

| | Rare, monogenic (Mendelian) disease | Common, complex, polygenic disease | Response to drugs |
|---|---|---|---|
| **Genetic cause** | One or a few loci | Many loci plus environmental effects | Often one or a few loci |
| **Allele frequency** | Rare | Common | Common |
| **Penetrance** | High | Low | High |
| **Covered in** | {doc}`chapter 5 <../monogenic-disease/index>` | this chapter | {doc}`chapter 10 <../pharmacogenomics/index>` |

**Penetrance** is the proportion of people carrying a variant (or genotype) who develop the trait. A
highly penetrant variant almost always causes disease. A variant with low penetrance raises risk only
slightly, so most carriers never develop the disease.

## Allele frequency and effect size

A useful way to organize disease-associated variants is by two numbers: how common the allele is in the
population (**allele frequency**) and how much it changes disease risk (**effect size**, here an odds
ratio, introduced in {doc}`odds-ratios`).

```{figure} images/spectrum_plain.png
:alt: Allele frequency (very rare to common) on the x-axis against effect size (low to high odds ratio) on the y-axis. Rare alleles causing Mendelian disease sit at top left (rare, high effect). Common variants found by GWAS sit at bottom right (common, low to modest effect). Low-frequency variants with intermediate effect sit in the middle. Rare variants of small effect at bottom left are very hard to identify, and there are few examples of high-effect common variants at top right.
:width: 680px

Disease-associated variants by allele frequency and effect size. Most known variants fall in the
shaded diagonal band. Redrawn after Manolio et al. (2009).
```

Most variants that have been linked to disease fall along a diagonal:

- **Rare alleles causing Mendelian disease** (top left) are rare but have large effects. Family studies
  found them first, because a variant that almost always causes disease is easy to follow through a
  pedigree.
- **Common variants implicated by GWAS** (bottom right) are common but each raises risk only a little,
  typically with odds ratios between about 1.1 and 1.5.
- **Low-frequency variants with intermediate effect** sit in between.

Two corners are nearly empty, for different reasons:

- **Few common variants of high effect.** An allele that strongly increases the risk of a serious disease
  tends to be removed by natural selection, so it rarely becomes common.
- **Rare variants of small effect** are very hard to find. They are too rare to test reliably in a
  population study and their effect is too small to see in a family.

## One gene versus many

The difference between the two architectures shows up clearly when you look at a whole genome. Everyone
carries millions of variants, and almost all of them have no effect on any given disease.

```{figure} images/karyotype_cf.png
:alt: Schematic genome of person one, with cystic fibrosis. Every chromosome carries many variants with no effect on the disease. Two disease-causing variants are highlighted, one on each copy of chromosome 7.
:width: 680px

Person one has cystic fibrosis. The disease is explained by two variants in one gene, *CFTR* on
chromosome 7, one on each copy. Variant positions are schematic.
```

**Cystic fibrosis** is an autosomal recessive disease: a person is affected when both copies of *CFTR*
carry a disease-causing variant. If you find those two variants, you have the diagnosis. In that
sense, the genetics is **deterministic**.

```{figure} images/karyotype_t2d.png
:alt: Schematic genome of person two, with type 2 diabetes. The same background variants are shown, plus many small-effect risk variants spread across most chromosomes.
:width: 680px

Person two has type 2 diabetes. Risk comes from many small-effect variants spread across the
genome. Variant positions are schematic.
```

**Type 2 diabetes** is a common, complex disease. GWAS have found more than 600 associated loci
(Suzuki et al., [2024](https://doi.org/10.1038/s41586-024-07019-6)), each of which changes risk only slightly, and obesity, diet, physical activity and age matter too. No single
variant decides whether a person develops the disease. Carrying more risk alleles
shifts the odds, so the genetics is **probabilistic**.

## Rare disease in clinical practice

The clinical path for rare, monogenic disease is well established (see {doc}`../clinical-genetic-testing/index`):

1. **Discovery.** Linkage studies in families, and later sequencing designs such as exome sequencing of
   an affected child and both parents (trios), identify genes and variants that cause disease.
2. **Databases.** Those findings are compiled in public resources. **OMIM** catalogues human genes and
   genetic disorders, **ClinVar** archives the classifications that laboratories submit for individual
   variants, and **ClinGen** curates which genes cause which diseases and how strong the evidence is.
3. **Interpretation.** The American College of Medical Genetics and Genomics (**ACMG**) sets standards of
   evidence for linking a patient's variant to disease. Variants are classified as pathogenic, likely
   pathogenic, of uncertain significance, likely benign or benign.

There is no equivalent of this pipeline for common disease. A single common variant with an odds ratio of
1.2 is not "pathogenic"; it is one small contribution to risk among many.

## Changing our minds: from diagnosis to risk

Moving from rare to common disease means moving from a deterministic to a probabilistic way of thinking:

- Genetic associations with disease, found by GWAS, are **probabilistic**: an associated allele is more
  common in cases than in controls, but most carriers are healthy.
- Estimates of disease risk, such as a PRS, are **probabilistic**: they shift a person's chance of disease
  up or down, but they do not say who will become ill.

The rest of this chapter builds the ideas needed to work with that kind of genetics:
{doc}`polygenic-model` explains how many small effects combine, {doc}`heritability` measures how much
genetics matters overall, and the GWAS modules show how the individual variants are found,
using type 2 diabetes as their running example.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/variation-spectrum.json")
```
