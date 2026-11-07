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

# Precision medicine and precision public health

```{admonition} Module at a glance
**Lecture:** Pharmacogenomic diversity &nbsp;·&nbsp; **Prerequisite:** {doc}`pharmacogenomics-defined`

**After this module you should be able to:** explain why the unit of intervention changes when genomics is applied at population scale, and what precision public health asks for instead.
```

The first half of this chapter treated pharmacogenomics one patient at a time. Pharmacogenomic variants
are common, so they also have **population-level allele frequencies** worth knowing — unlike the rare
variants of Mendelian disease, where only the individual genotype matters.

## Precision medicine

**Precision medicine** is an emerging approach for disease treatment and prevention that takes into
account individual variability in genes, environment and lifestyle for **each person**.

```{figure} images/precision_medicine.png
:alt: Three patients, each with their own genotype, each matched to a different treatment: drug A, drug B, or a lower dose.
:width: 100%

Precision medicine: the right treatment to the right patient at the right time.
```

The formulation is Collins and Varmus's, from the launch of the US Precision Medicine Initiative
(Collins FS & Varmus H, 2015, *N Engl J Med* 372:793–795). The slogan is useful, but notice what it
assumes: that every patient can be genotyped individually, and that several alternative drugs can be
stocked and chosen between.

## Precision public health

Those assumptions do not hold in much of the world. Sequencing every patient is out of reach for most
health systems, and so is maintaining alternatives for every prescription.

**Precision public health** keeps the genomic input but changes the unit of intervention from individuals
to **populations** (Khoury MJ, Iademarco MF & Riley WT, 2016, *Am J Prev Med* 50:398–401).

```{figure} images/precision_public_health.png
:alt: Three populations, each with its own allele frequencies, each matched to a different regimen.
:width: 100%

Precision public health: the right intervention to the right population at the right time.
```

| | Precision medicine | Precision public health |
|---|---|---|
| Unit of intervention | the individual | the population |
| Genomic input | individual genotype | population allele frequencies |
| Cost | prohibitive for many health systems | the same technologies, applied once per population |

The argument is not that individual care is wrong. It is that population-level genomic information is
affordable **now** and can be acted on **now**.

```{admonition} The efavirenz case, reconsidered
:class: important
{doc}`cyp2b6-efavirenz` ended with Zimbabwe's switch to an efavirenz-containing regimen. Read it again
with this module's vocabulary: allocating a drug regimen to a population without knowing its
pharmacogenomic allele frequencies is a precision public health failure. The information needed was a
population allele frequency, not 14 million individual genotypes.
```

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/precision-public-health.json")
```
