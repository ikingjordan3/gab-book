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

# Knowledge bases and guidelines

```{admonition} Module at a glance
**Lecture:** Pharmacogenomics &nbsp;·&nbsp; **Prerequisite:** {doc}`pharmagenes`

**After this module you should be able to:** find the curated variant–drug evidence and the prescribing guidelines, and read CPIC's standardized terminology for allele function and metabolizer phenotype.
```

Two resources carry almost all of clinical pharmacogenomics, and they do different jobs.

## ClinPGx: the evidence

[ClinPGx](https://www.clinpgx.org/) curates **variant–drug relationships**, each with a level of
evidence. It grew out of PharmGKB, a knowledge base of variant–drug annotations begun at Stanford, and
now integrates CPIC, PharmCAT and PharmVar in one resource. Every annotation is graded:

```{figure} images/evidence_levels.png
:alt: PharmGKB clinical annotation levels of evidence: levels 1a and 1b high, 2a and 2b moderate, level 3 low, level 4 preliminary.
:width: 560px

PharmGKB clinical annotation levels of evidence. Levels 1 and 2 are the clinically actionable subset.
```

These levels matter later: {doc}`pharmasnp-variation` uses the broad set (levels 1–4, n = 1,995 variants)
for genome-wide comparisons, and the strict, clinically actionable set (levels 1–2, n = 155) for the
Colombian comparison.

## CPIC: what to do about it

The barrier to pharmacogenomics in the clinic has rarely been the test. It is knowing **what to do with
the result**. The [Clinical Pharmacogenetics Implementation Consortium (CPIC)](https://cpicpgx.org/),
a joint project of PharmGKB and the NIH Pharmacogenomics Research Network begun in 2009, writes
peer-reviewed, freely available **prescribing guidelines** for gene–drug pairs.

CPIC guidelines deliberately do **not** tell a clinician whether to order a test. They say what to do
when a result already exists. That design decision is what made them usable: test results arrive from
many sources, and the clinician's question is always the same.

Each guideline follows one chain: **genotype → phenotype → therapeutic recommendation**, with a strength
attached to the recommendation.

## Standardized terminology

Before 2017 the same diplotype was reported as "extensive metabolizer" by one laboratory and "normal
metabolizer" by another — and "extensive" was routinely misread by clinicians as unusual rather than
typical. CPIC's standardized terms fixed this, and they are what laboratory reports and guideline tables
now share.

Allele function is assigned first:

| Final term | Functional definition | Example |
|---|---|---|
| Increased function | Function greater than normal function | *CYP2C19\*17* |
| Normal function | Fully functional / wild-type | *CYP2C19\*1* |
| Decreased function | Function less than normal function | *CYP2C19\*9* |
| No function | Nonfunctional | *CYP2C19\*2* |
| Unknown function | No literature describing function, or the allele is novel | *CYP2C19\*29* |
| Uncertain function | Literature supporting function is conflicting or weak | *CYP2C19\*12* |

The diplotype then determines the phenotype:

| Final term | Functional definition | Genetic definition | Example |
|---|---|---|---|
| Ultrarapid metabolizer | Increased enzyme activity compared to rapid metabolizers | Two increased function alleles, or more than two normal function alleles | *CYP2C19\*17/\*17* |
| Rapid metabolizer | Increased activity compared to normal, less than ultrarapid | Combinations of normal and increased function alleles | *CYP2C19\*1/\*17* |
| Normal metabolizer | Fully functional enzyme activity | Combinations of normal and decreased function alleles | *CYP2C19\*1/\*1* |
| Intermediate metabolizer | Decreased activity, between normal and poor metabolizer | Combinations of normal, decreased and/or no function alleles | *CYP2C19\*1/\*2* |
| Poor metabolizer | Little to no enzyme activity | Combination of no function and/or decreased function alleles | *CYP2C19\*2/\*2* |

Source: Caudle KE et al. (2017) Standardizing terms for clinical pharmacogenetic test results.
*Genet Med* 19:215–223.

The same five phenotype terms are used for **every** drug-metabolizing enzyme, which is what lets a
laboratory report and a guideline table line up. Transporters use the function terms rather than the
metabolizer terms, and HLA results are reported simply as positive or negative.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/knowledge-bases.json")
```
