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

# From variants to assays

```{admonition} Module at a glance
**Lecture:** Pharmacogenomic diversity &nbsp;·&nbsp; **Prerequisite:** {doc}`pharmasnp-variation`

**After this module you should be able to:** explain how a population-guided variant list becomes an affordable clinical test, and how such a test is validated.
```

The practical payoff of precision public health is this: once you know **which** variants matter in a
given population, you no longer need to sequence anybody.

## A targeted, low-cost assay

```{figure} images/pcr_assay.png
:alt: Design of an allele-specific PCR assay for rs4149056 in SLCO1B1, with effect and non-effect allele forward primers and a shared reverse primer, and gel images for the three genotypes.
:width: 100%

An allele-specific PCR assay for rs4149056 (*SLCO1B1\*5*). Reproduced from
[Nagar et al. (2019)](https://doi.org/10.3389/fgene.2019.00241),
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
```

The design is ordinary molecular biology. Two forward primers differ at their 3′ end — one matching the
effect allele, one the non-effect allele — and share a reverse primer. Whether a band appears tells you
which allele is present. Any clinical laboratory with a thermocycler and a gel rig can run it, for pennies
per reaction.

## Validation

A cheap assay is worth nothing unless its calls are right. The assay was checked two ways against the same
individuals: against aligned sequencing reads at the variant position, and against exome genotypes.

```{figure} images/pcr_concordance.png
:alt: Three confusion matrices comparing PCR genotypes with exome genotypes. rs4149056 simvastatin toxicity 97.7% accurate, rs1799853 warfarin dosage 98.5%, rs1057910 warfarin dosage 100%.
:width: 100%

Concordance between allele-specific PCR calls and exome genotypes, for three pharmaSNPs. Reproduced from
[Nagar et al. (2019)](https://doi.org/10.3389/fgene.2019.00241),
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
```

Accuracy of 97.7% to 100% against exome genotypes — accurate enough to use, and cheap enough to deploy at
scale. That combination is the whole argument.

## Precision public health, concretely

Put the chapter back together:

1. Pharmacogenomic variants are **common** and have **large effects** on drug response
   ({doc}`pharmacogenomics-defined`).
2. Because they are common, their frequencies **differ among populations**
   ({doc}`pharmasnp-variation`).
3. So the variants worth testing for differ among populations too — and in a resource-limited setting you
   can test for exactly those, cheaply, instead of sequencing everyone.

That is what population-guided pharmacogenomics means: allocate the effort where it will do the most good,
by finding the variants with special relevance to a specific population.

```{admonition} The argument, and its limits
:class: important
Closing the genomics research gap, and adopting precision public health, can help to ensure that genomic
medicine works for all people everywhere.

Be clear about what has been shown and what is being argued. The **evidence** in this chapter is that
pharmacogenomic allele frequencies differ among populations, that European-derived polygenic scores lose
accuracy elsewhere, and that a targeted assay can be both accurate and cheap. That a broader evidence base
would improve care everywhere **follows** from those findings — but it is a claim about what we should do,
and it is worth being able to separate the two.
```

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/from-variants-to-assays.json")
```
