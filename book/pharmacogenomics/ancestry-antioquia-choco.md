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

# Genetic ancestry: Antioquia and Chocó

```{admonition} Module at a glance
**Lecture:** Pharmacogenomic diversity &nbsp;·&nbsp; **Prerequisite:** {doc}`colombia`

**After this module you should be able to:** characterise two admixed Colombian populations by genetic ancestry, and explain why admixture makes them informative for pharmacogenomics.
```

Colombian populations are **admixed**: three-way mixtures of African, European and Native American
ancestry, in proportions that differ sharply from place to place. The formal methods for estimating
ancestry come later, in {doc}`../genetic-diversity-pca/index` and
{doc}`../genetic-ancestry-admixture/index`. Here we only need to read the output.

## Ancestry proportions across populations

```{figure} images/admixture_populations.png
:alt: ADMIXTURE bar plot across Yoruba, Mende, British, Spanish, Zapotec, Embera, Quechua, Antioquia, Chocó and African American samples, showing African, European and Native American ancestry fractions.
:width: 100%

ADMIXTURE ancestry fractions across reference and Colombian populations. Each vertical bar is one
individual. Reproduced from [Nagar et al. (2019)](https://doi.org/10.3389/fgene.2019.00241),
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
```

Each vertical bar is one individual, and its colours are that person's estimated ancestry fractions. The
reference populations at the left and middle are close to unmixed. The two Colombian cohorts are not:
they sit between the references, and they sit in different places.

- **Antioquia** is roughly three-quarters European, with most of the remainder Native American.
- **Chocó** is roughly three-quarters African.

Two adjacent departments, two very different genomes. The history is not subtle: Chocó's population
descends largely from Africans brought to the Pacific lowlands under slavery, while Antioquia's
descends largely from Spanish settlers and Indigenous peoples of the interior highlands.

## Why admixture is useful here

Admixture is not just a complication to be controlled for — it is what makes these cohorts informative.
Within an admixed population, individuals differ in their ancestry fractions. So for any variant whose
frequency differs among the source populations, you can ask a question **within** one cohort: do the
people carrying the variant have more of one ancestry than the people who do not? That is how
{doc}`pharmasnp-variation` links particular pharmaSNPs to particular ancestries, and it is the same logic
as the admixture mapping you will meet in {doc}`../genetic-ancestry-admixture/index`.

It also sets up the pharmacogenomic prediction. Any variant whose frequency differs among African,
European and Native American source populations will differ between Antioquia and Chocó. The question is
whether the variants that differ are ones anybody prescribes against.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/ancestry-antioquia-choco.json")
```
