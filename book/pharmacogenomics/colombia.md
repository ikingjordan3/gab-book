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

# Population pharmacogenomics in Colombia

```{admonition} Module at a glance
**Lecture:** Pharmacogenomic diversity &nbsp;·&nbsp; **Prerequisite:** {doc}`genomics-research-gap`

**After this module you should be able to:** describe the design of a population pharmacogenomics study, and explain which evidence levels it rests on and why.
```

The rest of this chapter follows one study, carried out by our own lab with collaborators in Colombia
([Nagar et al. 2019](https://doi.org/10.3389/fgene.2019.00241)). It is a worked example of the precision
public health argument in {doc}`precision-public-health`, and of the kind of work that closes the gap in
{doc}`genomics-research-gap`.

The question is deliberately narrow and answerable:

> Do pharmacogenomic variants differ in frequency between Colombian populations — enough to matter for
> prescribing?

## Study design

```{figure} images/pharmasnp_workflow.png
:alt: Three-step workflow. One, mine pharmaSNPs — variants annotated with patient drug responses: efficacy, dosage and toxicity. Two, genotype populations — worldwide populations plus Antioquia and Chocó in Colombia. Three, compare frequencies — pharmaSNP frequencies within and between populations.
:width: 100%

The analysis in three steps.
```

1. **Mine pharmaSNPs.** Pull from PharmGKB every variant with a curated clinical annotation linking it to
   a drug and to an aspect of response — efficacy, dosage or toxicity.
2. **Genotype populations.** Use whole-genome data for worldwide reference populations, plus two Colombian
   populations: Antioquia and Chocó.
3. **Compare frequencies.** Ask how pharmaSNP allele frequencies vary within and between populations, and
   whether they behave differently from genome-wide SNPs.

## Two SNP sets

The study uses the PharmGKB evidence levels from {doc}`knowledge-bases` to define two sets, and the
distinction matters when reading the results:

| Set | Levels | Count | Used for |
|---|---|---|---|
| Broad | 1–4 | 1,995 pharmaSNPs | the worldwide comparison, where numbers matter for genome-wide statistics |
| Strict | 1–2 | 155 pharmaSNPs | the Colombian comparison — variants a clinician could act on today |

Compared against **28,137,656** pruned genome-wide SNPs as a background.

The broad set is large enough to support genome-wide comparisons of population structure. The strict set
is the clinically actionable one, which is why *n* drops to 155 when the question becomes "would this
change a prescription?"

```{admonition} Why two Colombian populations?
:class: note
Antioquia and Chocó are adjacent departments in north-west Colombia with sharply different colonial
histories, and — as {doc}`ancestry-antioquia-choco` shows — sharply different genomes. If pharmacogenomic
allele frequencies differ even between neighbouring departments of one country, then a single national
prescribing guideline cannot be right for both.
```

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/colombia.json")
```
