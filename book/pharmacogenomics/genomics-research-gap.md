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

# The genomics research gap

```{admonition} Module at a glance
**Lecture:** Pharmacogenomic diversity &nbsp;·&nbsp; **Prerequisites:** {doc}`precision-public-health`, {doc}`../complex-traits-gwas/index`, {doc}`../polygenic-risk-scores/index`

**After this module you should be able to:** quantify the European bias of genomics research cohorts, and explain its consequences for genomic medicine, for polygenic risk prediction and for health disparities.
```

**Genomic medicine** is the use of genomic information about an individual as part of their clinical care,
in support of better diagnostic, prognostic and therapeutic decisions. It promises to change health care —
but only for the populations it has studied.

The vast majority of genomics research cohorts are made up of individuals with European ancestry. Clinical
insights based on the study of European-ancestry genomes will not necessarily replicate across diverse
populations.

## How skewed is the evidence base?

Popejoy and Fullerton audited the GWAS Catalog in 2016 and compared it with the same audit in 2009.

```{figure} images/diversity_gap.png
:alt: Stacked bars. In 2009, from 373 studies and 1.7 million samples, 96% of GWAS participants were of European ancestry and 4% were not. In 2016, from 2,511 studies and 35 million samples, 81% and 19%.
:width: 640px

GWAS participants by ancestry, 2009 and 2016. Data from Popejoy AB & Fullerton SM (2016) Genomics is
failing on diversity. *Nature* 538:161–164.
```

Seven years, a sevenfold increase in the number of studies, and a twentyfold increase in samples — and
European-ancestry participants fell only from 96% to 81%. That is real progress, but the absolute
shortfall grew, because the total grew so much faster.

The breakdown is the more damaging number.

```{figure} images/diversity_breakdown.png
:alt: Non-European GWAS participants by group as a share of all samples in 2009 and 2016. Asian 3% to 14%, African 0.57% to 3%, mixed 0% to 1%, Hispanic and Latin American 0.06% to 0.54%, Pacific Islander 0.15% to 0.28%, Native peoples 0.06% to 0.05%.
:width: 700px

Non-European participants as a share of all GWAS samples. Data from Popejoy & Fullerton (2016).
```

Almost all of the improvement between 2009 and 2016 came from Asian-ancestry cohorts. African-ancestry
individuals account for about 3% of samples, despite Africa holding more genetic diversity than any other
continent. Hispanic and Latin American populations — the subject of the rest of this chapter — account for
about half a percent. For Native peoples there was no increase at all.

## Consequences for polygenic risk prediction

This connects directly to {doc}`../polygenic-risk-scores/index`. A PRS is a weighted sum of GWAS effect
sizes, and both the tag SNPs and their linkage disequilibrium with the causal variants are
population-specific. Prediction accuracy therefore falls roughly with genetic distance from the discovery
cohort: scores derived in European samples perform worst in African-ancestry individuals — the populations
already least well served.

> "This disparity is an inescapable consequence of Eurocentric biases in GWAS…"
>
> — Martin AR et al. (2019) Clinical use of current polygenic risk scores may exacerbate health
> disparities. *Nat Genet* 51:584–591

Deploying current scores as they are would widen, not close, health disparities.

## Why it matters

The gap limits the reach of genomic medicine and threatens to exacerbate existing health disparities. The
promise of genomic medicine will not be fully realized until genomic studies become more broadly
representative of global populations.

> "Genome-based medicine is at the heart of a revolution in medical care. There is broad agreement on the
> need for genomic studies to be more representative of all populations so that people everywhere benefit
> from improved health outcomes."
>
> — Simon Gravel & Mark Lathrop, *Closing the genomics research gap*, McGill University, June 2018

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/genomics-research-gap.json")
```
