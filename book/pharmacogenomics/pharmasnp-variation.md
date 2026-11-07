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

# PharmaSNP variation

```{admonition} Module at a glance
**Lecture:** Pharmacogenomic diversity &nbsp;·&nbsp; **Prerequisite:** {doc}`ancestry-antioquia-choco`

**After this module you should be able to:** interpret $F_{ST}$, explain why pharmaSNPs are unusually variable and unusually differentiated among populations, and identify the Colombian frequency differences that matter clinically.
```

## How different are two populations?

To compare populations we need a number for how much of the genetic variation at a locus lies **between**
populations rather than within them. That number is the **fixation index**, $F_{ST}$.

```{figure} images/fst_intuition.png
:alt: Two scenarios at one SNP. When the two populations have similar allele frequencies — 0.5 and 0.4 — F_ST is 0.01 and the variation is mostly within populations. When they differ — 0.9 and 0.2 — F_ST is 0.49 and much of the variation is between them.
:width: 100%

$F_{ST}$ at one SNP, in two scenarios. The values are computed from the allele frequencies shown.
```

The intuition is all we need here:

- If two populations have **similar** allele frequencies, knowing which population someone came from tells
  you little about their genotype. $F_{ST}$ is near 0.
- If the frequencies **differ sharply**, population membership is informative. $F_{ST}$ approaches 1.

Human populations are mostly in the first regime: typical genome-wide $F_{ST}$ between continental groups
is on the order of 0.1. That is what makes the next result notable. The estimator itself, its variance
components, and its relationship to principal components analysis come in
{doc}`../genetic-diversity-pca/index`.

## PharmaSNPs are unusually variable and unusually differentiated

```{figure} images/pharmasnp_variation.png
:alt: Two bar charts. Average minor allele frequency is 0.02 for genome-wide SNPs and 0.25 for pharmaSNPs. Average fixation index F_ST is 0.01 for genome-wide SNPs and 0.07 for pharmaSNPs.
:width: 100%

PharmaSNPs compared with genome-wide SNPs across 26 populations. Redrawn from the values reported by
[Nagar et al. (2019)](https://doi.org/10.3389/fgene.2019.00241): avg. MAF 0.25 vs 0.02, avg. $F_{ST}$
0.07 vs 0.01; *n* = 1,995 pharmaSNPs and *n* = 28,137,656 genome-wide SNPs.
```

PharmaSNPs are more than ten times as common as the average genome-wide SNP, and several times more
differentiated between populations.

```{admonition} Read this result carefully
:class: warning
Part of the frequency difference is **ascertainment**. A variant gets a clinical annotation because
somebody studied it, and common variants are easier to study, so the pharmaSNP set is enriched for common
variants by construction. That does not explain away the $F_{ST}$ result, though: conditional on being
common, these variants are still more differentiated among populations than the genome-wide background.
Pharmacogenomic genotypes really do differ systematically between populations.
```

## PharmaSNPs recapitulate population structure

```{figure} images/mds_panels.png
:alt: Two multidimensional scaling plots, one from genome-wide SNPs and one from pharmaSNPs, both separating African, European, Native American and Colombian samples into the same arrangement.
:width: 100%

Multidimensional scaling of genetic distances, from genome-wide SNPs (left) and from pharmaSNPs alone
(right). Each point is one individual. Reproduced from
[Nagar et al. (2019)](https://doi.org/10.3389/fgene.2019.00241),
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
```

Fewer than 2,000 pharmaSNPs reproduce the same population structure as 28 million genome-wide SNPs. Two
consequences follow. Pharmacogenomic variation is **structured by ancestry** — which is why a prescribing
guideline tuned to one population may misfire in another. And a small panel of pharmaSNPs carries real
ancestry information, which matters for privacy as much as for interpretation.

## Which variants differ in Colombia?

Now the clinically actionable set: the 155 pharmaSNPs at PharmGKB evidence levels 1–2.

```{figure} images/colombia_pharmasnps.png
:alt: Effect allele frequency in Antioquia and Chocó for twelve ancestry-associated pharmaSNPs. Variants such as rs7793837 (salbutamol efficacy) and rs11615 (platin toxicity) are far more frequent in Antioquia; rs2740574 (tacrolimus metabolism), rs6977820 (antipsychotic drug toxicity) and rs11212617 (metformin efficacy) are far more frequent in Chocó.
:width: 100%

Effect allele frequencies for ancestry-associated pharmaSNPs in the two Colombian populations. Redrawn
from Table 1 of [Nagar et al. (2019)](https://doi.org/10.3389/fgene.2019.00241).
```

The differences are large and they run in both directions. Some are drugs prescribed every day:

- **rs4149056** (*SLCO1B1*, simvastatin toxicity) is more frequent in Antioquia. This is the transporter
  variant from {doc}`pharmagenes` that raises plasma statin concentrations and the risk of myopathy.
- **rs11212617** (metformin efficacy) is more frequent in Chocó. Metformin is first-line treatment for
  type 2 diabetes, the running example of {doc}`../complex-traits-gwas/index` and
  {doc}`../polygenic-risk-scores/index`.
- **rs2740574** (tacrolimus metabolism) is more frequent in Chocó. Tacrolimus dosing matters for
  transplant patients, where errors are dangerous.

## The differences track ancestry

These frequency differences are not an accident of sampling two towns. Within the admixed cohorts, the
genotype tracks the ancestry fraction of the individual — the admixture logic from
{doc}`ancestry-antioquia-choco`.

```{figure} images/ancestry_rs4149056.png
:alt: Panels for rs4149056, simvastatin toxicity. Genotype frequencies in Antioquia and Chocó, and the association of genotype with African, European and Native American ancestry fractions.
:width: 100%

rs4149056 (*SLCO1B1*, simvastatin toxicity) tracks European ancestry within the admixed cohorts.
Reproduced from [Nagar et al. (2019)](https://doi.org/10.3389/fgene.2019.00241),
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
```

```{figure} images/ancestry_rs11212617.png
:alt: Panels for rs11212617, metformin efficacy. Genotype frequencies in Antioquia and Chocó, and the association of genotype with African, European and Native American ancestry fractions.
:width: 100%

rs11212617 (metformin efficacy) tracks African ancestry within the admixed cohorts. Reproduced from
[Nagar et al. (2019)](https://doi.org/10.3389/fgene.2019.00241),
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
```

Consider what this means for a national treatment guideline. A single recommendation for Colombia would be
wrong for one of these two departments — and the information needed to do better is a population allele
frequency, which is cheap to measure once.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/pharmasnp-variation.json")
```
