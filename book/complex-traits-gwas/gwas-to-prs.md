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

# GWAS, SNP-heritability and polygenic risk

```{admonition} Module at a glance
**Lecture:** Genetics of complex, common disease &nbsp;·&nbsp; **Prerequisite:** {doc}`heritability`

**After this module you should be able to:** describe what a GWAS looks for, and show how GWAS results feed two summary measures: SNP-heritability for a trait and a polygenic risk score for a person.
```

The previous modules established that complex traits are heritable and that their genetic basis is many
variants of small effect. This module closes the first lecture by previewing the tool that finds those
variants, the **genome-wide association study (GWAS)**, and the two things GWAS results are used for. The
GWAS modules that follow, starting with {doc}`heritability-to-gwas`, explain how a GWAS is done.

## How do we study complex traits that are heritable?

**Genome-wide association studies (GWAS)** seek to identify the single nucleotide polymorphisms (SNPs) that
are common in the human genome and contribute to the variance in observed phenotypes. Very simply, this means
finding **SNPs that vary in frequency between cases and controls**.

At each of hundreds of thousands to millions of SNPs, a GWAS asks whether one allele is more common in people
with the disease (cases) than in people without it (controls). For each SNP it reports:

- an **effect size** ($\beta$, or an odds ratio for a disease), which says how much each copy of the allele
  changes the trait or the risk; and
- a **P-value**, which says how surprising an association that strong would be if the SNP had no effect.

{doc}`case-control-design` and {doc}`odds-ratios` show how these numbers are computed for one SNP.

## GWAS and heritability

**Heritability** quantifies the overall genetic contribution to trait variance. **GWAS** discover and quantify
the specific genetic factors that drive heritability, as SNP effect sizes and P-values. GWAS data can therefore
be used to estimate **SNP-heritability** ($h^2_\text{SNP}$): the proportion of trait variance explained by the
common variants included in the study.

SNP-heritability is estimated in the way described in {doc}`heritability`, by relating genetic similarity to
trait similarity, but with genetic similarity measured from genotypes rather than inferred from family
relationships. Two widely used approaches are:

- **GREML** (implemented in the software GCTA), which estimates how much of the trait variance is explained by
  the genome-wide genetic similarity between pairs of unrelated, genotyped individuals; and
- **LD score regression**, which estimates $h^2_\text{SNP}$ from GWAS summary statistics alone, using the fact
  that SNPs in high linkage disequilibrium with many other SNPs tag more causal variation.

**Missing heritability.** For most traits, SNP-heritability is lower than twin heritability, and early GWAS hits
explained far less still. For height, twin studies give $h^2 \approx 0.8$; the first GWAS-significant SNPs
explained only about 5% of the variance, while all common SNPs together explained about 45% (Yang et al.,
[2010](https://doi.org/10.1038/ng.608)). The gap between heritability and what GWAS explain became known as
**missing heritability**. Part of it is closed by larger GWAS, which detect more small effects; part reflects
rare variants that arrays do not tag, and part may reflect overestimation of twin heritability. This is why

$$
h^2_\text{SNP} \;\le\; h^2 \;\le\; H^2
$$

## GWAS and polygenic risk

**Polygenic risk** quantifies a person's overall genetic liability for a trait or disease, in the sense of the
liability-threshold model ({doc}`polygenic-model`). GWAS results can be used to calculate it as the
**weighted sum of risk alleles** in that person's genome, the **polygenic risk score (PRS)**:

$$
\text{PRS} = \sum_j \beta_j \, G_j
$$

where

- $\beta_j$ is the GWAS effect size of SNP $j$ (for a disease, usually the log odds ratio), and
- $G_j$ is the number of risk alleles (0, 1 or 2) that the person carries at SNP $j$.

This is the genetic score of {doc}`polygenic-model`, with each allele weighted by its estimated effect instead
of counted equally.

**A worked example.** The table uses six well-known type 2 diabetes SNPs. The effect sizes are illustrative,
the same values used in the summary-statistics table in {doc}`significance-and-plots`.

| SNP | Gene | Risk allele | $\beta$ (log OR) | Person 1: $G$ | Person 2: $G$ |
|---|---|---|---|---|---|
| rs7903146 | *TCF7L2* | T | 0.315 | 2 | 0 |
| rs7756992 | *CDKAL1* | G | 0.157 | 1 | 1 |
| rs9939609 | *FTO* | A | 0.122 | 1 | 0 |
| rs4402960 | *IGF2BP2* | T | 0.113 | 0 | 1 |
| rs13266634 | *SLC30A8* | C | 0.113 | 2 | 1 |
| rs1801282 | *PPARG* | C | 0.122 | 2 | 1 |

```{code-cell} ipython3
beta = {"rs7903146": 0.315, "rs7756992": 0.157, "rs9939609": 0.122,
        "rs4402960": 0.113, "rs13266634": 0.113, "rs1801282": 0.122}
person1 = {"rs7903146": 2, "rs7756992": 1, "rs9939609": 1,
           "rs4402960": 0, "rs13266634": 2, "rs1801282": 2}
person2 = {"rs7903146": 0, "rs7756992": 1, "rs9939609": 0,
           "rs4402960": 1, "rs13266634": 1, "rs1801282": 1}

def prs(genotypes):
    return sum(beta[snp] * g for snp, g in genotypes.items())

for name, g in [("Person 1", person1), ("Person 2", person2)]:
    print(f"{name}: PRS = {prs(g):.3f}")
```

Person 1 has the higher score, mainly because they carry two copies of the *TCF7L2* risk allele, the SNP with
the largest effect. The difference in their scores, about 0.87 on the log-odds scale, corresponds to roughly
2.4 times higher odds of disease (e<sup>0.87</sup> ≈ 2.4), if the effects simply add up.

**Many SNPs, a bell curve.** Real PRS combine thousands to millions of SNPs, each of tiny effect. Summed across
so many loci, the scores in a population follow the normal distribution predicted by the central limit theorem
({doc}`polygenic-model`). The simulation below gives 10,000 people genotypes at 1,000 SNPs with small random
effect sizes:

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt

rng = np.random.default_rng(1)
n_people, n_snps = 10_000, 1_000
freq = rng.uniform(0.05, 0.5, n_snps)            # risk-allele frequencies
beta = rng.normal(0, 0.02, n_snps)               # small effects (log OR)
G = rng.binomial(2, freq, size=(n_people, n_snps))
scores = G @ beta                                # PRS = sum of beta_j * G_j

fig, ax = plt.subplots(figsize=(6.5, 3.2))
ax.hist(scores, bins=60, color="#3B5B92")
ax.set_xlabel("Polygenic risk score"); ax.set_ylabel("People")
ax.spines[["top", "right"]].set_visible(False)
plt.show()
```

A PRS places a person on this distribution, which is to say on the genetic part of the liability scale. Like
everything else in common-disease genetics, it is **probabilistic**: a high score raises the chance of disease
but does not say who will become ill. {doc}`../polygenic-risk-scores/index` covers how PRS are built, evaluated and
used, and their limitations, including their reduced accuracy in populations under-represented in GWAS.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/gwas-to-prs.json")
```
