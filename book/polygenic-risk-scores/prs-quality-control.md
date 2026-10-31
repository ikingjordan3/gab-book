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

# Quality control

```{admonition} Module at a glance
**Lecture:** Polygenic risk scores 1 &nbsp;·&nbsp; **Prerequisite:** {doc}`base-target-data`

**After this module you should be able to:** list the quality-control checks base and target data need before a PRS is calculated, and why each matters.
```

A PRS depends on the **quality of the base and target data**. Both data sets get the same quality control (QC) as
in a GWAS — genotyping rate, minor allele frequency, Hardy–Weinberg equilibrium, imputation quality — plus some
checks that matter specifically for PRS. The QC procedures in the Choi et al. tutorial
([2020](https://doi.org/10.1038/s41596-020-0353-1)) are deliberately **conservative**: a PRS sums over thousands to
millions of variants, so small errors that would not matter for a single SNP can be amplified across the whole score.

## Base data

**Heritability check.** The predictive power of a PRS is limited by heritability, the fraction of trait variance
explained by genetic variance ({doc}`../complex-traits-gwas/heritability`). A PRS can never explain more than the
SNP-heritability of the trait ({doc}`../complex-traits-gwas/snp-heritability`); if heritability is low, PRS
predictive power will be low.

- Recommendation: build a PRS only from GWAS where **h²<sub>SNP</sub> > 0.05**.
- If no h²<sub>SNP</sub> is reported, you can estimate it yourself, for example with LD score regression. This will
  almost never be needed for large GWAS or validated scores.

**Effect allele.** Be sure that the effect allele is specified for the GWAS summary statistics. If the wrong
effect allele is assumed, the **direction of the effect will be opposite**. For example, at *SLC30A8* rs13266634 the
alternate allele T has OR = 0.89 and the risk allele C has OR = 1.12 — the same SNP, reported two ways
({doc}`base-target-data`). Read the wrong way round, every protective allele is counted as a risk allele and the
score is reversed.

## Target data

- QC can be **less stringent** when the target sample is small: removing people or SNPs costs power.
- Small target samples can also lead to **underpowered** PRS–trait association tests.
- Recommendation: calculate PRS in target samples of **at least 100** people.

## Base and target data

| Check | Why |
|---|---|
| Make sure files are not corrupted after transfer (compare `md5sum` checksums). | A truncated or altered file silently drops or changes variants. |
| Make sure base and target data are on the **same genome build**, or convert positions with liftOver. | On different builds (e.g. GRCh37 and GRCh38) the same SNP has different positions, so variants fail to match. |
| Remove **ambiguous** variants (A/T and C/G SNPs) and **duplicate** variants. | For A/T and C/G SNPs the two strands carry the same pair of alleles, so the strand — and the effect allele — cannot be checked. Duplicates cannot be matched uniquely. |
| Correct (or remove) **mismatching** variants. | Alleles reported on opposite strands (e.g. A/G in one file, T/C in the other) must be flipped before the weights are applied. |
| Remove individuals whose **reported sex** does not match their genetic sex (checked with PLINK). | A mismatch signals a sample mix-up. Sex also matters for type 2 diabetes, where it is a covariate ({doc}`../complex-traits-gwas/covariates`). |
| Make sure there is **no sample overlap** between base and target data. | People in both data sets helped estimate the weights, so the PRS over-fits them and looks more predictive than it is. |
| Make sure there is no **high degree of relatedness** (first- or second-degree relatives) between base and target data. | Relatives share genotypes and environments, which also inflates apparent performance. |

Checking for overlap and relatedness needs individual-level data from both data sets, which is rarely available for
a published GWAS. In practice, if base and target data come from different cohorts or populations, it is a safe
assumption that they do not overlap.

## Spotting a strand-ambiguous SNP

```{code-cell} ipython3
def ambiguous(a1, a2):
    """A/T and C/G SNPs look the same on both strands."""
    return {a1, a2} in ({"A", "T"}, {"C", "G"})

for snp, a1, a2 in [("rs7903146", "C", "T"), ("rs13266634", "C", "T"), ("rs1801282", "C", "G")]:
    print(snp, a1 + "/" + a2, "ambiguous – remove" if ambiguous(a1, a2) else "ok")
```

*PPARG* rs1801282 is a C/G SNP. Its alleles are strand-ambiguous, so a conservative pipeline would remove it unless
its strand can be confirmed — for example from the allele frequencies, or because both files are known to use the
forward strand of the same genome build.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/prs-quality-control.json")
```
