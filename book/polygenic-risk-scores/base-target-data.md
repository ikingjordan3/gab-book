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

# Base and target data

```{admonition} Module at a glance
**Lecture:** Polygenic risk scores 1 &nbsp;·&nbsp; **Prerequisite:** {doc}`prs-workflow`

**After this module you should be able to:** use the PGS Catalog and GWAS Catalog, read summary statistics with either allele convention, and calculate a PRS by hand.
```

## The PGS Catalog

The **Polygenic Score (PGS) Catalog** ([www.pgscatalog.org](https://www.pgscatalog.org/); Lambert et al.,
[2021](https://doi.org/10.1038/s41588-021-00783-5)) is an open database of published polygenic scores. For every
score it records:

- the **scoring file**: the variants, effect alleles and weights, which you can download;
- how the score was **developed**: the method, the GWAS it was built from, and the samples used to train it;
- how it **performed** when evaluated, including the ancestry of every sample.

The home page shows the current numbers of scores, traits and publications, which grow with every release.
The catalog also provides **pgsc_calc**, a reproducible workflow for applying catalog scores to your own
genotype data.

```{figure} images/pgs_home.png
:alt: PGS Catalog home page, with a search box and counts of polygenic scores, traits and publications.
:width: 100%

The PGS Catalog home page. Screenshot from [www.pgscatalog.org](https://www.pgscatalog.org/).
```

Searching for a trait lists every published score for it. The type 2 diabetes page
([MONDO_0005148](https://www.pgscatalog.org/trait/MONDO_0005148)) lists about 250 scores, built with anything from
62 variants to almost 7 million.

```{figure} images/pgs_trait_t2d.png
:alt: PGS Catalog trait page for type 2 diabetes mellitus listing polygenic scores; the first rows are PGS000014 (Khera 2018, 6,917,436 variants), PGS000020 and PGS000031.
:width: 100%

Polygenic scores for type 2 diabetes in the PGS Catalog.
```

### A published score: PGS000014

[PGS000014](https://www.pgscatalog.org/score/PGS000014) is Khera et al.'s type 2 diabetes score from
{doc}`prs-workflow`. Its entry shows how the score was constructed — LDpred with ρ = 0.01 and an LD panel of 503
Europeans from 1000 Genomes, giving weights for 6,917,436 variants — and the ancestry of the samples used at each
stage. Those stages are the three samples of {doc}`prs-workflow`:

| PGS Catalog field | PGS000014 | Data set |
|---|---|---|
| Source of variant associations (GWAS) | European, 318,416 individuals | **Base data** |
| Score development / training | European, 120,280 individuals (UK Biobank phase 1) | **Validation data** |
| PGS evaluation | 23 sample sets of many ancestries | **Testing data** |

```{figure} images/pgs000014.png
:alt: PGS Catalog entry PGS000014 (GPS_T2D) showing score construction (LDpred, rho = 0.01, 6,917,436 variants) and the ancestry distribution of the GWAS source, score development and evaluation samples.
:width: 100%

The PGS Catalog entry for PGS000014.
```

Scores in the PGS Catalog have already been validated, so they can be applied directly to new target data. The
ancestry columns matter: a score built from European-ancestry GWAS usually predicts less well in other ancestries.

## Base data: a type 2 diabetes GWAS

The base data for a new type 2 diabetes score should come from the largest, most diverse GWAS available. The
largest to date is the multi-ancestry meta-analysis by the Type 2 Diabetes Global Genomics Initiative (T2DGGI):
Suzuki et al., *Genetic drivers of heterogeneity in type 2 diabetes pathophysiology*
([Nature, 2024](https://doi.org/10.1038/s41586-024-07019-6)).

| | |
|---|---|
| Individuals | 2,535,601 |
| Type 2 diabetes cases | 428,452 |
| Not of European ancestry | 39.7% |
| Independent genome-wide significant signals | 1,289, at 611 loci |
| Loci not reported before | 145 |
| Clusters of signals with distinct mechanisms | 8 |

The eight clusters group signals by their effects on other cardiometabolic traits, for example signals acting
through beta-cell function, obesity or insulin resistance. Its summary statistics — one row per variant — are
base data for new type 2 diabetes scores.

## The GWAS Catalog

The **NHGRI-EBI GWAS Catalog** ([www.ebi.ac.uk/gwas](https://www.ebi.ac.uk/gwas/); Sollis et al.,
[2023](https://doi.org/10.1093/nar/gkac1010)) is the curated database of published GWAS, run by EMBL-EBI and the
US National Human Genome Research Institute. Curators extract every significant variant–trait association from
each study: the variant, its risk allele, P-value, OR or β, mapped gene and study accession (GCST…). The catalog
also hosts full summary statistics where authors share them. You can search by trait, variant (rsID), gene,
genomic region, study accession or PubMed ID.

```{figure} images/gwascat_home.png
:alt: GWAS Catalog home page, the NHGRI-EBI Catalog of human genome-wide association studies, with a search box and links to download, summary statistics and submission.
:width: 100%

The GWAS Catalog home page. Screenshot from [www.ebi.ac.uk/gwas](https://www.ebi.ac.uk/gwas/).
```

For rs7903146 in *TCF7L2*, the strongest common-variant type 2 diabetes signal, the catalog lists more than 70
type 2 diabetes associations. Every one reports the same risk allele, **T**, with an OR above 1 (1.30–1.44).

```{figure} images/gwascat_rs7903146.png
:alt: GWAS Catalog associations for rs7903146 filtered to type 2 diabetes: risk allele rs7903146-T, P-values from 1e-11 to 4e-94, odds ratios 1.30 to 1.44, mapped gene TCF7L2.
:width: 100%

Type 2 diabetes associations for [rs7903146](https://www.ebi.ac.uk/gwas/variants/rs7903146) in the GWAS Catalog.
```

## Two ways to report the same result

Summary statistics can report each SNP in two ways, and it is easy to mix them up. Here are three type 2 diabetes
lead SNPs from {doc}`../complex-traits-gwas/significance-and-plots`. The SNPs, positions (GRCh38) and alleles are
real; the ORs and P-values are the illustrative values used in Chapter 7.

**Effect allele = alternate allele.** Most GWAS software codes the alternate allele as the effect allele: the
allele that was tested in the model and counted as 0, 1 or 2 copies.

| Chromosome | Position | rsID (gene) | Ref allele | Alt allele (effect) | OR | P-value |
|---|---|---|---|---|---|---|
| 10 | 112998590 | rs7903146 (*TCF7L2*) | C | **T** | 1.37 | 2.1 × 10⁻³⁸ |
| 8 | 117172544 | rs13266634 (*SLC30A8*) | C | **T** | 0.89 | 2.7 × 10⁻¹⁵ |
| 3 | 12351626 | rs1801282 (*PPARG*) | C | **G** | 0.88 | 7.9 × 10⁻¹² |

At *TCF7L2* the alternate allele increases risk (OR > 1). At *SLC30A8* (T, the Arg325Trp variant) and *PPARG*
(G, the Pro12Ala variant) the alternate allele is protective (OR < 1).

**Risk allele convention.** Many databases, including the GWAS Catalog, report the risk-increasing allele instead,
so every OR is above 1:

| Chromosome | Position | rsID (gene) | Risk allele | Other allele | OR | P-value |
|---|---|---|---|---|---|---|
| 10 | 112998590 | rs7903146 (*TCF7L2*) | **T** | C | 1.37 | 2.1 × 10⁻³⁸ |
| 8 | 117172544 | rs13266634 (*SLC30A8*) | **C** | T | **1.12** | 2.7 × 10⁻¹⁵ |
| 3 | 12351626 | rs1801282 (*PPARG*) | **C** | G | **1.13** | 7.9 × 10⁻¹² |

For *SLC30A8* and *PPARG* the alleles swap and the OR inverts: 1 / 0.89 = 1.12 and 1 / 0.88 = 1.13. The P-value
does not change. This can be confusing — always check which allele an OR refers to before using it as a PRS
weight ({doc}`prs-quality-control`).

## Target data and a PRS by hand

Target data are individual-level genotypes at the same SNPs, coded as the number of alternate alleles (0, 1 or 2)
each person carries. The genotypes below are made up for the exercise.

| rsID | Ref | Alt | Individual 1 | Individual 2 | … | Individual n |
|---|---|---|---|---|---|---|
| rs7903146 | C | T | 2 | 1 | | 0 |
| rs13266634 | C | T | 0 | 1 | | 2 |
| rs1801282 | C | G | 1 | 2 | | 1 |

Because the genotypes count **alternate** alleles, the weights must be the log ORs of the **alternate** alleles
(the first table). Then each person's score is the weighted sum from {doc}`prs-defined`:

```{code-cell} ipython3
import math

alt_or = {"rs7903146": 1.37, "rs13266634": 0.89, "rs1801282": 0.88}   # OR of the alternate (effect) allele
weights = {snp: math.log(o) for snp, o in alt_or.items()}           # beta = ln(OR)

genotypes = {                       # number of alternate alleles
    "Individual 1": {"rs7903146": 2, "rs13266634": 0, "rs1801282": 1},
    "Individual 2": {"rs7903146": 1, "rs13266634": 1, "rs1801282": 2},
    "Individual n": {"rs7903146": 0, "rs13266634": 2, "rs1801282": 1},
}

for person, g in genotypes.items():
    prs = sum(g[snp] * weights[snp] for snp in weights)
    print(f"{person}: PRS = {prs:+.3f}")
```

Individual 1 has the highest score: two copies of the *TCF7L2* risk allele and the *SLC30A8* risk allele (C) on
both chromosomes. Using the risk-allele ORs with alternate-allele genotype counts would give the wrong answer — the
error that the quality-control step is designed to catch.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/base-target-data.json")
```
