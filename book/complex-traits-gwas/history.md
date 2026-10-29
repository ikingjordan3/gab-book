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

# A short history of complex-disease genetics

```{admonition} Module at a glance
**Lecture:** Genetics of complex, common disease &nbsp;·&nbsp; **Prerequisite:** {doc}`variation-spectrum`

**After this module you should be able to:** trace the move from linkage studies to candidate genes to GWAS, and explain why each transition happened.
```

The methods used to find disease genes have changed twice since the 1980s. Each change happened because
the previous method had reached its limits for common, complex disease, and new technology made something
better possible.

```{figure} images/history_timeline.png
:alt: Timeline from 1980 to the present with three eras. Linkage studies, 1980s to 1990s, including CFTR mapped by linkage in 1989 and the Huntington's disease gene HTT in 1993. Candidate genes, late 1990s to early 2000s, including the human genome draft sequence in 2001. GWAS, 2005 to present, including HapMap and the first GWAS in 2005, the WTCCC study of seven common diseases in 2007, and a 2024 type 2 diabetes GWAS of 2.5 million people that found more than 600 loci.
:width: 100%

Three eras of complex-disease genetics, with a few landmark events.
```

## Linkage studies (1980s–1990s)

**Linkage studies** map disease genes using families. A pedigree is genotyped at genetic markers spread
across the genome, and the analysis looks for a marker that is inherited together with the disease more
often than chance would allow. Markers and causal variants that are physically close on a chromosome are
rarely separated by recombination, so they travel together through the generations. (This co-inheritance
within families is **linkage**; the related population-level correlation between nearby variants is
**linkage disequilibrium**, which becomes central to GWAS.)

**Linkage was highly successful for rare Mendelian disorders.** Genes for cystic fibrosis (*CFTR*, 1989)
and Huntington's disease (*HTT*, 1993) were found this way. In these diseases a single variant of large
effect and high penetrance drives the phenotype, so it segregates cleanly through a pedigree.

**Linkage studies did not work for complex common diseases** such as diabetes and schizophrenia:

- They are **statistically underpowered** to detect small effect sizes spread across many loci. A variant
  that raises risk by 20% barely changes who is affected within a single family.
- They are also underpowered to detect **common variants**, whose patterns of linkage disequilibrium are
  shaped by population history rather than by inheritance within one family.
- Common diseases show **genetic and phenotypic heterogeneity**: different families can carry different
  risk variants, and the same diagnosis can cover several underlying conditions.

## Candidate genes (late 1990s–early 2000s)

The failure of linkage led to a shift towards **population-based association studies**: instead of
following variants through families, compare variant frequencies between unrelated people with and without
the disease (cases and controls; see {doc}`case-control-design`).

The first association studies were **hypothesis-driven**. Researchers chose a handful of **candidate genes**
that seemed biologically plausible (for example, an insulin-signalling gene for diabetes) and tested
variants in them with simple statistical methods. A typical study looked like this:

| Genes tested | Variants tested | Individuals per cohort |
|---|---|---|
| 1–10 | hundreds to low thousands | 100–500 |

**Many associations were reported, generating excitement, followed by massive failure to replicate.**
When other groups tested the same variants in new samples, most associations disappeared. The reasons are
now well understood:

- **low statistical power**: samples of a few hundred people cannot reliably detect small effects;
- **no correction for multiple testing**: testing many variants and reporting any P < 0.05 guarantees
  false positives;
- **publication bias**: positive results were published and negative ones were not;
- **ignoring confounders** such as population structure (see {doc}`covariates`); and
- **disregarding polygenic architecture**: real effects are small and spread across the genome, so a
  plausible gene is rarely where the signal is.

## Genome-wide association studies (2005–present)

Disillusionment with candidate gene studies, together with advances in genome technology, led to the
**GWAS era**. Three advances made it possible:

- **The Human Genome Project** (draft 2001, completed 2003) provided the genomic landscape and a catalogue of
  genes, the reference against which all variants are described.
- **The International HapMap Project** (2005) mapped common variants and characterized **linkage
  disequilibrium (LD) blocks**, showing that a few hundred thousand well-chosen SNPs could "tag" most common
  variation in the genome.
- **Affordable, high-throughput SNP arrays** genotyped hundreds of thousands of variants per person, at a
  cost that allowed much larger cohorts.

A **GWAS** is **hypothesis-free**: it tests every genotyped variant (hundreds of thousands to millions) in
large cohorts, without first guessing which genes matter. Because so many tests are run, it uses
**stringent statistical thresholds**, conventionally P < 5 × 10⁻⁸ ({doc}`significance-and-plots`), and
requires replication in independent samples. The approach has been very productive. One of the first
successes, in 2005, found a variant in the complement factor H gene (*CFH*) associated with age-related
macular degeneration. Today GWAS of millions of people find hundreds of loci for a single disease, such as
the more than 600 loci for type 2 diabetes.

| | Linkage studies | Candidate genes | GWAS |
|---|---|---|---|
| **Samples** | Families (pedigrees) | Unrelated cases and controls | Unrelated cases and controls |
| **Scope** | Genome-wide markers | 1–10 genes chosen in advance | Whole genome, hypothesis-free |
| **Best at finding** | Rare variants of large effect | — (poor replication) | Common variants of small effect |
| **Significance** | LOD score in families | P < 0.05, often uncorrected | P < 5 × 10⁻⁸, plus replication |

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/history.json")
```
