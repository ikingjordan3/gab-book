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

# SNP-heritability from GWAS

```{admonition} Module at a glance
**Lecture:** Genome-wide association studies &nbsp;·&nbsp; **Prerequisite:** {doc}`significance-and-plots`

**After this module you should be able to:** explain how SNP-heritability is estimated from GWAS data with GREML and with LD score regression, compare the strengths and weaknesses of the two methods, and explain why SNP-heritability is lower than twin heritability ("missing heritability").
```

**Heritability** is the proportion of the phenotypic variation in a population that is attributable to genetic
differences among individuals ({doc}`heritability`). Methods for estimating heritability work the same way: they
ask whether people who are more similar genetically are also more similar in the trait.

GWAS data can be used to estimate heritability: because every participant is genotyped ({doc}`case-control-design`),
genetic similarity can be **measured** from SNPs rather than inferred from family relationships. The result is
the **SNP-heritability**, $h^2_\text{SNP}$: the proportion of trait variance explained by the common variants
included in the study ({doc}`heritability-to-gwas`).

```{figure} images/variance_explained.png
:alt: Four nested bars. The top bar is all phenotypic variance. Below it, the variance splits into a genetic part (heritability, from twins) and an environmental part. Below the genetic part, common SNPs capture SNP-heritability, and the rest is not tagged by SNP arrays. Below that, only part of the SNP-heritability comes from genome-wide significant GWAS loci; the rest comes from SNPs whose effects are below the significance threshold.
:width: 100%

SNP-heritability is the part of heritability captured by common SNPs. Bar widths are illustrative.
```

This module covers the two most widely used methods, which need different kinds of GWAS data:

| | GREML | LD score regression |
|---|---|---|
| **Input** | Individual-level genotypes and phenotypes | GWAS summary statistics + reference LD scores |
| **Key idea** | Genetically more similar people are more similar in the trait | SNPs that tag more variation show stronger associations |
| **Also gives** | Genetic correlations; heritability by chromosome or annotation | Confounding (the intercept); genetic correlations; heritability by annotation |

## GREML: genetic similarity measured from genotypes

### Overview

GREML (**g**enomic-relatedness-based **r**estricted **m**aximum **l**ikelihood) estimates SNP-heritability from
the individual-level genotypes of a large sample of **unrelated** people. It fits a linear mixed model in which
the genetic similarity between every pair of people, measured from their SNPs, is used to explain how similar
their trait values are. It was introduced to explain the "missing heritability" of height (Yang et al.,
[2010](https://doi.org/10.1038/ng.608)) and is implemented in the software GCTA (Yang et al.,
[2011](https://doi.org/10.1016/j.ajhg.2010.11.011)).

### The intuition

Two people picked at random from a population are not relatives, but they still share slightly more or slightly
fewer alleles than average, purely by chance. If a trait is heritable, pairs of people who happen to be a little
more similar genetically should also be a little more similar in the trait. GREML uses these tiny differences in
genetic similarity among unrelated people the way a twin study uses the difference between MZ and DZ twins
({doc}`heritability`). Because the people are unrelated, they do not share a family environment, so their
similarity cannot be explained by growing up together.

At its core this is a regression across **pairs of people**. For each pair, the genetic similarity (their entry in
the genetic relationship matrix defined below) goes on the x-axis, and their trait similarity (the product of their
standardized trait values) goes on the y-axis. If the trait is heritable, the points rise from left to right, and
the slope of the line estimates $h^2_\text{SNP}$.

```{figure} images/greml_intuition.png
:alt: Schematic scatter plot. The x-axis is genetic similarity of a pair of unrelated people, from minus 0.02 to plus 0.02; the y-axis is their trait similarity. Points, each an average over many pairs, rise along a straight line through the origin; a slope triangle marks the slope as SNP-heritability.
:width: 520px

The regression behind GREML (simulated). Pairs of unrelated people who are slightly more alike genetically are
slightly more alike in the trait; the slope is $h^2_\text{SNP}$. Genetic similarity among unrelated people varies
only by a few hundredths, so very many pairs are needed.
```

### The mathematical model

**Genetic similarity.** For every pair of people $j$ and $k$, genotyped at $M$ SNPs, the **genetic relationship
matrix (GRM)** is

$$
A_{jk} = \frac{1}{M} \sum_{i=1}^{M} \frac{(x_{ij} - 2p_i)(x_{ik} - 2p_i)}{2p_i(1-p_i)}
$$

where $x_{ij}$ is the number of copies (0, 1 or 2) of the reference allele that person $j$ carries at SNP $i$ and
$p_i$ is that allele's frequency. $A_{jk}$ is a genome-wide correlation between the two people's genotypes: about
0.5 for full siblings, and close to 0 for unrelated people, with small positive and negative values from chance
sharing.

**The trait.** Each person's trait value is the sum of fixed covariate effects, a random genetic effect $g$ from
all the SNPs together, and a random environmental effect $e$:

$$
\mathbf{y} = \mathbf{X}\boldsymbol{\beta} + \mathbf{g} + \mathbf{e},
\qquad
\mathbf{g} \sim N(\mathbf{0},\, \mathbf{A}\,\sigma^2_g),
\qquad
\mathbf{e} \sim N(\mathbf{0},\, \mathbf{I}\,\sigma^2_e)
$$

where $\mathbf{X}$ holds the covariates (age, sex, ancestry PCs) and $\mathbf{I}$ is the identity matrix. People
who are genetically more similar have more similar genetic effects, so the covariance between trait values is

$$
\operatorname{Var}(\mathbf{y}) = \mathbf{A}\,\sigma^2_g + \mathbf{I}\,\sigma^2_e,
\qquad
h^2_\text{SNP} = \frac{\sigma^2_g}{\sigma^2_g + \sigma^2_e}
$$

For a standardized trait, this means the expected product of two people's trait values is
$A_{jk}\,h^2_\text{SNP}$: trait similarity rises in proportion to genetic similarity, with slope $h^2_\text{SNP}$.

### How it works, step by step

1. **Quality-control the genotypes.** Keep well-genotyped (or well-imputed) common SNPs and well-genotyped people.
2. **Compute the GRM** for all pairs of people from the genome-wide SNPs.
3. **Remove close relatives**, typically one person from each pair with $A_{jk} > 0.025$, so that shared family
   environment cannot masquerade as genetic similarity.
4. **Fit the mixed model** with covariates (age, sex and the top principal components, {doc}`covariates`) as fixed
   effects. Restricted maximum likelihood (REML) finds the values of $\sigma^2_g$ and $\sigma^2_e$ that best fit
   the observed trait values.
5. **Report** $h^2_\text{SNP} = \sigma^2_g / (\sigma^2_g + \sigma^2_e)$ with its standard error. For a disease,
   the estimate on the observed 0/1 scale is converted to the liability scale ({doc}`polygenic-model`) using the
   disease prevalence.

### Advantages and disadvantages

**Advantages**

- Uses **unrelated people**, so estimates are not inflated by shared family environment, a common criticism of
  twin studies.
- Gives a **standard error**, and the model can be extended: partitioning $h^2_\text{SNP}$ by chromosome or by
  functional category, or estimating the **genetic correlation** between two traits.
- Uses all of the information in the individual-level data, so for a given sample size it is more precise than
  methods that use summary statistics only.

**Disadvantages**

- Needs **individual-level genotypes**, which often cannot be shared between studies because of privacy and
  consent restrictions.
- **Computationally demanding**: the GRM has one entry for every pair of people, so memory grows with the square
  of the sample size, which limits standard GREML to tens of thousands of people (faster approximations exist).
- Needs **large samples**: because genetic similarities between unrelated people are tiny, thousands of people
  are needed for a useful estimate. The standard error shrinks roughly in proportion to $1/N$; with about 10,000
  unrelated people it is around 0.03.
- Captures only variation **tagged by the genotyped SNPs**, and rests on assumptions about how effect sizes
  relate to allele frequency and LD; if causal variants are rarer than the SNPs used, $h^2_\text{SNP}$ is
  underestimated.
- Residual **population stratification** can inflate the estimate if ancestry is not fully accounted for.

### Illustration: a simulation

The simulation below follows the intuition directly. It plots the phenotypic similarity of every pair of 2,000
unrelated people against their genetic similarity and fits a line through the origin. (This simple regression is
called Haseman–Elston regression; GREML uses the same information more efficiently.) The true SNP-heritability is
0.5:

```{code-cell} ipython3
import numpy as np
import matplotlib.pyplot as plt

rng = np.random.default_rng(7)
N, M, h2 = 2_000, 3_000, 0.5
p = rng.uniform(0.05, 0.5, M)            # allele frequencies
X = rng.binomial(2, p, size=(N, M))      # genotypes 0/1/2
Z = (X - 2 * p) / np.sqrt(2 * p * (1 - p))   # standardized
g = Z @ rng.normal(0, np.sqrt(h2 / M), M)   # genetic effects
y = g + rng.normal(0, np.sqrt(1 - h2), N)    # plus environment
y = (y - y.mean()) / y.std()             # standardized trait

A = Z @ Z.T / M                          # the GRM
pairs = np.triu_indices(N, k=1)          # each pair once
a = A[pairs]                             # genetic similarity
yy = np.outer(y, y)[pairs]               # trait similarity
slope = np.sum(a * yy) / np.sum(a * a)   # slope through 0
print(f"{len(a):,} pairs; estimated h2_SNP = {slope:.2f} (true value {h2})")

edges = np.quantile(a, np.linspace(0, 1, 21))   # 20 bins
which = np.digitize(a, edges[1:-1])
mid = np.array([a[which == b].mean() for b in range(20)])
mean_yy = [yy[which == b].mean() for b in range(20)]
fig, ax = plt.subplots(figsize=(6.5, 3.6))
ax.scatter(mid, mean_yy, color="#3B5B92", zorder=2, label="pairs (binned means)")
ax.plot(mid, slope * mid, color="#6B3FA0", lw=2.5, label=f"slope = {slope:.2f}")
ax.axhline(0, color="#BDBDBD", lw=1); ax.axvline(0, color="#BDBDBD", lw=1)
ax.set_xlabel("Genetic similarity $A_{jk}$")
ax.set_ylabel("Trait similarity $y_j y_k$")
ax.spines[["top", "right"]].set_visible(False); ax.legend(frameon=False)
plt.show()
```

Even though every pair is unrelated, so the genetic similarities are all close to zero, nearly two million pairs
give enough information to recover the true value.

## LD score regression: heritability from summary statistics

### Overview

**LD score regression (LDSC)** (Bulik-Sullivan et al., [2015](https://doi.org/10.1038/ng.3211)) estimates
SNP-heritability using only GWAS **summary statistics**, the per-SNP association results ({doc}`association-models`),
together with LD information from a reference panel. It also separates true polygenic signal from **confounding**
such as population stratification, which the genomic inflation factor λ ({doc}`significance-and-plots`) cannot do.

### The intuition

A SNP shows an association either because it is a causal variant or because it is in linkage disequilibrium (LD)
with one. If a trait is **polygenic**, causal variants are scattered throughout the genome. A SNP in a region of
high LD is correlated with many other variants, so it is more likely to tag a causal one, and on average it will
show a **stronger association** than a SNP in a region of low LD.

**Confounding behaves differently.** Population stratification or cryptic relatedness inflates the test statistic
of every SNP by about the same amount, whether the SNP is in high or low LD. So when association strength is
plotted against LD:

- polygenic signal makes the line **slope upward**; and
- confounding **lifts the whole line**, raising the intercept.

In this regression each point is a **SNP**: its LD score (how much variation it tags, defined below) goes on the
x-axis, and its association test statistic, χ², from the GWAS ({doc}`association-models`) goes on the y-axis.

```{figure} images/ldsc_intuition.png
:alt: Schematic scatter plot. The x-axis is the LD score of a SNP, from 0 to 250; the y-axis is mean chi-square. Points, each an average over many SNPs, rise along a solid line from an intercept of 1; the slope is labelled N times SNP-heritability over M. A dashed red parallel line starts higher, at 1.6, showing that confounding raises the intercept.
:width: 520px

The regression behind LD score regression (simulated). SNPs that tag more of the genome (higher LD score) have
larger χ² on average; the slope gives $h^2_\text{SNP}$. Without confounding the line starts at χ² = 1, the
value expected for a SNP with no effect; confounding shifts the whole line up (dashed).
```

### The mathematical model

**The LD score** of SNP $j$ measures how much other variation it tags:

$$
\ell_j = \sum_k r^2_{jk}
$$

the sum of its squared correlations ($r^2$) with all nearby SNPs $k$, including itself ($r^2_{jj} = 1$).

**The regression.** For a GWAS of $N$ people and $M$ SNPs, the expected association statistic of SNP $j$ is

$$
\operatorname{E}\!\left[\chi^2_j\right] = \frac{N\,h^2_\text{SNP}}{M}\,\ell_j \;+\; N a \;+\; 1
$$

- The **slope**, $N h^2_\text{SNP} / M$, gives the SNP-heritability.
- The **intercept**, $1 + Na$, measures confounding: $a$ is the inflation from population stratification or
  cryptic relatedness. Without confounding, the intercept is 1, the expected value of a null $\chi^2$ statistic.

### How it works, step by step

1. **Obtain GWAS summary statistics** for the trait: for each SNP, its alleles, its $z$- or $\chi^2$-statistic
   ($\chi^2 = z^2$) and the sample size. Restrict to common, well-imputed SNPs.
2. **Obtain LD scores** for the same SNPs, computed once from a reference panel whose ancestry matches the GWAS,
   such as the European samples of the 1000 Genomes Project for a European-ancestry GWAS.
3. **Regress** each SNP's $\chi^2$ on its LD score. (The real method weights SNPs to account for the fact that
   neighbouring SNPs are correlated and that $\chi^2$ is noisier at high LD scores.)
4. **Read off** $h^2_\text{SNP}$ from the slope and confounding from the intercept. Standard errors are estimated
   by a block jackknife: dropping blocks of the genome one at a time and refitting.
5. **Extend** the same idea, if needed: regressing the product of two traits' $z$-statistics on LD score
   estimates their **genetic correlation** (cross-trait LDSC), and splitting LD scores by functional annotation
   estimates which kinds of DNA contribute most heritability (stratified LDSC).

### Advantages and disadvantages

**Advantages**

- Needs **summary statistics only**, which are routinely shared, so it can be applied to the largest GWAS
  meta-analyses without access to anyone's genotypes.
- **Fast**: a single regression across the genome, run in minutes.
- **Distinguishes polygenicity from confounding**, which λ cannot do. This tells researchers whether an
  inflated Q-Q plot reflects a problem with the study or a highly polygenic trait.
- **Extends naturally** to genetic correlations between traits and to partitioning heritability by functional
  category.

**Disadvantages**

- **Less precise** than GREML for the same sample size, so it works best with large GWAS (tens of thousands of
  people or more).
- Needs a **reference LD panel that matches the GWAS ancestry**. This is difficult for admixed or multi-ancestry
  GWAS, and LD reference data are less complete for under-represented populations.
- Relies on **model assumptions**, for example that effect sizes are unrelated to LD; if these are violated, the
  estimates can be biased.
- Like GREML, it captures only variation **tagged by common SNPs**.

### Illustration: a simulation

The simulation below generates $\chi^2$ statistics for a million SNPs from a GWAS of 50,000 people, with and
without confounding, and fits the LD score regression:

```{code-cell} ipython3
rng = np.random.default_rng(3)
M, N, h2 = 1_000_000, 50_000, 0.3
ld = rng.gamma(2.0, 30.0, M) + 1          # LD scores from a reference panel
edges = np.quantile(ld, np.linspace(0, 1, 21))
which = np.digitize(ld, edges[1:-1])      # 20 bins of LD score

fig, ax = plt.subplots(figsize=(6.5, 3.6))
for Na, colour in [(0.0, "#3B5B92"), (0.15, "#B03A2E")]:
    label = "with confounding" if Na else "no confounding"
    chi2 = (1 + N * h2 * ld / M + Na) * rng.chisquare(1, M)
    slope, intercept = np.polyfit(ld, chi2, 1)
    lam = np.median(chi2) / 0.455         # genomic inflation factor
    print(f"{label}: h2_SNP = {slope * M / N:.2f}, "
          f"intercept = {intercept:.2f}, lambda = {lam:.2f}")
    x = [ld[which == b].mean() for b in range(20)]
    y = [chi2[which == b].mean() for b in range(20)]
    ax.scatter(x, y, color=colour, s=25, zorder=2)
    xs = np.array([0, max(x)])
    ax.plot(xs, intercept + slope * xs, color=colour, lw=2.5, label=label)
ax.axhline(1, color="#BDBDBD", lw=1, ls="--")
ax.set_xlabel("LD score"); ax.set_ylabel(r"Mean $\chi^2$")
ax.spines[["top", "right"]].set_visible(False); ax.legend(frameon=False)
plt.show()
```

Both analyses recover the true SNP-heritability of 0.3 from the slope. Without confounding the intercept is 1;
with confounding it is about 1.15. Notice that in both cases the average $\chi^2$ is well above 1: in a large GWAS
of a polygenic trait, many true small effects inflate the test statistics even without any confounding. LD score
regression separates the two, the question raised by the Q-Q plots in {doc}`significance-and-plots`.

## Missing heritability

For most traits, SNP-heritability is lower than twin heritability, and early GWAS hits
explained far less still. For height, twin studies give $h^2 \approx 0.8$; the first GWAS-significant SNPs
explained only about 5% of the variance, while all common SNPs together explained about 45% (Yang et al.,
[2010](https://doi.org/10.1038/ng.608)). The gap between heritability and what GWAS explain became known as
**missing heritability**. Part of it is closed by larger GWAS, which detect more small effects; part reflects
rare variants that arrays do not tag, and part may reflect overestimation of twin heritability. So, in general,

$$
h^2_\text{SNP} \;\le\; h^2 \;\le\; H^2
$$

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/snp-heritability.json")
```
