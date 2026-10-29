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

### GREML: genetic similarity measured from genotypes

Two people picked at random from a population are not relatives, but they still share slightly more or slightly
fewer alleles than average, purely by chance. GREML (**g**enomic-relatedness-based **r**estricted **m**aximum
**l**ikelihood) uses those tiny differences in genetic similarity among unrelated people the way a twin study uses
MZ and DZ twins.

**Step 1: the genetic relationship matrix (GRM).** For every pair of people $j$ and $k$, genotyped at $M$ SNPs,

$$
A_{jk} = \frac{1}{M} \sum_{i=1}^{M} \frac{(x_{ij} - 2p_i)(x_{ik} - 2p_i)}{2p_i(1-p_i)}
$$

where $x_{ij}$ is the number of copies (0, 1 or 2) of the reference allele that person $j$ carries at SNP $i$ and
$p_i$ is that allele's frequency. $A_{jk}$ is a genome-wide correlation between the two people's genotypes: about
0.5 for full siblings, and close to 0 for unrelated people, with small positive and negative values from chance
sharing. Close relatives are removed (typically pairs with $A_{jk} > 0.025$) so that shared family environment
cannot masquerade as genetic similarity.

**Step 2: a model for the phenotypes.** Each person's trait value is the sum of a genetic effect $g$, from all the
SNPs together, and an environmental effect $e$. People who are genetically more similar have more similar genetic
effects, so

$$
\operatorname{Var}(\mathbf{y}) = \mathbf{A}\,\sigma^2_g + \mathbf{I}\,\sigma^2_e,
\qquad
h^2_\text{SNP} = \frac{\sigma^2_g}{\sigma^2_g + \sigma^2_e}
$$

where $\mathbf{A}$ is the GRM and $\mathbf{I}$ the identity matrix. GREML finds the values of $\sigma^2_g$ and
$\sigma^2_e$ that best fit the observed phenotypes, by restricted maximum likelihood, as implemented in the
software GCTA (Yang et al., [2011](https://doi.org/10.1016/j.ajhg.2010.11.011)).

**The intuition.** For a standardized trait, the expected product of two people's trait values is
$A_{jk} \, h^2_\text{SNP}$. So if you plot the phenotypic similarity of every pair of people against their genetic
similarity, the slope of the line is the SNP-heritability. (This simpler regression is called Haseman–Elston
regression; GREML uses the same information more efficiently.) The simulation below does exactly that for 2,000
unrelated people genotyped at 3,000 SNPs, with a true SNP-heritability of 0.5:

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

### LD score regression: heritability from summary statistics

GREML needs every individual's genotypes, which large consortia often cannot share. **LD score regression**
(Bulik-Sullivan et al., [2015](https://doi.org/10.1038/ng.3211)) estimates SNP-heritability from GWAS summary
statistics alone.

**The LD score.** A SNP is associated with a trait if it is a causal variant or if it is in linkage
disequilibrium (LD) with one. The **LD score** of SNP $j$ measures how much other variation it tags:

$$
\ell_j = \sum_k r^2_{jk}
$$

the sum of its squared correlations ($r^2$) with all nearby SNPs $k$, including itself. LD scores are computed
once from a reference panel with the same ancestry as the GWAS.

**The regression.** If a trait is polygenic, causal variants are scattered throughout the genome, so a SNP that
tags more variation (a higher LD score) is expected to show a stronger association. For a GWAS of $N$ people
and $M$ SNPs, the expected association statistic of SNP $j$ is

$$
\operatorname{E}\!\left[\chi^2_j\right] = \frac{N\,h^2_\text{SNP}}{M}\,\ell_j \;+\; N a \;+\; 1
$$

Regressing each SNP's $\chi^2$ on its LD score gives two useful numbers:

- the **slope**, which estimates $h^2_\text{SNP}$; and
- the **intercept**, $1 + Na$, which measures confounding. Population stratification or cryptic relatedness
  ($a$) inflates the $\chi^2$ of every SNP equally, whatever its LD score, so it raises the intercept above 1
  but leaves the slope alone.

The simulation below generates $\chi^2$ statistics for a million SNPs from a GWAS of 50,000 people, with and
without confounding:

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
regression separates the two, a point that returns in {doc}`significance-and-plots`.

### Missing heritability

For most traits, SNP-heritability is lower than twin heritability, and early GWAS hits
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
