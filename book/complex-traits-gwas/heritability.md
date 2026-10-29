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

# Heritability

```{admonition} Module at a glance
**Lecture:** Genetics of complex, common disease &nbsp;·&nbsp; **Prerequisite:** {doc}`polygenic-model`

**After this module you should be able to:** define heritability, partition genetic variance into additive, dominance and epistatic parts, distinguish broad- from narrow-sense heritability, estimate h² from twins with the ACE model, and correct common misreadings of heritability.
```

Fisher's model ({doc}`polygenic-model`) says that a complex trait is shaped by many genes plus the
environment. The natural next question is *how much* of the variation we see is due to genes. Heritability
answers that question, and it is one of the most useful and most misunderstood numbers in genetics.

## Definition

**Heritability** is the proportion of the phenotypic variation in a population that is attributable to
genetic differences among individuals. It measures how relevant genetic variation is to phenotypic
variation, which is loosely the genetic influence on a trait.

$$
H^2 = \frac{V_G}{V_P}, \qquad \text{where } V_P = V_G + V_E
$$

- $V_P$ is the **phenotypic variance**: how much the trait varies across the population.
- $V_G$ is the **genetic variance**: the part of that variation due to differences in genotype.
- $V_E$ is the **environmental variance**: everything else, including measurement error.

Heritability runs from 0 (none of the variation is genetic) to 1 (all of it is). Note that it is defined
in terms of **variance**, a measure of differences between people, and not in terms of the trait value of
any one person.

## Partitioning genetic variance

Fisher also showed that genetic variance can itself be split into three parts:

$$
V_G = V_A + V_D + V_I
$$

```{figure} images/vc_panels.png
:alt: Three panels. Additive: trait value rises by the same amount from aa to Aa to AA. Dominance: Aa has the same value as AA, above the dashed additive expectation by an amount d. Epistasis: the trait rises across bb, Bb, BB when locus A is AA, but stays flat when locus A is aa.
:width: 100%

The three components of genetic variance, illustrated with genotypic values (arbitrary units, not data).
```

- **Additive variance ($V_A$).** Each allele adds a fixed amount to the trait, and effects **sum** across
  alleles and loci. In the left panel, each copy of A raises the trait by the same amount, *a*.
- **Dominance variance ($V_D$).** This comes from interaction between the two alleles **at one locus**: the
  heterozygote is not midway between the homozygotes. In the middle panel, Aa has the same value as AA,
  a deviation *d* above the additive expectation.
- **Epistatic (interaction) variance ($V_I$).** This comes from interaction between alleles **at different
  loci**: the effect of one locus depends on the genotype at another. In the right panel, locus B raises
  the trait only in people who are AA at locus A.

Only the additive part is reliably **passed from parent to offspring**. A parent transmits single alleles,
not genotypes, and dominance and epistatic effects depend on particular combinations of alleles that are
broken up each generation by segregation and recombination.

## Broad- vs. narrow-sense heritability

That distinction gives two kinds of heritability:

| | Broad-sense heritability | Narrow-sense heritability |
|---|---|---|
| **Definition** | $H^2 = \dfrac{V_G}{V_P} = \dfrac{V_A + V_D + V_I}{V_P}$ | $h^2 = \dfrac{V_A}{V_P}$ |
| **Counts** | All genetic variance | Additive variance only |
| **Tells you** | How well genotype predicts phenotype in this population | How much relatives resemble each other, and how a population responds to selection |

Narrow-sense heritability is the one used most in human genetics. It sets the slope of the regression of
offspring on the average of their parents, it is what twin and family studies usually estimate, and it is
what GWAS and polygenic risk scores try to capture, because they model additive effects of individual
alleles. In general

$$
h^2_\text{SNP} \;\le\; h^2 \;\le\; H^2
$$

where $h^2_\text{SNP}$ is the SNP-heritability estimated from GWAS data ({doc}`gwas-to-prs`).

## Estimating h² from relatives

Every method for estimating heritability works the same way: it asks whether people who are more similar
genetically are also more similar in the trait. The methods differ in how they measure genetic similarity.

One way is from the **resemblance between relatives**, because relatives share known proportions of their
genetic variance. Another is from genotype data: SNP arrays measure genetic similarity directly, even between
people who are not related, which gives the SNP-heritability described in {doc}`gwas-to-prs`.
The table shows what relatives share:

| Relatives | Genetic variance shared |
|---|---|
| Identical (monozygotic, MZ) twins | $V_A + V_D$ (all genes) |
| Fraternal (dizygotic, DZ) twins and full siblings | $\tfrac{1}{2}V_A + \tfrac{1}{4}V_D$ |
| Parent and offspring | $\tfrac{1}{2}V_A$ |
| Half siblings | $\tfrac{1}{4}V_A$ |

(Epistatic variance is shared too, in smaller fractions, and is left out here.) Twins are especially useful
because MZ and DZ twins are the same age and are usually raised together, yet differ in how many genes
they share.

### The ACE twin model

The classical twin design uses the **ACE model**, which splits phenotypic variance into three sources, each
expressed as a proportion of the total:

- **A**: additive genetic variance, $a^2$ (this is $h^2$);
- **C**: common, or shared, environment, $c^2$: everything that makes twins raised together alike, such as
  family, home and school;
- **E**: unique environment, $e^2$: everything that makes twins differ, including measurement error.

$$
a^2 + c^2 + e^2 = 1
$$

**What twins share.** The correlation between twins equals the proportion of variance that comes from the
sources they share. MZ twins share all of their genes and their family environment. DZ twins share, on
average, half of the additive genetic variance and all of their family environment:

$$
r_{MZ} = a^2 + c^2, \qquad r_{DZ} = \tfrac{1}{2}a^2 + c^2
$$

**Solve.** Subtracting the second equation from the first gives $r_{MZ} - r_{DZ} = \tfrac{1}{2}a^2$, and
substituting back gives the other two components:

$$
a^2 = 2\,(r_{MZ} - r_{DZ}), \qquad c^2 = 2\,r_{DZ} - r_{MZ}, \qquad e^2 = 1 - r_{MZ}
$$

The first of these, $h^2 \approx 2(r_{MZ} - r_{DZ})$, is known as **Falconer's formula**. Try it with
different twin correlations:

```{code-cell} ipython3
def ace(r_mz, r_dz):
    a2 = 2 * (r_mz - r_dz)      # additive genetic (h2)
    c2 = 2 * r_dz - r_mz        # shared environment
    e2 = 1 - r_mz               # unique environment
    return a2, c2, e2

for r_mz, r_dz in [(0.80, 0.50), (0.80, 0.40), (0.50, 0.40)]:
    a2, c2, e2 = ace(r_mz, r_dz)
    print(f"rMZ={r_mz}, rDZ={r_dz}: a2={a2:.2f}, c2={c2:.2f}, e2={e2:.2f}")
```

With $r_{MZ} = 0.80$ and $r_{DZ} = 0.50$, 60% of the variance is additive genetic, 20% is shared
environment and 20% is unique environment. When $r_{DZ}$ is exactly half of $r_{MZ}$ there is no shared
environment at all, and when MZ twins are only a little more alike than DZ twins most of the familial
resemblance is environmental.

**Assumptions.** The ACE model assumes that MZ and DZ twins share their environments to the same degree
(the *equal environments assumption*), that genetic effects are additive, that mating is random and that
twins are representative of the population. If $r_{MZ}$ is more than twice $r_{DZ}$, dominance is suspected
and an ADE model is fitted instead. In practice these models are fitted by maximum likelihood rather than
with hand formulas, but the logic is the same.

## A big misconception: heritability is not determination

```{admonition} Misconception
:class: error
High heritability implies genetic determination.
```

```{admonition} Reality
:class: tip
- A high heritability means that *most of the variation that is observed in the present population is
  caused by variation in genotypes*.
- It means that, in the current population, the genotype of an individual is a good predictor of the
  phenotype.
- However, it **does not mean** that the phenotype is determined once we know the genotype, because the
  environment can change or can be manipulated to alter the phenotype.
```

Height is a classic example. Its heritability is about 0.8 in many populations, yet average height has
increased substantially over the past century as nutrition and health improved. The genes did not change;
the environment did. A high heritability describes variation *within* a population under its current
conditions; it does not limit what a change in environment can do.

## Heritability is a population-level statistic

Heritability **cannot be applied to individuals**. If the heritability of trait X is 0.5, it **does not**
mean that half the reason someone exhibits trait X is genetic, while the other half is environmental. It
means that, in this population, half of the *differences between people* in trait X are associated with
genetic differences.

For one person, genes and environment are not separable shares. Both are needed for every trait, and
their effects interact. Heritability is a property of a trait **in a particular population at a particular
time**, and it can differ between populations and change when the environment changes.

These two misconceptions are discussed in detail by Visscher, Hill and Wray in *Heritability in the genomics
era — concepts and misconceptions* ([Nature Reviews Genetics, 2008](https://doi.org/10.1038/nrg2322)).

## How heritable are human traits?

One of the largest answers to that question comes from a meta-analysis of fifty years of twin studies by Polderman
and colleagues ([Nature Genetics, 2015](https://doi.org/10.1038/ng.3285)). They combined twin correlations and
reported variance components for:

| Traits | Publications | Twin pairs |
|---|---|---|
| 17,804 | 2,748 | 14,558,903 (partly dependent) |

covering virtually all published twin studies of complex traits. Their two headline findings were:

- "Our results provide compelling evidence that **all human traits are heritable**: not one trait had a
  weighted heritability estimate of zero."
- "… and across all traits the reported heritability is **49%**."

For a majority of traits (69%), the twin correlations were consistent with a simple model in which twin
resemblance is due solely to additive genetic variation. So genes matter for essentially every trait, but
on average they account for about half of the variation, leaving plenty of room for the environment.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/heritability.json")
```
