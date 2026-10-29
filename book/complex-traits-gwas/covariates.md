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

# Covariates and confounding

```{admonition} Module at a glance
**Lecture:** Genome-wide association studies &nbsp;·&nbsp; **Prerequisite:** {doc}`odds-ratios`

**After this module you should be able to:** explain why age, sex and population structure must be included as covariates in a GWAS, show how population structure can create a false association, and describe how genetic principal components correct for it.
```

The χ² test in {doc}`odds-ratios` compares allele counts in cases and controls and nothing else. **The problem
with using χ² tests is that they do not let you add covariates**: other variables that affect the trait and
must be accounted for. This module explains why covariates are needed. The three that appear in almost every
GWAS are **age**, **sex** and **population structure**.

## Age and sex

Many diseases become more common with age, and many differ between the sexes. Type 2 diabetes does both.

```{figure} images/diabetes_by_age.png
:alt: Bar chart of diabetes prevalence by age in US adults. 3.6% at ages 20 to 39, 17.7% at 40 to 59, and 27.3% at 60 and older.
:width: 600px

Total diabetes (diagnosed and undiagnosed) by age, US adults aged 20 and over, NHANES 2021–2023.
Source: CDC/NCHS Data Brief No. 516 (2024).
```

```{figure} images/diabetes_by_sex.png
:alt: Bar chart of diabetes prevalence by sex in US adults. Men 18.0%, women 13.7%.
:width: 600px

Total diabetes by sex, same survey. Overall prevalence was 15.8%.
```

Diabetes prevalence rises steeply with age, from 3.6% in adults aged 20–39 to 27.3% in those aged 60 and over,
and it is higher in men (18.0%) than in women (13.7%). Type 2 diabetes accounts for about 90–95% of these cases.

Why does this matter for a GWAS? If cases are on average older than controls, or include more men, then any
genetic variant whose frequency happens to differ by age or sex in the sample can appear associated with
diabetes. Even when no variant is affected, age and sex explain part of the variation in disease risk, so
including them in the model makes the genetic effects easier to detect. Either way, **age and sex are included
as covariates**.

## Population structure

Population structure is the most important confounder in GWAS, and the least obvious.

**Wherever there is a difference in allele frequency as well as in the prevalence of a disease or trait between
two groups, false associations can arise.** Those differences must be controlled for.

Consider a population made up of two subpopulations, *blue* and *red*, that have not been separated in the
analysis. A SNP has alleles A and G.

| Blue subpopulation (9.1% diabetic, f(A) = 0.2) | AA | AG | GG |
|---|---|---|---|
| Case | 80 | 640 | 1,280 |
| Control | 800 | 6,400 | 12,800 |
| Case/control | 0.1 | 0.1 | 0.1 |

| Red subpopulation (16.7% diabetic, f(A) = 0.5) | AA | AG | GG |
|---|---|---|---|
| Case | 200 | 400 | 200 |
| Control | 1,000 | 2,000 | 1,000 |
| Case/control | 0.2 | 0.2 | 0.2 |

Within each subpopulation the SNP has **no effect**: every genotype has the same case/control ratio, cases and
controls have the same allele frequency, and the allelic odds ratio is exactly 1.0. But the red subpopulation
has both more A alleles **and** more diabetes. Pool the two, as a study that ignores ancestry would:

| Total population (blue + red) | AA | AG | GG |
|---|---|---|---|
| Case | 280 | 1,040 | 1,480 |
| Control | 1,800 | 8,400 | 13,800 |
| Case/control | 0.155 | 0.124 | 0.107 |

Now A looks like a risk allele:

```{code-cell} ipython3
import numpy as np
from scipy.stats import chi2_contingency

# genotype counts [AA, AG, GG] for cases and controls
blue = {"case": np.array([80, 640, 1280]), "control": np.array([800, 6400, 12800])}
red  = {"case": np.array([200, 400, 200]), "control": np.array([1000, 2000, 1000])}
pooled = {k: blue[k] + red[k] for k in blue}

def allele_counts(g):              # [A, G] from [AA, AG, GG]
    return [2 * g[0] + g[1], g[1] + 2 * g[2]]

for name, grp in [("blue", blue), ("red", red), ("pooled", pooled)]:
    table = [allele_counts(grp["case"]), allele_counts(grp["control"])]
    (a, g), (ac, gc) = table
    chi2, p, _, _ = chi2_contingency(table, correction=False)
    print(f"{name:6s}: OR = {(a / g) / (ac / gc):.2f}, chi2 = {chi2:4.1f}, P = {p:.1e}")
```

In the pooled data the odds ratio is 1.2 and P ≈ 6 × 10⁻⁹, which passes the genome-wide significance threshold
of 5 × 10⁻⁸ ({doc}`significance-and-plots`). Yet the association is **entirely spurious**: the SNP has no
effect in either group. It is an artifact of population structure. (Example adapted from Gibson G., *A Primer of
Human Genetics*, Figure 6.5.)

### The chopsticks gene

A well-known parable makes the same point (Hamer & Sirota, [2000](https://doi.org/10.1038/sj.mp.4000662), after
Lander & Schork, 1994):

1. **The study.** An ethnogeneticist searches for variants linked to **manual dexterity**: university
   students are asked how often they use chopsticks, and then they are genotyped.
2. **The discovery.** A marker explains **nearly half** of the variance in chopstick use! The
   *SUSHI* gene (*s*uccessful-*u*se-of-*s*elected-*h*and-*i*nstruments) is announced.
3. **The catch.** The marker is an **HLA (histocompatibility) gene**. Its allele is more common in students of
   East Asian ancestry, who also use chopsticks more, for cultural reasons.

```{figure} images/confounding_diagram.png
:alt: Diagram. Ancestry points to both HLA allele frequency and chopstick use; a dashed line between HLA allele frequency and chopstick use is labelled spurious association.
:width: 560px

Ancestry affects both the HLA allele frequency and chopstick use, creating a spurious association between them.
```

Ancestry is a **confounder**: it influences both the allele frequency and the trait, so the two are associated
even though neither affects the other.

## Correcting for ancestry: principal components

The standard solution is to measure each person's genetic ancestry from their genome-wide SNP data, using
**principal components analysis (PCA)**, and to include it in the association model.

```{figure} images/pca_ancestry.png
:alt: PCA of simulated genotypes. PC1 separates African from European and Native American reference individuals; PC2 separates Native American from European. Admixed African American and Hispanic/Latino individuals fall between the reference groups.
:width: 480px

PCA of simulated genotypes (4,000 SNPs) for reference and admixed individuals.
```

PCA summarizes the main axes of genetic variation among the people in a study. Here the first principal
component (PC1) separates African from non-African ancestry and PC2 separates European from Native American
ancestry, while admixed individuals fall in between. **Each person's PC values (PC1, PC2, …) summarize their
genetic ancestry**, computed from genome-wide SNPs. Chapters {doc}`11 <../genetic-diversity-pca/index>` and
{doc}`12 <../genetic-ancestry-admixture/index>` cover PCA and ancestry in depth.

The PCs are then added to the association model as covariates:

$$
Y \sim \beta_0 + \beta_1 x_1 + \beta_2 \,\text{PC}_1 + \beta_3 \,\text{PC}_2 + \dots + \beta_{k+1}\,\text{PC}_k
$$

- $x_1$ is the SNP allele dosage (0, 1 or 2), and $\text{PC}_1 \dots \text{PC}_k$ are the ancestry covariates
  (often the top 10).
- $\beta_1$ now measures the SNP effect **within ancestry**, so spurious signals like *SUSHI* disappear.
- Age and sex enter the model the same way.

This kind of model is a regression, the subject of the next module, {doc}`association-models`.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/covariates.json")
```
