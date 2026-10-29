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

# Association models

```{admonition} Module at a glance
**Lecture:** Genome-wide association studies &nbsp;·&nbsp; **Prerequisite:** {doc}`covariates`

**After this module you should be able to:** write the linear and logistic regression models used in GWAS, interpret the SNP coefficient β₁ (including e<sup>β₁</sup> as an odds ratio), and describe the summary statistics a GWAS produces.
```

{doc}`covariates` showed that a GWAS must account for age, sex and ancestry. The χ² test cannot do that, so
GWAS use **regression** instead, which can include any number of covariates. The kind of regression depends on
the trait:

| Trait | Example | Model |
|---|---|---|
| **Continuous** (quantitative) | height, blood glucose, BMI | **Linear regression** |
| **Binary** (case-control) | type 2 diabetes yes/no | **Logistic regression** |

In both models, each SNP is coded by its **allele dosage** $x_1$: the number of copies (0, 1 or 2) of the effect
allele a person carries. Coding the genotype as 0, 1, 2 assumes an **additive** model, in which each extra copy
of the allele has the same effect, which matches the additive variance of {doc}`heritability`.

## Linear regression

A continuous, quantitative trait value ($Y$) is modelled **linearly** for a SNP:

$$
Y \sim \beta_0 + \beta_1 x_1 + \beta_2 C_1 + \dots + \beta_n C_{n-1}
$$

where $x_1$ is the SNP allele dosage {0, 1, 2} and $C_a$ is a covariate (age, sex, PC1, PC2, …).

```{figure} images/regression_linear.png
:alt: Simulated trait values for people with genotypes TT, CT and CC (dosage 0, 1, 2), with a fitted regression line rising across genotypes; the slope is beta 1.
:width: 460px

Linear regression of a simulated quantitative trait on allele dosage. The slope is β₁.
```

- $\beta_1$ is the **change in the trait for each additional copy** of the effect allele, holding the covariates
  constant. For example, β₁ = 0.05 mmol/L for fasting glucose means each copy of the allele raises glucose by
  0.05 mmol/L on average.
- $\beta_0$ is the intercept, and $\beta_2 \dots \beta_n$ are the covariate effects, which are estimated but not
  usually reported.

## Logistic regression

In a logistic regression we model a **binary** response variable ($Y$), specifying $p = P(Y = 1)$, the probability
of being a case. Probabilities are bounded between 0 and 1, so instead of modelling $p$ directly, we work with the
odds,

$$
\text{Odds} = \frac{p}{1-p}
$$

and assume a linear relationship between the predictor variables and the **log-odds** of the event that $Y = 1$:

$$
\ln \frac{p}{1-p} \sim \beta_0 + \beta_1 x_1 + \beta_2 C_1 + \dots + \beta_n C_{n-1}
$$

where $x_1$ is the SNP allele dosage and $C_a$ is a covariate. The log-odds uses the natural logarithm, as in all
standard GWAS software.

```{figure} images/regression_logistic.png
:alt: Logistic curve of the probability of being a case against the log-odds. Points mark TT, CT and CC genotypes; each C allele moves one step along the curve, multiplying the odds by the odds ratio.
:width: 460px

Logistic regression: each copy of the C allele adds β₁ to the log-odds, which multiplies the odds by e<sup>β₁</sup>.
```

The key to interpreting $\beta_1$: adding β₁ to the log-odds **multiplies the odds by $e^{\beta_1}$**. So

$$
e^{\beta_1} = \text{OR per copy of the effect allele}
$$

the same odds ratio as in {doc}`odds-ratios`, but now adjusted for the covariates. A β₁ of 0 means no association
(OR = 1); a positive β₁ means the allele increases risk.

**Check it with the example data.** The genotype counts from {doc}`case-control-design` gave an allelic odds ratio
of 1.70. Fitting a logistic regression of case status on allele dosage gives the same answer. The few lines below
fit the model with Newton's method (GWAS software such as PLINK does the same thing, very fast, at every SNP):

```{code-cell} ipython3
import numpy as np
from scipy.stats import norm

# genotype counts from the case-control example: dosage = number of C alleles
counts = {"cases":    {2: 384, 1: 472, 0: 144},
          "controls": {2: 240, 1: 500, 0: 260}}
x, y = [], []
for status, group in [(1, "cases"), (0, "controls")]:
    for dosage, n in counts[group].items():
        x += [dosage] * n
        y += [status] * n
X = np.column_stack([np.ones(len(x)), x])        # intercept + dosage
y = np.array(y, dtype=float)

beta = np.zeros(2)
for _ in range(25):                              # Newton-Raphson for logistic regression
    p = 1 / (1 + np.exp(-X @ beta))
    H = X.T @ (X * (p * (1 - p))[:, None])
    beta += np.linalg.solve(H, X.T @ (y - p))

se = np.sqrt(np.diag(np.linalg.inv(H)))          # standard errors
z = beta[1] / se[1]
print(f"beta1 = {beta[1]:.3f}, SE = {se[1]:.3f}, OR = exp(beta1) = {np.exp(beta[1]):.2f}, "
      f"P = {2 * norm.sf(abs(z)):.1e}")
```

$e^{\beta_1} = 1.70$, matching the χ² calculation. The regression also gives a **standard error** for β₁, and from
the ratio $z = \beta_1 / \text{SE}$ a P-value. Adding covariates would simply add columns to `X`.

## GWAS statistical analysis

A GWAS consists of running a statistical association analysis of **each individual SNP** with a trait (disease)
of interest, controlling for covariates. The model is the same at every SNP; only $x_1$ changes. GWAS results
are the **statistical coefficient values** for the association models run for each SNP.

## Summary statistics

The collected results are the GWAS **summary statistics**: a summary of all associations between variants (SNPs)
and the trait, with the results of the statistical analysis for each SNP genome-wide. They are used to derive
polygenic risk scores ({doc}`gwas-to-prs`). Each row describes one SNP:

| Variant (SNP) location information | GWAS analysis results (statistical coefficients) |
|---|---|
| Chromosome | Effect size: β (or OR = e<sup>β</sup>) |
| Position | Standard error |
| rsID | P-value |
| Effect allele / other (reference/alternate) allele | |

Real summary statistics also usually include the effect-allele frequency and the sample size.
{doc}`significance-and-plots` shows what such a table looks like for type 2 diabetes and how it is turned into
Manhattan and Q-Q plots.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/association-models.json")
```
