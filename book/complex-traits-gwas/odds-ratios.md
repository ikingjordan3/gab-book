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

# Effect size: odds ratios and P-values

```{admonition} Module at a glance
**Lecture:** Genome-wide association studies &nbsp;·&nbsp; **Prerequisite:** {doc}`case-control-design`

**After this module you should be able to:** compute odds, odds ratios and a χ² P-value for one SNP.
```

A case-control GWAS asks the same question at every genotyped SNP: *is one allele more common in cases
than in controls?* For each SNP we want two numbers:

- an **effect size**, which says how strongly the allele is associated with the disease, and
- a **P-value**, which says how surprising an association that strong would be if the allele had no effect.

For a binary trait (case or control), the effect size is usually reported as an **odds ratio**.

## Odds

The **odds** of an event is the probability that it happens divided by the probability that it does not:

$$
\text{odds} = \frac{p}{1-p}
$$

For a biallelic SNP with alleles C and T, the odds of carrying C on a chromosome is the ratio of the two
allele frequencies:

$$
\text{odds} = \frac{\text{freq}(C)}{\text{freq}(T)}
$$

## Odds ratio

The **odds ratio (OR)** compares the odds in two groups. In a case-control GWAS it is the odds of the
effect allele in cases divided by the odds of the same allele in controls:

$$
\text{OR} = \frac{\text{freq}_\text{cases}(C) \,/\, \text{freq}_\text{cases}(T)}
                 {\text{freq}_\text{controls}(C) \,/\, \text{freq}_\text{controls}(T)}
$$

An OR of 1 means no association. An OR above 1 means the allele is more common in cases, and an OR below
1 means it is more common in controls.

```{figure} images/allele_frequencies.png
:alt: Allele frequencies at one SNP. In 1,000 cases, C is 62% of alleles and T is 38%. In 1,000 controls, C is 49% and T is 51%.
:width: 520px

Allele frequencies at one SNP in 1,000 cases and 1,000 controls.
```

**Worked example.** In cases, C is at 62% and T at 38%, so the odds are 0.62 / 0.38 = 1.63.
In controls, C is at 49% and T at 51%, so the odds are 0.49 / 0.51 = 0.96. Then

$$
\text{OR} = \frac{0.62 / 0.38}{0.49 / 0.51} = \frac{1.63}{0.96} = 1.70
$$

so the odds of carrying C are 1.7 times higher in cases than in controls. C is the **effect allele**.

## Is the difference significant? The χ² test

An OR of 1.70 could still arise by chance in a small study. A **χ² (chi-square) test** compares the allele
counts we observed with the counts we would expect if cases and controls had the same allele frequency:

$$
\chi^2 = \sum \frac{(\text{Observed} - \text{Expected})^2}{\text{Expected}}
$$

Each person carries two alleles, so 1,000 cases and 1,000 controls give 2,000 alleles per group.
Expected counts come from the table margins. For example, the expected number of C alleles in cases is
2,000 × 2,220 / 4,000 = 1,110.

| | C | T | Total |
|---|---|---|---|
| Cases | 1,240 (1,110) | 760 (890) | 2,000 |
| Controls | 980 (1,110) | 1,020 (890) | 2,000 |
| Total | 2,220 | 1,780 | 4,000 |

*Observed (expected) allele counts.*

$$
\chi^2 = 2\left(\frac{130^2}{1110}\right) + 2\left(\frac{130^2}{890}\right) = 68.4
$$

With 1 degree of freedom, χ² = 68.4 gives P ≈ 1 × 10⁻¹⁶. You can check this yourself:

```{code-cell} ipython3
from scipy.stats import chi2_contingency

table = [[1240, 760],    # cases:    C, T
         [980, 1020]]    # controls: C, T
chi2, p, dof, expected = chi2_contingency(table, correction=False)
odds_ratio = (1240 / 760) / (980 / 1020)
print(f"OR = {odds_ratio:.2f}, chi2 = {chi2:.1f}, df = {dof}, P = {p:.1e}")
```

`correction=False` turns off Yates' continuity correction so the result matches the hand calculation.

## The limitation: no covariates

The χ² test has one big drawback for GWAS: **it cannot include covariates**. Age, sex and population
structure can all change allele frequencies or disease risk, and if they differ between cases and
controls they can create associations that are not real. The next module, {doc}`covariates`, explains
why, and {doc}`association-models` replaces the χ² test with regression, which can adjust for them.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/odds-ratios.json")
```
