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

# Significance and visualization

```{admonition} Module at a glance
**Lecture:** Genome-wide association studies &nbsp;·&nbsp; **Prerequisite:** {doc}`association-models`

**After this module you should be able to:** read a Q-Q plot and the genomic inflation factor λ, explain the genome-wide significance threshold of 5 × 10⁻⁸, and read GWAS summary statistics, Manhattan plots and LocusZoom plots.
```

A GWAS produces a P-value for every SNP, often millions of them. This module covers how to check that those
P-values can be trusted, how to decide which ones are significant, and how to display the results. The examples
use a simulated type 2 diabetes GWAS with peaks at real, well-known T2D loci.

## Q-Q plots: checking the whole study

A **quantile-quantile (Q-Q) plot** compares the observed GWAS P-values (y-axis) with the P-values expected under
the null hypothesis of no association (x-axis), both on a −log₁₀ scale. If no SNP were associated, the points
would follow the diagonal.

```{figure} images/qq_single.png
:alt: Simulated Q-Q plot. Points follow the diagonal within a grey 95% band, then a tail of points rises above the diagonal and past the genome-wide significance line. Lambda equals 1.00.
:width: 460px

Q-Q plot of a simulated, well-controlled GWAS with true associations (λ = 1.00).
```

- Q-Q plots are used to identify **systematic biases**, such as population structure or genotyping error.
- **No deviation** indicates a well-controlled study; deviations at the top end are due to true associations
  (λ ≈ 1).
- **Deviation with λ ≫ 1** indicates a high number of false positives.
- **Ideal:** no deviation at the low end, with a tail of deviations at the high end.

### The genomic inflation factor, λ

The overall deviation is summarized by the **genomic inflation factor**:

$$
\lambda = \frac{\text{median}(\chi^2)}{0.455}
$$

the median observed χ² statistic divided by its expected median under the null. For a χ² with 1 degree of
freedom, that expected median is 0.455. The median is used because most SNPs are null, so a handful of true
associations barely move it, whereas confounding shifts every SNP.

```{code-cell} ipython3
import numpy as np
from scipy.stats import chi2

rng = np.random.default_rng(11)
null = rng.chisquare(1, 500_000)                       # 500,000 null SNPs
print(f"expected median of chi2 (1 df) = {chi2.ppf(0.5, 1):.3f}")
print(f"well-controlled: lambda = {np.median(null) / 0.455:.2f}")
print(f"stratified:      lambda = {np.median(1.2 * null) / 0.455:.2f}")   # every statistic inflated by 1.2
```

### Inflation vs. true signal

```{figure} images/qq_panels.png
:alt: Four simulated Q-Q plots. Panel 1, well-controlled with no signal, lambda 1.00, points on the diagonal. Panel 2, bona fide hits, lambda 1.00, points on the diagonal with a tail rising at the top. Panel 3, inflation from stratification, lambda 1.20, the whole curve lifts off the diagonal early. Panel 4, after PC correction, lambda 1.00, recovers panel 2.
:width: 100%

Simulated Q-Q plots for 500,000 SNPs. Panels 2–4 include 150 truly associated SNPs; panel 3 multiplies every
test statistic by 1.2 to mimic population stratification.
```

**Inflation lifts the whole curve (λ > 1); true hits lift only the tail (λ ≈ 1).**

1. No true associations: all points lie within the 95% band.
2. The same null SNPs plus 150 truly associated SNPs: the curve stays on the diagonal and only the extreme tail
   rises past the significance threshold; λ is unchanged because the median SNP is null.
3. Population stratification inflates every test statistic, so the curve leaves the diagonal early and
   λ = 1.20. Even null SNPs drift toward significance: false positives.
4. Adding ancestry principal components as covariates ({doc}`covariates`) removes the inflation and recovers
   panel 2.

One caveat: in very large GWAS of highly polygenic traits, λ can exceed 1 **without any confounding**, because
thousands of SNPs carry small true effects. LD score regression ({doc}`gwas-to-prs`) separates polygenic signal
(the slope) from confounding (the intercept).

## Statistical significance thresholds

P-value thresholds are chosen to balance the need to detect true associations (avoiding Type II errors) against
the need to avoid false positive associations (Type I errors).

- Hundreds of thousands to millions of variants are tested, so the probability of Type I errors increases
  dramatically: at P < 0.05, one null SNP in 20 would be called significant.
- Not all variants (SNPs) are independent, owing to linkage disequilibrium (LD); there are estimated to be
  **about 1 million independent tests** in the genome of a population of European ancestry.
- A **Bonferroni correction** divides the overall error rate by the number of tests:

$$
p = \frac{0.05}{1 \times 10^6} = 5 \times 10^{-8}
$$

This is the conventional **genome-wide significance** threshold. Less stringent **suggestive** cutoffs
(P = 10⁻⁵ or 10⁻⁶) are used in discovery cohorts, to pick SNPs for follow-up, or in polygenic risk score
development, where including many weaker SNPs can improve prediction.

## GWAS summary statistics: what the results look like

Each row of the summary statistics ({doc}`association-models`) is one SNP. Here are the top lead SNPs, sorted by
P-value:

| CHR | POS | rsID | Gene | EA | OA | EAF | β | SE | P |
|---|---|---|---|---|---|---|---|---|---|
| 10 | 112,998,590 | rs7903146 | *TCF7L2* | T | C | 0.29 | 0.315 | 0.0243 | 2.1 × 10⁻³⁸ |
| 6 | 20,679,478 | rs7756992 | *CDKAL1* | G | A | 0.28 | 0.157 | 0.0162 | 3.4 × 10⁻²² |
| 16 | 53,786,615 | rs9939609 | *FTO* | A | T | 0.41 | 0.122 | 0.0139 | 1.6 × 10⁻¹⁸ |
| 3 | 185,793,899 | rs4402960 | *IGF2BP2* | T | G | 0.32 | 0.113 | 0.0140 | 4.8 × 10⁻¹⁶ |
| 8 | 117,172,544 | rs13266634 | *SLC30A8* | C | T | 0.70 | 0.113 | 0.0143 | 2.7 × 10⁻¹⁵ |
| 3 | 12,351,626 | rs1801282 | *PPARG* | C | G | 0.88 | 0.122 | 0.0179 | 7.9 × 10⁻¹² |

*The full table continues for about 1.9 million more rows, one per SNP tested.*

*Illustrative: real lead SNPs and alleles (GRCh38; EA = effect allele, the T2D risk allele; OA = other allele;
EAF = effect-allele frequency). EAF, β, SE and P are simulated to match the Manhattan plot below.*

- **CHR + POS** give the position on the x-axis of a Manhattan plot.
- **−log₁₀(P)** gives the height on the y-axis.
- **β (ln OR) and SE** give the effect size and its precision. For example, *TCF7L2* β = 0.315 corresponds to
  OR = e<sup>0.315</sup> ≈ 1.37 per T allele, and z = β/SE = 13.

## Manhattan plots

A **Manhattan plot** shows the −log₁₀-transformed P-values for each SNP–trait association genome-wide; higher
values are more significant.

```{figure} images/manhattan.png
:alt: Simulated type 2 diabetes Manhattan plot. Points for each SNP across chromosomes 1 to 22, coloured alternately by chromosome. Peaks rise above the dashed genome-wide significance line at TCF7L2 (the tallest, on chromosome 10), CDKAL1, FTO, IGF2BP2, SLC30A8, PPARG and JAZF1.
:width: 100%

Simulated T2D GWAS with peaks placed at well-known T2D loci (approximate GRCh38 positions). The dashed line is
the genome-wide significance threshold, 5 × 10⁻⁸.
```

Each point is one SNP, plotted by genomic position (x) and −log₁₀ P (y). Each "tower" is a **locus**: many
correlated SNPs in LD with a causal variant rise together, which is what gives the plot its skyline. The labelled
peaks are real, well-replicated T2D loci, including *TCF7L2*, the strongest common-variant T2D signal. Real T2D
GWAS such as the multi-ancestry meta-analysis of 2.5 million people find more than 600 loci (Suzuki et al.,
[2024](https://doi.org/10.1038/s41586-024-07019-6)).

## LocusZoom plots

A **LocusZoom plot** is a local Manhattan plot. It shows the **lead SNP** at a given locus, here rs7903146 in
*TCF7L2*, and the other SNPs in linkage disequilibrium with it.

```{figure} images/locuszoom.png
:alt: Simulated LocusZoom plot of the TCF7L2 region on chromosome 10. The lead SNP rs7903146 is a purple diamond at the top. Nearby SNPs are coloured by r squared with the lead SNP, red for high LD through dark blue for low LD, and their significance falls off with distance. A blue line shows recombination rate, with peaks at the edges of the associated region. Gene tracks for VTI1A and TCF7L2 are shown below.
:width: 620px

Simulated LocusZoom plot at *TCF7L2*. Gene positions are approximate (GRCh38).
```

- **Colour** is the r² of each SNP with the lead SNP. Highly correlated SNPs (red and orange) carry almost the
  same signal as the lead SNP, which is why GWAS identify **loci** rather than single causal variants.
- The **blue line** is the recombination rate. Recombination hotspots mark the edges of the LD block, beyond
  which the association drops off.
- The **gene track** shows which genes lie in the region. The lead SNP here sits in an intron of *TCF7L2*.

Going from an associated locus to the causal variant and gene (**fine-mapping** and functional follow-up) is a
research field in its own right. The lab, {doc}`lab-single-snp-association`, lets you run an association test
yourself.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/significance-and-plots.json")
```
