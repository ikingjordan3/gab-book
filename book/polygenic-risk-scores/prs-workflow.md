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

# Calculating PRS: workflow and data

```{admonition} Module at a glance
**Lecture:** Polygenic risk scores 1 &nbsp;·&nbsp; **Prerequisite:** {doc}`prs-defined`

**After this module you should be able to:** describe the PRS workflow and the independent training, validation and testing data it needs.
```

## A guide to PRS analysis

The workflow in this chapter follows the tutorial by Choi, Mak and O'Reilly,
*A guide to performing polygenic risk score analyses* ([Nature Protocols, 2020](https://doi.org/10.1038/s41596-020-0353-1)).
The same group develops **PRSice** and **PRSice-2** ([Choi & O'Reilly, 2019](https://doi.org/10.1093/gigascience/giz082)),
widely used programs for calculating PRS by clumping and thresholding ({doc}`clumping-thresholding`).

The tutorial divides a PRS analysis into five stages:

```{figure} images/prs_workflow.png
:alt: Five stages of a PRS analysis — Data (base and target), Quality control, PRS calculation (LD clumping, shrinkage, P-value thresholds), Test, Validate.
:width: 100%

The PRS analysis workflow, redrawn after Choi, Mak & O'Reilly (2020), Fig. 1.
```

1. **Data.** Base data (GWAS summary statistics, with β or OR weights) and target data (individual-level genotypes,
   and phenotypes to test the score). See {doc}`base-target-data`.
2. **Quality control (processing).** Both data sets are checked as in a GWAS, with extra care for PRS. See {doc}`prs-quality-control`.
3. **PRS calculation.** LD adjustment (clumping), shrinkage of the effect sizes, and P-value thresholding turn
   GWAS results into weights. This is the subject of the second lecture.
4. **Test.** The score is calculated for every person in the target data and tested for association with the trait.
5. **Validate.** The chosen score is tested again, out of sample, in data separate from both base and target.

## Data for calculating PRS

A PRS needs at least two **independent samples**, with no people in common:

| | Base data | Target data |
|---|---|---|
| Sample | Sample 1 | Sample 2 |
| Role | **Training** data | **Testing** data |
| Contents | GWAS summary statistics: variant associations with the trait of interest, variant locations and identifiers, variant association weights | Individual-level genotype (and phenotype) data |
| Size | As large as possible | Can be much smaller (more than 100 people) |

Because the best way to build the score is not known in advance (which P-value threshold? how much shrinkage?),
many candidate scores are usually built and compared. Choosing the best one in the same data used to report its
performance would make it look better than it is, so the target data are often split again:

| | Base data | Target data | |
|---|---|---|---|
| Sample | Sample 1 | Sample 2 | Sample 3 |
| Role | Training data | **Validation** data | **Testing** data |
| Used to | estimate effect sizes | choose the best candidate score | report the chosen score's performance once |

## An example: genome-wide polygenic scores for five diseases

Khera et al. ([2018](https://doi.org/10.1038/s41588-018-0183-z)),
*Genome-wide polygenic scores for common diseases identify individuals with risk equivalent to monogenic mutations*,
built scores for coronary artery disease, atrial fibrillation, type 2 diabetes, inflammatory bowel disease and
breast cancer using exactly this three-sample design:

- **Derivation (training).** Summary statistics from published GWAS, with LD estimated from 503 Europeans in the
  1000 Genomes Project, were used to derive 31 candidate scores per disease: 24 by pruning and thresholding and
  7 with the LDpred algorithm.
- **Validation.** The best candidate was chosen by its AUC in UK Biobank phase 1 (120,280 people).
- **Testing.** The chosen score was assessed once in UK Biobank phase 2 (288,978 people).

| Disease | Discovery GWAS (cases; controls) | Variants in score | Method (tuning parameter) | AUC, validation | AUC, testing |
|---|---|---|---|---|---|
| Coronary artery disease | 60,801; 123,504 | 6,630,150 | LDpred (ρ = 0.001) | 0.81 | 0.81 |
| Atrial fibrillation | 17,931; 115,142 | 6,730,541 | LDpred (ρ = 0.003) | 0.77 | 0.77 |
| **Type 2 diabetes** | **26,676; 132,532** | **6,917,436** | **LDpred (ρ = 0.01)** | **0.72** | **0.73** |
| Inflammatory bowel disease | 12,882; 21,770 | 6,907,112 | LDpred (ρ = 0.1) | 0.63 | 0.63 |
| Breast cancer | 122,977; 105,974 | 5,218 | Pruning and thresholding (r² < 0.2; P < 5 × 10⁻⁴) | 0.68 | 0.69 |

*From Khera et al. (2018), Table 1. AUCs come from logistic models that also include age, sex, genotyping array and
four ancestry principal components.*

The AUC in the testing data matches the AUC in the validation data, so the chosen scores were not over-fitted.
The type 2 diabetes score in this table is in the PGS Catalog as PGS000014 ({doc}`base-target-data`), and its
tuning parameter ρ appears again in {doc}`shrinkage-methods`.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/prs-workflow.json")
```
