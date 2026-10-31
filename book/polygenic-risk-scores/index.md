# Polygenic risk scores

A genome-wide association study (GWAS) finds hundreds of variants that each change disease risk a little.
A **polygenic risk score (PRS)** adds those small effects up into one number per person. This chapter follows
the two lectures on polygenic risk scores, using **type 2 diabetes** as the running example, the same disease
as the GWAS examples in {doc}`../complex-traits-gwas/index`.

The first lecture, *Polygenic risk scores 1*, defines a PRS and shows what it can and cannot tell us. It then covers the
data a PRS needs (base data from a GWAS and target data from the people being scored), the public catalogs that hold them, and
quality control. The second lecture, *Polygenic risk scores 2*, covers the step in the middle: how GWAS effect sizes are turned
into PRS weights by shrinkage, thresholding and clumping.

Work through the modules in order. Each one ends with a short self-check quiz.

| Lecture | Module | After this module you should be able to… |
|---|---|---|
| PRS 1 | {doc}`From GWAS to polygenic risk scores <gwas-to-prs>` | Explain how GWAS summary statistics become the input for a polygenic risk score. |
| PRS 1 | {doc}`Polygenic risk scores defined <prs-defined>` | Define a PRS, write it as a weighted sum of effect alleles, and distinguish the effect allele from the risk allele. |
| PRS 1 | {doc}`Risk stratification, not diagnosis <risk-stratification>` | Explain why a PRS stratifies risk rather than diagnoses disease, and distinguish relative risk, absolute risk and AUC. |
| PRS 1 | {doc}`Calculating PRS: workflow and data <prs-workflow>` | Describe the PRS workflow and the independent training, validation and testing data it needs. |
| PRS 1 | {doc}`Base and target data <base-target-data>` | Use the PGS Catalog and GWAS Catalog, read summary statistics with either allele convention, and calculate a PRS by hand. |
| PRS 1 | {doc}`Quality control <prs-quality-control>` | List the quality-control checks base and target data need before a PRS is calculated, and why each matters. |
| PRS 2 | {doc}`Shrinkage and the winner's curse <shrinkage-winners-curse>` | Explain why GWAS effect sizes are over-estimated and why they must be shrunk before use in a PRS. |
| PRS 2 | {doc}`Hard and soft thresholding <hard-soft-thresholding>` | Contrast P-value thresholding (hard) with shrinking all effects (soft). |
| PRS 2 | {doc}`Shrinkage methods <shrinkage-methods>` | Name the main shrinkage methods (LASSO, ridge, LDpred) and explain why GWAS power matters more than the method. |
| PRS 2 | {doc}`Clumping and thresholding <clumping-thresholding>` | Describe how clumping and thresholding use linkage disequilibrium and P-values to choose the variants in a PRS. |
