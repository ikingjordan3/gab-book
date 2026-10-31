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

# Clumping and thresholding

```{admonition} Module at a glance
**Lecture:** Polygenic risk scores 2 &nbsp;·&nbsp; **Prerequisite:** {doc}`hard-soft-thresholding`

**After this module you should be able to:** describe how clumping and thresholding use linkage disequilibrium and P-values to choose the variants in a PRS.
```

**Clumping and thresholding** (C+T, also called pruning and thresholding) is the classic way to build a PRS, and
the method used by PRSice ({doc}`prs-workflow`). It combines hard thresholding ({doc}`hard-soft-thresholding`) with
a step that removes redundant, correlated SNPs. Khera et al. built 24 of their 31 candidate scores per disease this
way; for type 2 diabetes, an LDpred score did better ({doc}`shrinkage-methods`).

## The algorithm

1. **Rank** all variants by a statistic, usually the GWAS P-value (−log₁₀ P).
2. Take the SNPs above the −log₁₀ P threshold as **index SNPs**, starting with the most significant SNP that has
   not yet been clumped.
3. Form a **clump** around the index SNP from the SNPs that are
   - above the r² threshold for LD with it, and
   - within the physical distance threshold.

   These SNPs are removed; only the index SNP stays.
4. **Repeat** until there are no more index SNPs.

Then **threshold**: keep the index SNPs with P below a chosen cut-off. Typical clumping settings are r² < 0.1 and a
250 kb window; the P-value threshold is chosen empirically, by building scores at several thresholds and comparing
them in validation data.

## Why clump? Linkage disequilibrium

SNPs near each other are often inherited together, so their genotypes are correlated. This is **linkage
disequilibrium (LD)**, measured by r². Correlated SNPs fall into **LD blocks** (haplotype blocks), separated by
recombination hotspots.

```{figure} images/ld_blocks.png
:alt: LD triangle for 20 SNPs in four haplotype blocks; each diamond gives the pairwise r squared between two SNPs, high (red) within blocks and low (white) between them.
:width: 100%

An LD plot in the style of Haploview. Each diamond is the r² between two SNPs: follow the diagonals down from two SNP
numbers to the diamond where they meet. Simulated haplotypes for 20 SNPs; the triangles mark four LD blocks.
```

Every SNP in LD with a causal variant shows an association signal, so a GWAS peak is a tower of correlated SNPs.
Adding all of them to a PRS would count the same signal many times. Clumping keeps one SNP — the most significant —
per signal.

```{figure} images/ld_clump.png
:alt: LocusZoom-style plot of the 20 simulated SNPs by position, coloured by r squared with index SNP 7, above the LD triangle for the same SNPs; recombination peaks separate the four blocks; block 2 SNPs are red (high LD with SNP 7) and block 3 SNPs are dark blue but significant, peaking at SNP 13.
:width: 100%

The same region as a LocusZoom-style plot ({doc}`../complex-traits-gwas/significance-and-plots`) above its LD triangle.
Points are coloured by their r² with the most significant SNP, SNP 7. Recombination-rate peaks (blue line, right
axis) mark the edges of the LD blocks. Simulated: P-values follow from LD with two causal SNPs, 7 and 13.
```

Reading the two panels together:

- **SNP 7** is the most significant SNP, so it becomes the first index SNP.
- The other **block 2** SNPs are in high LD with SNP 7 (red) and carry its signal. With an r² threshold of 0.2 they
  join SNP 7's clump and are dropped.
- The **block 3** SNPs are also significant, but they are in low LD with SNP 7 (dark blue). They reflect a second,
  independent signal, so the most significant of them, **SNP 13**, becomes the next index SNP.
- Only the index SNPs, 7 and 13, go into the PRS.

Here is the same greedy procedure on a small example:

```{code-cell} ipython3
import numpy as np

snps = ["A", "B", "C", "D", "E", "F"]
logp = np.array([9.0, 7.5, 3.0, 6.0, 5.2, 1.0])        # -log10(P) from the GWAS
r2 = np.array([                                         # pairwise LD (r^2)
    [1.0, 0.8, 0.6, 0.0, 0.0, 0.1],
    [0.8, 1.0, 0.5, 0.0, 0.0, 0.0],
    [0.6, 0.5, 1.0, 0.1, 0.0, 0.0],
    [0.0, 0.0, 0.1, 1.0, 0.7, 0.0],
    [0.0, 0.0, 0.0, 0.7, 1.0, 0.0],
    [0.1, 0.0, 0.0, 0.0, 0.0, 1.0]])

def clump(logp, r2, p_threshold, r2_threshold=0.2):
    clumped = np.zeros(len(logp), bool)
    index_snps = []
    for s in np.argsort(-logp):                         # most significant first
        if clumped[s] or logp[s] < -np.log10(p_threshold):
            continue
        index_snps.append(s)
        clumped |= r2[s] > r2_threshold                 # s and every SNP in LD with it
    return index_snps

for threshold in [5e-8, 1e-5, 1e-2]:
    kept = clump(logp, r2, threshold)
    print(f"P < {threshold:g}: index SNPs {[snps[k] for k in kept]}")
```

At P < 5 × 10⁻⁸ only SNP A is kept (B and C are in its clump). Relaxing the threshold to 10⁻⁵ adds SNP D, the lead
SNP of a second signal (E is in its clump). At 10⁻² nothing new is added: C passes the threshold but is already in A's clump, and F (P = 0.1) does not pass it.

## Clumping and thresholding in practice

The figure below applies C+T to the simulated type 2 diabetes GWAS from {doc}`gwas-to-prs`. After clumping (one index
SNP per ±1 Mb window), each circle is an index SNP. The P-value threshold then decides how many enter the score.

```{figure} images/ct_manhattan.png
:alt: Simulated type 2 diabetes Manhattan plot after clumping; 26 index SNPs circled in orange pass P < 5 × 10⁻⁸ and 33 more circled in purple pass P < 10⁻⁵.
:width: 100%

Clumping and thresholding on the simulated type 2 diabetes GWAS. At P < 5 × 10⁻⁸, 26 variants enter the PRS (orange);
at P < 10⁻⁵, 59 do (orange and purple).
```

The extra variants at the more lenient threshold include real signals that were just below genome-wide significance,
but also some chance peaks. Whether adding them improves the score is an empirical question, answered in the
validation data — which is why the P-value threshold is said to be **empirically determined**.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/clumping-thresholding.json")
```
