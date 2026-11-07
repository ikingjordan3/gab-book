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

# CFTR and ivacaftor

```{admonition} Module at a glance
**Lecture:** Pharmacogenomics &nbsp;·&nbsp; **Prerequisites:** {doc}`pharmacogenomics-defined`, {doc}`pharmagenes`

**After this module you should be able to:** work a pharmacodynamic example end to end, where the variant is in the drug's own target and genotype decides whether the drug is prescribed at all.
```

Cystic fibrosis is the monogenic disease you met in {doc}`../monogenic-disease/index`. Here it returns as
a pharmacogenomic problem, because the gene that causes the disease is also the drug's target.

## CFTR is a gated chloride channel

**CFTR** — the cystic fibrosis transmembrane conductance regulator — is an ATP-gated chloride channel in
the apical membrane of epithelial cells.

```{figure} images/cftr_channel.png
:alt: Two epithelial membrane panels. A typical CFTR channel opens and chloride leaves the cell. A G551D channel is present in the membrane but its gate stays shut, and chloride accumulates in the cytoplasm.
:width: 100%

CFTR with and without the G551D variant. Schematic.
```

The **G551D** variant is a **gating** defect. The protein folds, traffics, and sits in the membrane
correctly — it simply opens only rarely. Hold on to that distinction, *present but shut* versus *never
arriving*: it is what decides whether a drug can help.

Without chloride transport, water does not follow, the airway surface liquid is depleted, and the mucus
layer becomes thick and sticky. Cilia cannot clear it, bacteria are trapped, and chronic infection
damages the lung.

```{figure} images/cf_airway.png
:alt: Two airway panels. In a typical airway a watery periciliary layer keeps the cilia upright and thin mucus is swept clear. In a cystic fibrosis airway the layer is collapsed, the cilia lie flat, and thick mucus traps bacteria.
:width: 100%

The consequence of a shut chloride channel at the airway surface. Schematic.
```

Cystic fibrosis also affects the pancreas, intestine and sweat glands. Sweat chloride remains the standard
diagnostic measurement, and it is the outcome measure used in the trials below.

## Six classes of defect need six kinds of drug

| Class | Nature of defect | Example | Treatment strategy |
|---|---|---|---|
| Class I | Biosynthesis | G542X | Suppress premature termination codon |
| Class II | Folding and trafficking | F508del | Restore folding |
| **Class III** | **Gating** | **G551D** | **Restore gating** |
| Class IV | Conductance | R117H | Increase pore size |
| Class V | Reduced splicing efficiency | 5T, 2789 + 5 G>A | Increase splicing efficiency |
| Class VI | Shortened time at plasma membrane | F508del (corrected) | Stabilize protein at plasma membrane |

Source: Clancy JP et al. (2014) *Clin Pharmacol Ther* 95:592–597, Table 1.

## Ivacaftor is a potentiator

**Ivacaftor** binds CFTR directly and increases the probability that a channel already at the cell
surface is open. It addresses class III — gating.

```{figure} images/ivacaftor_gating.png
:alt: An untreated G551D channel sits in the membrane with its gate shut and chloride accumulating inside the cell. With ivacaftor bound, the gate is held open and chloride leaves the cell.
:width: 100%

Ivacaftor as a potentiator. Schematic, redrawn after the figure in the CPIC guideline.
```

It follows directly that ivacaftor cannot help a channel that never reaches the membrane. F508del, the
commonest CF variant, is primarily a folding defect, and ivacaftor alone does nothing for those patients.
The prescribing table says exactly that:

| CFTR genotype | Implications for ivacaftor effects | Therapeutic recommendation | Strength |
|---|---|---|---|
| Homozygous or heterozygous G551D | Significant improvement in lung function, weight and sweat chloride, through enhanced channel activity | Use ivacaftor according to the product label | Strong |
| Noncarrier of G551D | Not studied clinically; in vitro effects vary by variant and the likelihood of response is not known | Ivacaftor is not recommended | Moderate |
| Homozygous for F508del | No significant reduction in sweat chloride or change in other clinical measurements; unlikely to respond | Ivacaftor is not recommended | Moderate |

Abridged from Clancy JP et al. (2014), Table 2.

## The two routes, side by side

Compare this table with the efavirenz table in {doc}`cyp2b6-efavirenz`. There, genotype moved the **dose**.
Here, genotype decides whether the drug is **given at all** — because the variant is in the drug's own
target, not in the machinery that clears it.

```{admonition} Since the guideline
:class: note
Ivacaftor was approved in 2012 for G551D, the first drug approved for a specific CFTR genotype, at a cost
of roughly $300,000 per patient per year. Its label has since broadened to further responsive variants,
and combination therapies — elexacaftor/tezacaftor/ivacaftor, pairing correctors with the potentiator —
now cover F508del homozygotes. The CPIC logic is unchanged: match the drug to the molecular defect.
```

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/cftr-ivacaftor.json")
```
