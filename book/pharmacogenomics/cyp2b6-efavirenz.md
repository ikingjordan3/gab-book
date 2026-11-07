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

# CYP2B6 and efavirenz

```{admonition} Module at a glance
**Lecture:** Pharmacogenomics &nbsp;·&nbsp; **Prerequisites:** {doc}`pharmagenes`, {doc}`knowledge-bases`

**After this module you should be able to:** work a pharmacokinetic example end to end — enzyme, star allele, phenotype, starting dose — and explain why allele frequencies make this a population-scale problem.
```

**Efavirenz** is a non-nucleoside reverse transcriptase inhibitor, used for years as first-line therapy
for HIV-1 and still widely used where newer agents are costly. It is cleared predominantly by **CYP2B6**.

```{figure} images/efavirenz_pathway.png
:alt: Efavirenz is converted by the hepatic enzyme CYP2B6 to inactive 8-hydroxy-efavirenz. Normal metabolizers take the standard 600 mg/day dose; poor metabolizers accumulate efavirenz and risk CNS toxicity, so the dose is lowered.
:width: 100%

Efavirenz is cleared by CYP2B6 8-hydroxylation. Loss-of-function alleles slow that clearance, so plasma
efavirenz rises.
```

This is the pharmacokinetic picture from {doc}`pharmacogenomics-defined` applied to one real drug: the
variant does not touch efavirenz's target, it changes how much efavirenz is in the blood. Patients with
slow clearance are over-exposed, and the adverse effects are central-nervous-system effects — vivid
dreams, dizziness, and in some patients hallucinations and depression.

## Assigning the phenotype

CPIC's guideline for CYP2B6 and efavirenz translates diplotypes into the standard phenotype terms from
{doc}`knowledge-bases`.

| Likely phenotype | Genotype | Example diplotypes |
|---|---|---|
| Ultrarapid metabolizer | Two increased function alleles | \*4/\*4, \*22/\*22, \*4/\*22 |
| Rapid metabolizer | One normal and one increased function allele | \*1/\*4, \*1/\*22 |
| Normal metabolizer | Two normal function alleles | \*1/\*1 |
| Intermediate metabolizer | One normal function allele and one decreased or no function allele | \*1/\*6, \*1/\*18, \*4/\*6 |
| Poor metabolizer | Two decreased function alleles, two no function alleles, or one of each | \*6/\*6, \*18/\*18, \*6/\*18 |

Abridged from Desta Z et al. (2019) *Clin Pharmacol Ther* 106:726–733, Table 1.

## The recommendation

| CYP2B6 phenotype | Implications | Therapeutic recommendation | Strength |
|---|---|---|---|
| Ultrarapid or rapid metabolizer | Slightly lower dose-adjusted trough concentrations | Initiate with standard dosing (600 mg/day) | Strong |
| Normal metabolizer | Normal efavirenz metabolism | Initiate with standard dosing (600 mg/day) | Strong |
| Intermediate metabolizer | Higher trough concentrations; increased risk of CNS adverse events | Consider initiating with a decreased dose of 400 mg/day | Moderate |
| Poor metabolizer | Higher trough concentrations; significantly increased risk of CNS adverse events and discontinuation | Consider initiating with a decreased dose of 400 or 200 mg/day | Moderate |

Abridged from Desta Z et al. (2019), Table 2 (children ≥ 40 kg and adults).

Read the poor-metabolizer row carefully. The recommendation is a **lower dose**, not avoidance: the drug
still works, at an exposure that a normal metabolizer would reach on 600 mg. That is the characteristic
shape of a pharmacokinetic recommendation — the genotype moves the dose.

## Why the allele frequency matters

The *CYP2B6\*6* allele, which carries no function, is common in African populations. In 2015 Zimbabwe
switched its national HIV programme from a three-drug cocktail to a cheaper, once-daily single-pill
regimen containing efavirenz. The pharmacokinetics had been published years before: a large fraction of
the population metabolizes efavirenz slowly, and those patients were over-exposed at the standard dose.
Many experienced CNS adverse effects, and stopped taking the drug.

```{admonition} The lesson
:class: important
A dose that is right **on average** can be wrong for a population with elevated pharmacogenomic risk
allele frequencies. The decision was taken on cost and convenience, without allele frequencies in view.
```

This is the bridge into the second half of the chapter. The science was not missing; the population-level
application of it was. Linda Nordling's 2017 feature in *Nature* ("Putting genomes to work in Africa",
*Nature* 544:20–22) reports the Zimbabwean case, and the researcher whose earlier work had shown the
variant's effect.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/cyp2b6-efavirenz.json")
```
