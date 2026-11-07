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

# Pharmacogenomics defined

```{admonition} Module at a glance
**Lecture:** Pharmacogenomics &nbsp;·&nbsp; **Prerequisite:** {doc}`../complex-traits-gwas/variation-spectrum`

**After this module you should be able to:** place drug response on the frequency × effect-size spectrum, define pharmacogenomics, and distinguish pharmacokinetic from pharmacodynamic mechanisms.
```

## A corner of the spectrum we have treated as empty

{doc}`../complex-traits-gwas/variation-spectrum` plotted disease-associated variants by allele frequency
against effect size. Variants trade one against the other along a diagonal band: the rare alleles that
cause Mendelian disease have large effects, and the common variants found by GWAS have small ones. The
top-right corner — common **and** high-effect — was nearly empty.

Pharmacogenomic variants sit in that corner.

```{figure} images/spectrum_pgx.png
:alt: Allele frequency against effect size. The common, high-effect class, labelled "response to drugs (medication)", is highlighted; the rare Mendelian, low-frequency intermediate, common GWAS and rare small-effect classes are shown in grey.
:width: 680px

Drug response on the frequency × effect-size spectrum. Conceptual; redrawn after Manolio et al. (2009).
```

No-function alleles of *CYP2D6*, *CYP2C19* and *SLCO1B1* reach frequencies of a few percent to twenty
percent in many populations, and they still change drug exposure several-fold. High frequency *and* high
penetrance is unusual for disease risk; for drug response it is the norm. Two consequences follow, and
they organise the rest of this chapter:

1. Because the effects are large, a genotype can change a prescribing decision for an individual patient.
2. Because the variants are common, their frequencies differ among populations — the subject of the
   second half of this chapter.

## Drugs and medications

A **drug** is a chemical substance that produces a biological effect when administered to a living
organism, other than a nutrient or essential dietary ingredient. A **medication** is a drug used to
diagnose, treat, cure or prevent disease. Every medication is a drug; not every drug is a medication.
Pharmacogenomics concerns medications: drugs prescribed to a patient, at a specific dose, with an
intended effect.

**Pharmacogenomics (PGx)** is the study of how a person's DNA affects their response to drugs. It
combines the study of drugs (pharmacology) with the study of genes and their functions (genomics).

## The pharmacogenomics paradigm

A clinical population that looks homogeneous — one diagnosis, one prescription — is in fact a mixture of
response groups.

```{figure} images/pgx_paradigm.png
:alt: A group of patients with the same condition, each drawn in the colour of the response they will have. DNA profiling sorts them into three groups of unequal size: a large good-response group, a small no-response group and a small adverse-reaction group.
:width: 100%

The pharmacogenomics paradigm. The response groups exist before anyone is treated; genotyping says which
group a patient is in. Schematic; group sizes are illustrative.
```

Notice what the figure claims and what it does not. The subgroups are already there, mixed through the
starting population and invisible to the clinician. What genotyping adds is the ability to tell them
apart **in advance**, rather than discovering the answer from the outcome. Notice the sizes too: the
majority usually do respond, which is why the drug reached the market at all. The clinical value of
pharmacogenomics lies in identifying the minorities who will not respond, or who will be harmed.

## Two mechanisms

Genetic variation can change drug response in two distinct ways, and the distinction runs through the
rest of this chapter.

**Pharmacokinetics** — what the body does to the drug. Variants affect the enzymes and transporters that
regulate the rate at which drugs are absorbed and metabolized.

```{figure} images/pharmacokinetics.png
:alt: Four-step flow — absorption, distribution, metabolism and elimination — with transporters acting at absorption, distribution and elimination, and CYP enzymes at metabolism.
:width: 100%

Pharmacokinetics: absorption, distribution, metabolism and elimination. Variation in transporters changes
how much drug reaches the circulation and the tissues; variation in metabolizing enzymes changes how fast
it is cleared.
```

**Pharmacodynamics** — what the drug does to the body. Variants affect the drug's target itself, so the
same drug concentration produces a different effect.

```{figure} images/pharmacodynamics.png
:alt: Two membrane panels. At a typical drug target a triangular drug nests into a matching V-shaped cleft and a signal follows. At a variant target the cleft is flat-bottomed, the same drug cannot seat, and there is no response.
:width: 100%

Pharmacodynamics: a variant in the drug target changes the response. Schematic.
```

Pharmacokinetic variation acts on drug **exposure**; pharmacodynamic variation acts on drug **effect**.
{doc}`cyp2b6-efavirenz` works through an example of the first, {doc}`cftr-ivacaftor` an example of the
second.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/pharmacogenomics-defined.json")
```
