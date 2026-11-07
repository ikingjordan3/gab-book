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

# Pharmagenes

```{admonition} Module at a glance
**Lecture:** Pharmacogenomics &nbsp;·&nbsp; **Prerequisite:** {doc}`pharmacogenomics-defined`

**After this module you should be able to:** name the four classes of pharmacogene, read star-allele nomenclature, translate a diplotype into a metabolizer phenotype, and classify effects as efficacy, dosage or toxicity.
```

Genes with variants that change drug response are called **pharmacogenes**, or pharmagenes. PharmGKB
sorts them into four classes.

| Class | What they do | Examples |
|---|---|---|
| **Drug transporters** | Cell-surface proteins that move biotics and xenobiotics across the membrane. Variation raises or lowers the concentration of a compound inside and outside the cell. | influx: *SLCO1B1*; efflux: *ABCB1*, *ABCG2* |
| **Drug metabolizing enzymes** | Control the chemical conversions that degrade drugs and allow them to be excreted, mainly in the liver. They can inactivate a drug — or activate a prodrug. | phase 1 (cytochrome P450): *CYP2C9*, *CYP2D6*, *CYP2C19*, *CYP3A5*; phase 2: *TPMT*, *UGT1A1* |
| **Drug targets** | The molecules or pathways a drug is designed to affect. Variation changes how well a drug works, or restricts therapy to particular genotypes. | *CFTR*, *VKORC1* |
| **MHC / HLA genes** | Cell-surface proteins that present antigens to T cells. Some alleles carry a high risk of hypersensitivity to particular drugs. | *HLA-B\*57:01* (abacavir), *HLA-B\*15:02* (carbamazepine) |

The first two classes act on drug exposure — pharmacokinetics. The last two act on drug effect and on the
immune response — pharmacodynamics.

```{admonition} A link to cancer genomics
:class: seealso
HLA typing is the same step that underlies neoantigen prediction in {doc}`../cancer-genomics/index`.
The genotype matters there for a different reason: not hypersensitivity, but which mutant peptides a
patient's own immune system can be shown.
```

## Star alleles: naming the haplotypes of a pharmagene

Pharmacogenomics has a nomenclature of its own, and it is unlike anything else in this course. A
**star allele** names a **haplotype** of a pharmacogene: a particular combination of variants carried
together on one copy of the gene. By convention `*1` is the reference haplotype.

```{figure} images/star_alleles.png
:alt: Three CYP2C19 haplotypes drawn along the gene. *1 is the reference with normal function; *2 carries c.681G>A (rs4244285), a splice defect, and has no function; *17 carries the promoter variant c.-806C>T (rs12248560) and has increased function. Each haplotype also carries other variants, shown in grey. Below, the diplotype CYP2C19*1/*2 gives an intermediate metabolizer.
:width: 100%

Star alleles are haplotypes, not single SNPs. Each allele is assigned a function; because we carry two
copies of each gene, a person has a **diplotype**, which is what guideline tables are written in.
```

Three points are worth keeping straight:

- A star allele is **not** a single SNP. Each haplotype carries many variants; one or a few of them alter
  function, and the rest come along for the ride.
- Each allele is assigned a **function**: normal, decreased, no function, or increased.
- The same star number means different things in different genes. An allele is only ever named together
  with its gene: *CYP2C19\*2* and *CYP2B6\*6* have nothing to do with each other.

## From diplotype to metabolizer phenotype

A diplotype is translated into a **metabolizer phenotype** by combining the functions of the two alleles.
The phenotype is what a prescribing guideline acts on.

```{figure} images/metabolizer_curves.png
:alt: Plasma drug concentration against time after one standard dose for ultrarapid, normal, intermediate and poor metabolizers. The poor metabolizer curve stays above the therapeutic window; the ultrarapid curve falls below it.
:width: 680px

The same dose produces different exposure in different metabolizer phenotypes. Illustrative
one-compartment model ($k_a$ = 1.2 /h), with clearance varying by phenotype; the therapeutic window is
schematic.
```

For a drug that the enzyme **inactivates**, a poor metabolizer accumulates the drug and risks toxicity,
while an ultrarapid metabolizer clears it quickly and risks treatment failure. For a **prodrug**, which
the enzyme activates, the picture reverses: the poor metabolizer gets too little active compound.
Clopidogrel and *CYP2C19* is the classic prodrug example.

## Three aspects of drug response

Pharmacogenomic variants in pharmagenes are associated with three aspects of drug response, and these are
the categories used to classify them in the second half of this chapter:

- **Efficacy** — does the drug work in this patient?
- **Dosage** — how much drug does this patient need?
- **Toxicity** — will this patient be harmed?

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/pharmagenes.json")
```
