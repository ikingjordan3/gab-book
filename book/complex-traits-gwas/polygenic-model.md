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

# The polygenic model of complex traits

```{admonition} Module at a glance
**Lecture:** Genetics of complex, common disease &nbsp;·&nbsp; **Prerequisite:** {doc}`variation-spectrum`

**After this module you should be able to:** explain how many small genetic effects add up to a continuous, normally distributed trait, and how a liability threshold turns that into a binary disease that clusters in families.
```

Mendel's peas were either yellow or green, round or wrinkled. Height, blood pressure and blood glucose
are not like that: they vary continuously, and most people sit near the middle of a bell-shaped
distribution. How can discrete genes produce continuous traits? Answering that question, more than a
century ago, gave us the model that underlies GWAS and polygenic risk scores today.

## Mendelians vs. biometricians (1900–1918)

After Mendel's work was rediscovered in 1900, two camps disagreed about heredity.

| | Mendelians | Biometricians |
|---|---|---|
| **Leading figures** | William Bateson, Hugo de Vries, Reginald Punnett | Francis Galton, Karl Pearson, W.F.R. Weldon |
| **Traits studied** | Discrete (e.g., pea colour, seed shape) | Continuous (e.g., height), normally distributed |
| **Key evidence** | Particulate inheritance; segregation ratios (3:1) | Correlation and regression between relatives |
| **View of evolution** | Discontinuous jumps (mutation) | Gradual, continuous change (Darwin) |

The biometricians measured traits in families and found that relatives resemble each other in proportion
to how closely they are related, with offspring regressing toward the population mean. They doubted that
discrete Mendelian factors could produce these smooth patterns.

**R.A. Fisher resolved the debate in 1918.** His paper *The correlation between relatives on the
supposition of Mendelian inheritance* showed that:

- **many Mendelian loci**, each with a small additive effect, plus environment, produce a continuous,
  normally distributed trait;
- the model predicts the **correlations between relatives** that the biometricians had measured; and
- the trait's **variance** (a word Fisher introduced) can be partitioned into additive, dominance and
  epistatic parts, the foundation of {doc}`heritability`.

The Royal Society declined to publish the paper (Pearson was one of the referees), so it appeared in
the *Transactions of the Royal Society of Edinburgh*.

## The Galton board: binary events → normal curve

Galton, the founder of the biometric school, built a device that captures Fisher's idea.

```{figure} images/galton_board.png
:alt: Galton board. A ball bounces left or right at each of 12 rows of pegs. The bins at the bottom fill up into a bell-shaped histogram that matches a normal curve.
:width: 420px

A Galton board with 12 rows of pegs. Each ball's final bin is the number of times it bounced right.
```

▶ **Watch one in action:** [a Galton board on YouTube](https://www.youtube.com/shorts/TwctT3Ncm1w) (short video).

- At each peg a ball bounces left or right: a **50:50 binary event**.
- Its final bin is the number of rightward bounces, the **sum of many binary outcomes**, so the bin counts
  follow a binomial distribution.
- With many rows, the bins trace a **normal (bell) curve**.

**The genetic analogy.** Treat each peg as one allele a person inherits, either a risk allele (bounce
right) or a non-risk allele (bounce left). The bin is the person's total number of risk alleles. One
peg gives two outcomes, which is discrete, Mendelian variation; many pegs give a bell curve, which is
continuous, biometric variation. Galton described the board in *Natural Inheritance* (1889) to
illustrate the normal law. It is a nice irony that it also illustrates how Mendelian factors produce
the continuous variation his followers studied.

## From binomial loci to a normal trait

Now make the analogy exact. Suppose a trait is influenced by $M$ biallelic loci, each with risk-allele
frequency $p$, and every risk allele adds the same small amount. A person's **genetic score** is the
number of risk alleles they carry across all loci. Each person has two alleles at each locus, so the score is
binomial with $2M$ trials:

$$
\text{score} \sim \text{Binomial}(2M,\ p), \qquad \text{mean} = 2Mp, \qquad \text{variance} = 2Mp(1-p)
$$

```{figure} images/binomial_panels.png
:alt: Distribution of the number of risk alleles carried for 1, 2, 5, 20 and 100 loci. One locus gives three genotype classes in 1 to 2 to 1 proportions; with more loci the distribution approaches a normal curve.
:width: 100%

Exact binomial distributions of the genetic score for 1, 2, 5, 20 and 100 loci (risk-allele frequency
0.5, equal effect per allele). The purple line is the normal curve with the same mean and variance.
```

- **One locus:** a person carries 0, 1 or 2 risk alleles (genotypes aa, Aa, AA) in 1 : 2 : 1 proportions
  under Hardy–Weinberg equilibrium. Three discrete classes: Mendelian variation.
- **Five loci:** the bars already follow the normal curve closely.
- **20 to 100 loci:** the distribution is indistinguishable from a normal curve.

This is the **central limit theorem**: the sum of many small, independent effects is approximately
normally distributed. It holds when the effects are

- **many, each small**, with no single locus of large effect;
- **independent**, for example loci that are not linked; and
- **additive**, so that effects simply add up.

Environmental influences, themselves the sum of many small factors, smooth out the remaining steps, so
the phenotype is continuous. A single large-effect variant breaks the first condition: a pathogenic
*BRCA1* or familial hypercholesterolemia variant can make a trait distribution skewed or bimodal, which
connects back to {doc}`variation-spectrum`.

The genetic score here is essentially what a **polygenic risk score** estimates, with each allele
weighted by its effect size from a GWAS ({doc}`heritability-to-gwas`, {doc}`../polygenic-risk-scores/index`).

## Common disease: the liability-threshold model

Height is continuous, but most diseases are present or absent. How can a normal model apply to a binary
trait? The **liability-threshold model**, formalized by Falconer in 1965, proposes that:

- every person has an unobserved **liability** to the disease, the sum of many genetic and environmental
  contributions, so by the central limit theorem it is normally distributed;
- people whose liability is above a **threshold** are affected, so the area of the tail above the
  threshold equals the disease prevalence; and
- relatives of an affected person share many of their risk alleles, so their liability distribution is
  **shifted to the right**. The threshold stays where it is, so more relatives exceed it.

```{figure} images/liability_threshold.png
:alt: Normal distribution of liability in the general population with a threshold; the small tail above the threshold is affected. The distribution for relatives of affected individuals is shifted right, so a larger area lies above the threshold.
:width: 640px

The liability-threshold model, illustrated for a disease with 5% prevalence. Relatives of affected
individuals (dashed) have a shifted liability distribution, so more of them are above the threshold.
```

This explains a long-standing puzzle: common diseases **cluster in families without following
Mendelian segregation ratios**. It also means that a polygenic risk score is an estimate of where a
person sits on the genetic part of the liability scale. That is a probability, not a diagnosis.

**How much does family history raise risk?** Under the model, the mean liability of affected people is
$i = \varphi(T)/K$ standard deviations above the population mean, where $K$ is the prevalence, $T$ the
threshold and $\varphi$ the normal density. A first-degree relative shares half their genome with an
affected person, so their expected liability is shifted by about $\tfrac{1}{2}\,h^2 \, i$, where $h^2$ is
the narrow-sense heritability of liability (see {doc}`heritability`). The cell below uses this
approximation, which keeps the relatives' variance at 1.

```{code-cell} ipython3
from scipy.stats import norm

K = 0.05                # prevalence (5%)
T = norm.isf(K)         # threshold on the liability scale
i = norm.pdf(T) / K     # mean liability of affected people

for h2 in (0.3, 0.6, 0.8):
    shift = 0.5 * h2 * i              # relatives share half their genome
    risk = norm.sf(T - shift)         # relatives above the threshold
    print(f"h2 = {h2}: relatives' risk = {risk:.1%} ({risk / K:.1f}x)")
```

With a heritability of 0.6, first-degree relatives of cases have about three times the population risk,
even though no single gene is responsible.

## Self-check

Pick an answer to see whether it is right and why. You can change your answer as often as you like.

```{code-cell} ipython3
:tags: [remove-input]
import sys; sys.path.insert(0, "..")
from gab_quiz import quiz
quiz("quizzes/polygenic-model.json")
```
