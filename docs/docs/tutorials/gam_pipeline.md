---
layout: default
title: Genetic Assortative Mating pipeline tutorial
nav_order: 3
parent: Tutorials
---
# Genetic Assortative Mating pipeline tutorial
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }


---



## Rationale

Assortative mating, that is, the tendency for individuals to select partners with similar characteristics, leaves detectable imprints in an individual's genome: when parents are genetically similar for a trait, the alleles they transmit tend to align more than expected under random mating, resulting in similar distributions of trait-associated alleles on the maternal and paternal haplotypes. Consequently, PGS computed from the two parental haplotypes of an individual are expected to correlate. Exploiting this phenomenon, we can estimate GAM by correlating maternal and paternal haplotype-based PGS across individuals. 

To estimate genetic assortative mating (GAM) from biobank individuals, in the absence of mate-pairs, we first reconstructed and separated maternally and paternally inherited haplotypes. For this, we followed [our pipeline](https://rjhfmstr.github.io/THORIN/docs/tutorials/interchrs.html) to perform inter-chromosomal phasing from inferred surrogate parents.

We then compute a PGS for each haplotype of each individual, effectively proxying partial parental PGS, which we correlate to estimate GAM.

In this tutorial, we provide a detailed explanation on how to estimate GAM from:
- A) mate-pairs
- B) biobank individuals, using inferred parental haplotypes (i.e inter-chromosomally phased data).

All scripts specific to this projects are hosted on a dedicated github: [https://github.com/RJHFMSTR/GAM](https://github.com/RJHFMSTR/GAM)

---




