---
layout: default
title: 4. Parental side inference
nav_order: 4
parent: Parent-of-origin inference pipeline tutorial
has_children: true

---
# Parental side inference from chromosome X, mtDNA and cross-over inference in siblings

{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---



## Rationale

Inferring close relative groups and using them as surrogate parent allows us to perform inter-chromosomal phasing, effectively segregating the focal individual's haplotypes into two sets of 22 haplotypes, one inherited from each parent. However, we need additional inference to determine which haplotype set is paternally inherited and which one is maternally inherited.

To do so, we use two key approaches:

- First, for second- to fourth-degree relatives, the PofO was determined by identifying whether the surrogate parent was on the maternal or paternal side. Once the PofO of a single haplotype segment was established through its IBD status with the used surrogate parents, this assignment could be confidently extended to the entire parental haplotype set. To determine whether the used surrogate parents are on the maternal or on the paternal side, we leveraged chromosome X genotype data and mtDNA whole-genome sequence data, while also demonstrating the potential value of including chromosome Y whole-genome sequences. Chromosome X is maternally inherited in males, and mtDNA is maternally inherited in both sexes, allowing us to identify maternal relationships. Although chromosome Y data can reliably identify paternal relationships, we found this approach to be less cost-effective and therefore chose not to implement it broadly. By analyzing IBD sharing with surrogate parents on these specialized chromosomes, we were able to determine whether the surrogate parent (and, consequently, all IBD segments shared on other chromosomes) belonged to the maternal or paternal side.

- For siblings, who cannot be assigned exclusively to a single parental side since they inherit genomes from both parents, we applied a different strategy. Specifically, we inferred crossover positions in siblings and probabilistically assigned them a PofO using sex-specific recombination maps. For this approach, we prioritized inter-chromosomally phased data that improved accuracy by aggregating crossovers events across chromosomes.


