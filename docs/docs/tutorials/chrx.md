---
layout: default
title: Chromosome X
nav_order: 1
parent: Parental side inference
has_children: true

---
# Parental side inference from chromosome X

{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---



## Rationale

In this step, we aim to assign the parental side of relatives using chromosome X. To so do, we use the simple idea that male individuals inherited their chromosome X from the maternal side of the family. Hence, if we can find IBD on the chromosome X between a male target and a given relative, this relative would be on the maternal side of the family.

In a first step, we use the validation cohort (i.e, individuals with both available parent and close relatives, for which the parental side of the relatives can be determined by looking at it kinship estimate with the available parents). 

In a second step, we use the validation cohort to derive a probability of being on the maternal or paternal side from the length of the segment shared IBD with the target individual.

All scripts are available in folder `pipeline/step1.2_parental_side_determination/step0_chrX_ibd`.


