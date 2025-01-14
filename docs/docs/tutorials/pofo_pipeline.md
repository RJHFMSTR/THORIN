---
layout: default
title: Parent-of-origin inference pipeline tutorial
nav_order: 2
parent: Tutorials
---
# Parent-of-origin inference pipeline tutorial
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---



## Rationale
THORIN is a tool to map IBD between a focal individual and reference individuals. We developped this tool to map IBD between a focal individual and its surrogate parents, allowing us to infer the parent-of-origin of alleles of this focal individual.

We detail here the entire pipeline for inferring the parent-of-origin of alleles a focal individuals. In addition of using THORIN, this pipeline also uses several state-of-the-art software for phasing and imputation, such as SHAPEIT5 and IMPUTE5.

This pipeline was initially developped on the UK Biobank cohort, and then also applied to the Estonian Biobank cohort. However, these biobanks data are under restricted access only.

For the purpose of this tutorial, we adapted our scripts to fit the 1,000 Genome Project, which provides publicly available genotype data. However, given the small number of individuals in this cohort compared to large-scale population-based biobanks, we can not expect to have a similar accuracy as on the UK Biobank and Estonian Biobank.

---




