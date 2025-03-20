---
layout: default
title: 6. Parental haplotype imputation
nav_order: 6
parent: Parent-of-origin inference pipeline tutorial
has_children: true

---
# haploid imputation of parental haplotypes

{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---



## Rationale
In this section, we will perform haploid imputation of our inter-chromosomally phased data, which separately impute each parental haplotype.

---


## Pipeline

To perform genotype imputation, we will use the software [IMPUTE5](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1009049). For any question/information on the software, please refer to the [official documentation](https://jmarchini.org/software/#impute-5).
Please also note that IMPUTE5 is freely available for academic use only.
You can download a binary of the software on [the official download link](https://www.dropbox.com/scl/fo/ukwimchnvp3utikrc3hdo/AKqYvE6-9C5kLpKDSfhR8xQ?rlkey=n2zty39bdst5j5tycd0sf89ee&e=1&dl=0).

All scripts for this step are in folder `pipeline/step4_impute_haplotypes/`.

---

### 0. Reference panel
The genotype imputation process impute missing variant from a reference panel of haplotypes. For this tutorial, since the aim is to rely solely on publicly available data, we will use the 1,000 genome project WGS files as reference panel. Note, however, that use the same individual as input AND reference panel is not correct at all. When imputing you data, please use traditional reference panel, such as the Haplotype Reference Consortium (HRC).

To download the example reference panel, you can run:

<div class="code-example" markdown="1">
```bash
bash step0_download_fake_ref_panel.sh
```
</div>


---


### 1. Prepare imputation input files
Here, we will first remove from our inter-chromosomally phased data all individual who also have a parent available (i.e our validation cohort). Then, we will add all individuals with available parental genomes, for which we performed pedigree-based phasing.

To do this, you can simply run:

<div class="code-example" markdown="1">
```bash 
bash step1_prepare_input_data.sh
```
</div>

---


### 2. Haploid genotype imputation
In this last step, we will perform haploid genotype imputation of our inter-chromosomally phased data. For this tutorial, we impute only a small chunk of the chromosome 20. When running on your own data, please adapt the genomic region you wish to impute.

To perform genotype imputation, run the command:

<div class="code-example" markdown="1">
```bash 
bash step2_imputation.sh
```
</div>










