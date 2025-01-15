---
layout: default
title: 1. Chromosome X
nav_order: 1
parent: 3. Parental side inference
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

---


## IBD mapping on chromosome X


---

**1. Validation cohort**

We first use THORIN to map IBD on chromosome X for individuals in our validation cohort. A first look at the difference in IBD segment length between maternal and paternal relative indicates us whether this information can be use on the remaining individuals (i.e those not included in the validation cohort).

To map IBD on chromosome X using THORIN for individuals in the validation cohort, you can use the script `pipeline/step1.2_parental_side_determination/step0_chrX_ibd/step0\_map\_IBD\_benchmark\_set.sh`, also included below:


<div class="code-example" markdown="1">
```bash
CHR=X

IN=../../step0_download_genotype/data/genotype/vcf/KGP.chr${CHR}.gsa.diploidized.bcf
UNR=../../step1_surrogate_parents/data/relatedness/unrelated/KGP.chr${CHR}.gsa.unrelated.bcf
MAP=../../../maps/chrX.b38.gmap.gz
GRP=../../step1_surrogate_parents/data/benchmark/Relatives.benchmark.group

BIN=../../../bin/thorin_v1.2_static

ODIR=data/THORIN/benchmark
mkdir -p ${ODIR}
OUT=${ODIR}/KGP.chr${CHR}.benchmark.thorin.prob

./${BIN} -I ${IN} -H ${UNR} -M ${MAP} -R chr${CHR} -G ${GRP} -O ${OUT} --ibd ${OUT}.ibd
```
</div>


You can then simply plot IBD sharing with paternal and maternal relatives to have a first overview of the IBD sharing differences between relative groups. The code provided in `pipeline/step1.2_parental_side_determination/step0_chrX_ibd/plot_chrX_IBD_validation_cohort.R` should give you the following plots:


![](https://github.com/RJHFMSTR/THORIN/blob/main/pipeline/step1.2_parental_side_determination/step0_chrX_ibd/chrX_IBD_validation_cohort.png?raw=true)

---


**2. Entire data**
Now that we verified that paternal and maternal surrogate parent have significant differences in IBD sharing on chromosome X in our cohort, we can extend this approach to the entire cohort. 

For this, we will use the two scripts `pipeline/step1.2_parental_side_determination/step0_chrX_ibd/step1_map_IBD_all.sh` and `pipeline/step1.2_parental_side_determination/step0_chrX_ibd/step2_assign_parental_side_and_plot.R`.

In the first script, `pipeline/step1.2_parental_side_determination/step0_chrX_ibd/step1_map_IBD_all.sh`, we use THORIN to map IBD on chromosome X for the entire cohort with groups of surrogate parents.

In the second script, `pipeline/step1.2_parental_side_determination/step0_chrX_ibd/step2_assign_parental_side_and_plot.R`, we use the validation cohort to derive a probability of a relative being on the paternal or maternal side given it's chromosome X IBD sharing with the focal individual. We compute the accuracy of this process and we probabilistically assign each relative group to a parental side.

After running those two script on the KGP example data, you should have the following plots. This is the equivalent of our supplementary figure on chrX in our [manuscript](https://www.medrxiv.org/content/10.1101/2024.12.03.24318392v1). 





![](https://github.com/RJHFMSTR/THORIN/blob/main/pipeline/step1.2_parental_side_determination/step0_chrX_ibd/chrx_accuracy_derived_prediction?raw=true)





---



