---
layout: default
title: 2. Cluster relatives
nav_order: 2
parent: Parent-of-origin inference pipeline tutorial
---
# Relatedness inference and clustering of close relatives
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---



## Rationale
In this step, we will identify parent-offspring duos and trios, siblings, and more distant relatives.

In addition, we will cluster more distant relatives into groups, that represent the group that we use as surrogate parents when inferring the parent-of-origin of haplotypes.

All scripts of this step are available in folder `pipeline/step1_surrogate_parents/`.

## Pipeline

### Infer relatedness

Here, we use the 1,000 Genome Project plink genome-wide files that we prepared [earlier](https://rjhfmstr.github.io/THORIN/docs/tutorials/data_download.html) to infer global relatedness estimates. For this, we use the [KING software](https://www.kingrelatedness.com/). 

<div class="code-example" markdown="1">
```bash
wget https://www.kingrelatedness.com/Linux-king.tar.gz
tar -xvzf Linux-king.tar.gz
chmod +x king

BED=../step0_download_genotype/data/genotype/plink/KGP.merged_chromosomes.bed
PFX=data/relatedness/KGP.king_relatedness
mkdir -p data/relatedness

./king -b ${BED} --related --degree 5 --cpus 4 --prefix ${PFX}
```
</div>



---


### Cluster relatives
In this step, we identify parent-offspring duos and trios, siblings and more distant relatives. We also cluster more distant relatives into surrogate parent groups. Finally, we identify individual for the validation cohort, that is, individuals with both available parental genomes and surrogate parents groups. We use those individuals to assess accuracy of the method and derive specific metrics for sex-chromosome parental side integration.


You can find the main script for this step is `pipeline/step1\_surrogate_parents/src/grouping.R'.
This script relies on function that you can find in script `pipeline/step1\_surrogate\_parents/src/grouping\_functions.KGP.R`. Please note that this script has been adapted to fit the available 1,000GP data, for which we don't have the age of participants. To run on you own cohort, please use script `pipeline/step1\_surrogate\_parents/src/grouping\_functions.R` instead.


After this step, we should have data in folder `pipeline/step1\_surrogate\_parents/data/`:
- Trios.ped : parent-offspring trios. Format : offspring father mother. This file can directly be used as input for SHAPEIT5 to perform inter-chromosomal phasing and parent-of-origin inference from parental genomes. See SHAPEIT% [documentation](https://odelaneau.github.io/shapeit5/docs/documentation/phase_common/#usage2-phasing-related-samples).
- Duos.ped : : parent-offspring trios. Format : offspring father mother. Missing parent is indicated with NA. This file can directly be used as input for SHAPEIT5 to perform inter-chromosomal phasing and parent-of-origin inference from parental genomes. See SHAPEIT% [documentation](https://odelaneau.github.io/shapeit5/docs/documentation/phase_common/#usage2-phasing-related-samples).
- Sibs.list : list of sibling pairs. Format: sibling1 sibling2.
- Relatives.group : Cluster of close relatives. Format: target group1 group2. Missing group is empty.
- Relatives.male_targets.group : Same file as `Relatives.group` restricted to male individuals. This is use to assign parental side for male individuals using chromosome X (see next steps).
- benchmark/Relatives.benchmark.group : Same file as `Relatives.group` restricted to individual having both close relative clusters and parental genomes.
- benchmark/Relatives.benchmark.side : Assignment of close relative to parental side using the direct relatedness between the parent and the close relative. Format : target relative group side target_sex relative_sex degree


---


### Unrelated individuals

For THORIN, we need a set of individual fully unrelated to all our focal individual. To build this set of individual, use the following code (also available in script `pipeline/step1\_surrogate\_parents/step2\_get\_100\_unrelated.sh``.


<div class="code-example" markdown="1">
```bash
mkdir -p data/relatedness/unrelated
threads=8

# list related individuals.
REL=data/relatedness/KGP.king_relatedness.kin0
awk '$14 != "UN" {print $1 "\n" $3}' ${REL} | grep -v FID | sort | uniq > data/relatedness/related.samples

# list all samples
IN=../step0_download_genotype/data/genotype/vcf/KGP.chr22.gsa.vcf.gz
bcftools query -l ${IN} > data/all.samples

# get 100 random unrelated
UNR=data/relatedness/unrelated/random_unrelated_samples.txt
comm -23 <(sort data/all.samples) <(sort data/relatedness/related.samples) | shuf -n 100 > ${UNR}

# subset each genotype file to keep only unrelated samples
for CHR in {1..22}; do
	IN=../step0_download_genotype/data/genotype/vcf/KGP.chr${CHR}.gsa.vcf.gz
	OUT=data/relatedness/unrelated/KGP.chr${CHR}.gsa.unrelated.bcf
	bcftools view -S ${UNR} -Ob -o ${OUT} ${IN} --threads ${threads}
	bcftools index ${OUT} --threads ${threads}
done

for CHR in X; do
	IN=../step0_download_genotype/data/genotype/vcf/KGP.chr${CHR}.gsa.diploidized.bcf
	OUT=data/relatedness/unrelated/KGP.chr${CHR}.gsa.unrelated.bcf
	bcftools view -S ${UNR} -Ob -o ${OUT} ${IN} --threads ${threads}
	bcftools index ${OUT} --threads ${threads}
done
```
</div>




