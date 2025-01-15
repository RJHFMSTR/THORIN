---
layout: default
title: 1. Example data
nav_order: 1
parent: Parent-of-origin inference pipeline tutorial
---
# Example data preparation
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---



## Rationale
Since large-scale population-based biobanks are available under restricted access only, we based our tutorials on publicly available data from the 1,000 Genome project. In the following steps, we will download the genotype data and pre-process it to simulate SNP-array data typically available in large-scale population-based biobanks.

All these steps are directly available on our github in folder [pipeline/step0\_download\_genotype/](https://github.com/RJHFMSTR/THORIN/tree/main/pipeline/step0_download_genotype).

---

## Requirements

In addition of using C++ software such as THORIN, SHAPEIT5 and IMPUTE5, our pipeline also uses standard tools for genomic analysis and R packages that you may need to install by yourself. For this part, we will need `bcftools`, `plink1.9` and `plink2`.

<div class="code-example" markdown="1">
```bash
sudo apt install -y bcftools pink1.9 plink2
```
</div>


## Data download

### Environment settings

<div class="code-example" markdown="1">
```bash
threads=16
ODIR=data
mkdir -p ${ODIR}
```
</div>


---

### Infinium Global Screening Array sites, to simulate SNP-array sites from whole-genome sequencing data

<div class="code-example" markdown="1">
```bash
wget https://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/global-screening-array-24/infinium-global-screening-array-24-v1-0-c2-annotated.zip -O ${ODIR}/gsa.zip
unzip data/gsa.zip -d ${ODIR}/ && rm ${ODIR}/gsa.zip
awk 'NR > 1 {print "chr"$2"\t"$3}' ${ODIR}/GSA-24v1-0_C2.hg38.annotated.txt > ${ODIR}/GSA.variant_positions.txt
rm ${ODIR}/GSA-24v1-0_C2.hg38.annotated.txt
```
</div>


---

### Download genotype files

<div class="code-example" markdown="1">
```bash
DIR="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"
PRE="CCDG_14151_B01_GRM_WGS_2020-08-05"
SUF="filtered.shapeit2-duohmm-phased.vcf.gz"
ODIR=data/genotype
mkdir -p ${ODIR}/vcf

for CHR in {1..22} X; do
	
	if [ ${CHR} == "X" ]; then
	        SUF="filtered.eagle2-phased.v2.vcf.gz"
	fi



	if [ ! -f "${ODIR}/vcf/${PRE}_chr${CHR}.gsa.${SUF}.csi" ]; then
	
	        # download chromosome-wide vcf.gz file
	        wget ${DIR}/${PRE}_chr${CHR}.${SUF} -O ${ODIR}/vcf/${PRE}_chr${CHR}.${SUF}
	        wget ${DIR}/${PRE}_chr${CHR}.${SUF}.tbi -O ${ODIR}/vcf/${PRE}_chr${CHR}.${SUF}.tbi
	
	        # subset to GSA sites
	        bcftools view -T data/GSA.variant_positions.txt ${ODIR}/vcf/${PRE}_chr${CHR}.${SUF} -Oz -o ${ODIR}/vcf/${PRE}_chr${CHR}.gsa.${SUF} --threads ${threads}
	        bcftools index -f ${ODIR}/vcf/${PRE}_chr${CHR}.gsa.${SUF} --threads ${threads}
	        rm ${ODIR}/vcf/${PRE}_chr${CHR}.${SUF}*

	fi



	if [ ! -f "${ODIR}/plink/KGP.chr${CHR}.gsa.bed" ]; then
	
	        # convert to plink file (to use with the relatedness inference software KING later)
	        mkdir -p ${ODIR}/plink
	        plink2 --vcf ${ODIR}/vcf/${PRE}_chr${CHR}.gsa.${SUF} --make-bed --out ${ODIR}/plink/KGP.chr${CHR}.gsa

	        mv ${ODIR}/vcf/${PRE}_chr${CHR}.gsa.${SUF} ${ODIR}/vcf/KGP.chr${CHR}.gsa.vcf.gz
	        mv ${ODIR}/vcf/${PRE}_chr${CHR}.gsa.${SUF}.csi ${ODIR}/vcf/KGP.chr${CHR}.gsa.vcf.gz.csi
	fi
done
```
</div>

---

### Merge chromosome level files into a genome-wide file

<div class="code-example" markdown="1">
```bash
printf "%s\n" ${ODIR}/genotype/plink/KGP.chr{1..22}.gsa > merge_list.txt
plink1.9 --merge-list merge_list.txt --make-bed --out data/genotype/plink/KGP.merged_chromosomes 
rm ${ODIR}/genotype/plink/KGP.chr*.gsa.*
rm merge_list.txt
```
</div>

---

### Re-format .fam files

We modify the fam file so that they don't all have the same family ID == 0, and set FID=IID.

<div class="code-example" markdown="1">
```bash
awk '{print $2" "$2" "$3" "$4" "$5" "$6}' data/genotype/plink/KGP.merged_chromosomes.fam > data/genotype/plink/KGP.merged_chromosomes.v2.fam
mv data/genotype/plink/KGP.merged_chromosomes.v2.fam data/genotype/plink/KGP.merged_chromosomes.fam
```
</div>


---

### Download sample information
For the clustering of close relatives, we will need age and sex of participants. Here, the age is not available for the 1,000 GP participants.

<div class="code-example" markdown="1">
```bash
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt -O data/1kGP.3202_samples.pedigree_info.txt
```
</div>

---

### Chromosome X specific processing
We apply additional Quality-Control steps to the chromosome X. This is due to the fact that because of for example genotyping errors, males can have heterozygous sites. We therefore force homozygosity for males.

<div class="code-example" markdown="1">
```bash

BIN=diploidize_static # available from https://github.com/odelaneau/otools
SEX=data/1kGP.3202_samples.pedigree_info.txt
VCF=data/genotype/vcf/KGP.chrX.gsa.vcf.gz
threads=8

## Get male individuals
awk '$4 == 2 {print $1}' data/1kGP.3202_samples.pedigree_info.txt > data/KGP.males.txt

## QC male individuals --> transform mixed diploid and haploid call into fully diploid calls.

IN=data/genotype/vcf/KGP.chrX.gsa.vcf.gz
OUT_M=data/genotype/vcf/KGP.chrX.gsa.males.bcf
bcftools view -S data/KGP.males.txt -Ob -o ${OUT_M} ${IN} --threads ${threads}
bcftools index ${OUT_M} --threads ${threads}


OUT_F=data/genotype/vcf/KGP.chrX.gsa.females.bcf
bcftools view -S ^data/KGP.males.txt -Ob -o ${OUT_F} ${IN} --threads ${threads}
bcftools index ${OUT_F} --threads ${threads}

OUT_D=data/genotype/vcf/KGP.chrX.gsa.males.diploidized.bcf
./${BIN} --input ${OUT_M} --output ${OUT_D}
bcftools index ${OUT_D} --threads ${threads}

OUT=data/genotype/vcf/KGP.chrX.gsa.diploidized.bcf
bcftools merge ${OUT_D} ${OUT_F} -Ob -o ${OUT}
bcftools index ${OUT} --threads ${threads}

rm ${OUT_D}* && rm ${OUT_M}* && rm ${OUT_F}*

```
</div>










