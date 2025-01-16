---
layout: default
title: 3. Inter-chromosomal phasing
nav_order: 3
parent: Parent-of-origin inference pipeline tutorial
has_children: true

---
# Inter-chromosomal phasing tutorial
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---



## Rationale
We previusly used kinship estimates, age and sex to identify parent-offspring trios, duos and siblings. For individuals with more distant relatives, we utilized a clustering approach to segregate relatives by parental sides, thereby classifying them into "one family side" versus "the other family side".

Here, we will use identity-by-descent (IBD) information from surrogate parents to conduct inter-chromosomal phasing, effectively identifying sets of variants inherited together from the same parent across different chromosomes. By analyzing haplotype segments shared IBD with the same set of surrogate parents across the 22 autosomes for a given target individual, we can determine which haplotype segments are inherited from the same parent. This enables the construction of partial parental haplotype sets, which can then be used to perform inter-chromosomal phasing—segregating haplotypes inherited from each parent across the genome. This approach goes beyond traditional phasing methods, which are limited to resolving haplotypes within individual chromosomes (intra-chromosomal phasing).

In practice, we use THORIN to map IBD with surrogate parent groups and identify IBD segments, and we then filter out IBD segments smaller than 3 centimorgans (cM). Larger IBD segment are automatically included into a scaffold file. This scaffold was then used as input for the [SHAPEIT5 phase\_common tool](https://odelaneau.github.io/shapeit5/docs/documentation/phase_common/). This step allowed us to refine and re-estimate haplotypes from genotype data, while simultaneously correcting intra-chromosomal phasing switch errors and performing inter-chromosomal phasing by assigning all haplotypes shared with the same surrogate parents to the same parental haplotype (e.g., first or second) across all 22 autosomes.

Scripts for this step are in `pipeline/step2_interchromosomal_phasing/`.

---

## Pipeline


### IBD-based scaffold construction
In this step, we use THORIN to map IBD with surrogate parent groups and to construct the IBD-based haplotype scaffold. The core script for this test is included below and located in `pipeline/step2_interchromosomal_phasing/src/thorin_scaffold.sh`. To directly run this script in parallel on the 22 autosomes, use the bash script `pipeline/step2_interchromosomal_phasing/step0_build_scaffold_from_IBD.sh` and adapt the number of threads and the software path.


<div class="code-example" markdown="1">
```bash
CHR=20
threads=1

IN=../step0_download_genotype/data/genotype/vcf/KGP.chr${CHR}.gsa.vcf.gz
UNR=../step1_surrogate_parents/data/relatedness/unrelated/KGP.chr${CHR}.gsa.unrelated.bcf
MAP=/home/rhofmeis/Dropbox/Ressources/Git_repository/THORIN/maps/chr${CHR}.b38.gmap.gz
GRP=../step1_surrogate_parents/data/Relatives.group

BIN=THORIN/bin/thorin_v1.2_static

ODIR=data/THORIN
mkdir -p ${ODIR}
OUT=${ODIR}/KGP.chr${CHR}.thorin.prob

if [ ! -f "${OUT}.ibd.scaffold.bcf.csi" ]; then
	${BIN} -I ${IN} -H ${UNR} -M ${MAP} -R chr${CHR} -G ${GRP} -O ${OUT} --ibd ${OUT}.ibd --scaffold ${OUT}.ibd.scaffold.bcf -T ${threads}
	bcftools index ${OUT}.ibd.scaffold.bcf --threads ${threads}
fi
```
</div>

When running this script, you will have one output file in VCF/BCF format for each of the autosome. This file is the IBD-based haplotype scaffold that we will use in the next step to perform inter-chromosomal phasing.


If you want to have a look at your data and how the IBD mapping and scaffold construction work, we also provide an R script to plot your data (`pipeline/step2_interchromosomal_phasing/src/plot_ibd_tracks.R`). In the figure below, shown for the individual `HG00544` of the 1,000 Genome Project example data, the top and bottom plots shows the probability of sharing IBD between the focal individual and its surrogate parent group(s) or unrelated individuals (Unr) on haplotype 0 (top) or haplotype (1). The middle plot shows the boundaries of IBD segment identified, with their respective labels (A=H0 shared with G1 and/or H1 with G2;BA=H1 shared with G1 and/or H0 with G2; C= both H0 and H1 shared with Unr; D=both H0 and H1 shared with G1 or G2). For the scaffold construction, we used only labels A or B, with segment longer than 3cM. The size of the segment used is indicated below their respective labels. Unused IBD segment have no size indicated.

![](https://github.com/RJHFMSTR/THORIN/blob/main/pipeline/step2_interchromosomal_phasing/Plots/thorin_IBD_plot.HG00544.png?raw=true)




### Inter-chromosomal phasing from IBD-based scaffolds
In this step, we will use the IBD-based haplotype scaffold as input in the phasing software SHAPEIT5 to perform inter-chromosomal phasing. Static binairies for the SHAPEIT5 software can be found [here](https://github.com/odelaneau/shapeit5/releases).

The core script for this test is included below, and located in `pipeline/step2_interchromosomal_phasing/src/shapeit5_interchromosomal_phasing.sh`. To directly run this script in parallel on the 22 autosomes, use the bash script `pipeline/step2_interchromosomal_phasing/step1_interchrs_phasing.sh` and adapt the number of threads and the software path.

<div class="code-example" markdown="1">
```bash
CHR=20
threads=1

IN=../step0_download_genotype/data/genotype/vcf/KGP.chr${CHR}.gsa.vcf.gz
MAP=/home/rhofmeis/Dropbox/Ressources/Git_repository/THORIN/maps/chr${CHR}.b38.gmap.gz
SCAF=data/THORIN/KGP.chr${CHR}.thorin.prob.ibd.scaffold.bcf

BIN=shapeit5/phase_common/bin/phase_common_static

ODIR=data/SHAPEIT5
mkdir -p ${ODIR}
OUT=${ODIR}/KGP.chr${CHR}.thorin_shapeit5.vcf.gz

if [ ! -f "${OUT}.csi" ]; then
        ${BIN} -I ${IN} -M ${MAP} -O ${OUT} -R ${CHR} -T ${threads} -S ${SCAF}
        bcftools index ${OUT} --threads ${threads}
fi
```
</div>


When running this script, you will have one output file in VCF/BCF format for each of the autosome. These files are inter-chromosomally phased data.

### Inter-chromosomally phased data filtering
Since we included in the scaffold file only haplotype segment longer than 3cM shared IBD with surrogate parents, it is possible that some chromosome of a given individual not sharing IBD with surrogate parents were not scaffolded.

To know which chromosomes of each individuals were scaffolded and can be used in downstream analysis, you need to use the `.ibd` output file of THORIN and keep only segments in class A or class B longer than 3cM. For simplycity, the code below was directly included within the THORIN script `pipeline/step2_interchromosomal_phasing/src/thorin_scaffold.sh`. As an example, on chromosome 20, you should have a total of 215 individuals when using the 1,000GP example data.

<div class="code-example" markdown="1">
```bash
CHR=20
ODIR=data/THORIN
OUT=${ODIR}/KGP.chr${CHR}.thorin.prob
SAMP=${OUT}.ibd.scaffolded_samples.txt
awk -F'\t' 'NR > 1 && $5 >= 3 && ($4 == "A" || $4 == "B") { print $6 }' "${OUT}.ibd" | sort | uniq > ${SAMP}
```
</div>















