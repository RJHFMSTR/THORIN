---
layout: default
title: 3. Mitochondrial DNA
nav_order: 3
parent: 4. Parental side inference
has_children: true

---
# Parental side inference from mitochondrial DNA

{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---



## Rationale

In this step, we assign the parental side of relatives using the mitochondrial DNA, which is exclusively inherited from the mother, and serves as a valuable tool for identifying maternal relatives.

However, typical IBD mapping software are derived from the Li \& Stephens Hidden Markov Model, which excels at modeling the human recombinant genomes architecture by accounting for mutation and recombination between haplotypes. Due to mtDNA's lack of recombination and its higher mutation rate compared to autosomes, these software proved unsuitable. As a result, we adopted an alternative approach to evaluate non-recombinant DNA sharing between pairs of individuals. This approach is inspired from the Jaccard index and termed here Minor Variant Sharing (MVS). It involves evaluating the proportion of shared minor alleles between relative pairs. This allows us to overcome the limitations of traditional IBD mapping software for mtDNA by basing our analysis on IBS (Identity-By-State) rather than IBD.

In addition, the UK Biobank Axiom array data proved inadequate for this purpose due to the limited availability of genotyped variants (N=265). To overcome this limitation, we used the whole-genome sequencing GraphTyper cram file available on the UK Biobank Research and Analysis Platform (RAP) for ~500,000 individuals. We called variants from the mtDNA whole-genome sequencing (WGS) cram files for 274,525 UK Biobank individuals (those with inter-chromosomal phasing) and their surrogate parents using the MitoHPC software.

All scripts for this step are available in folder `pipeline/step3_parental_side_determination/step2_mtDNA`.


Note: For the example on the 1,000 Genomes Project data, we will start with the publicly available VCF files to avoid re-mapping variant using MitoHPC.



## Pipeline

### 1. Download data
Instead of mapping the mtDNA variant using MitoHPC, we'll start - for this example - with the publicly available VCF files. You can download the files using the above code, also included in the script `pipeline/step3_parental_side_determination/step2_mtDNA/step0_data_download.sh`.



<div class="code-example" markdown="1">
```bash
mkdir -p data/
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz
mv ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz data/
```
</div>



---


### 2. Compute mtDNA Minor Variant Sharing (MVS)
Here, we use a Jaccard index to assess the relatedness of a pair of individual by the number of variant they share. To run this part, use the script `pipeline/step3_parental_side_determination/step2_mtDNA/src/mtdna_mvs.R` and adapt the input VCF data if not using the KGP example data.

This script will produce the output data , with the following columns:
- target : ID of the focal individual
- relative: ID of the relative
- group : group of the relative (used when mapping IBD and performing inter-chromosomal phasing)
- degree : relatedness degree between target and relative
- side : parental side of the relative, if any (available only for individual included in the validation cohort)
- mtdna : mitochondrial DNA MVS coefficient


---

### 3. Probabilistically assign parental side
Here, we use the mtDNA MVS from the validation cohort to assign the parental side of the remaining relatives. To run this part, use the script `pipeline/step3_parental_side_determination/step2_mtDNA/src/mtdna_predictions.R`.

The main output file produced by this script is `data/mtdna_accuracy_derived_prediction.uniq_per_degree.txt`. It contains the following columns:
- target : ID of the focal individual
- degree : degree of the relative(s)
- maternal\_group : maternal group determination
- maternal\_probability : probability of maternal group determination



In addition, the script outputs one plot per degree of relatedness summarizing the benchmark and prediction of this approach. It also outputs one plot summarizing the mtDNA MVS predictions across all relatedness degree. You can find these in the `Plots/` folder. An example of each plot is included below.




**Plots/mtdna_accuracy_derived_prediction.degree_3rd.jpg**: validation and prediction from 3rd degree relatives

![](https://github.com/RJHFMSTR/THORIN/blob/main/pipeline/step3_parental_side_determination/step2_mtDNA/Plots/mtdna_accuracy_derived_prediction.degree_3rd.jpg?raw=true)


**Plots/KGP_mtdna_predictions.jpg**: prediction distribution across all relatedness degrees

![](https://github.com/RJHFMSTR/THORIN/blob/main/pipeline/step3_parental_side_determination/step2_mtDNA/Plots/KGP_mtdna_predictions.jpg?raw=true)




