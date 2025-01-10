---
layout: default
title: Basic usage
nav_order: 1
parent: Documentation
---

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

### Description
We explain here the main utility and usage of the THORIN tool to map IBD between a focal individual and reference individual(s).

---

### Reference individual(s)
The THORIN tool maps IBD between a focal individual and reference individual(s). The reference individual(s) for each focal individual are listed in a separate `.group` file. This file is a tab-separated file, each line corresponding to a different focal individual and it's reference individuals. The first column lists the focal individual ID. Additional columns list the reference individual IDs. The format of the additional columns is such that the first element corresponds to the group ID `GID`, and the following elements correspond to the reference individual(s) IDs. `GID` is the element that is reported in the output file. If we want to tracks the IDs of reference individuals, we can specify `GID` as being the reference individual ID (particularly useful when using single reference individual, as seen in line 1 of the example file below). When grouping reference individuals, we can specify `GID` as being the reference group 1, reference group 2, etc.


An example of `.group` file is presented below (in practice, this file takes no header line) :


| Focal individual ID   | Reference 1                | Reference 2                  |
|:----------------------|:---------------------------|:-----------------------------|
| sample\_1             | sample\_2=sample\_2        | sample\_3=sample\_3          |
| sample\_4             | GID\_1=sample\_5;sample\_6 | GID\_2=sample\_7             |
| sample\_8             | GID\_1=sample\_9;sample\_10| GID\_2=sample\_11;sample\_12 |



The file lists here three focal individuals (`sample\_1`, `sample\_4` and `sample\_8`). Each of these focal individuals have a different group structure. For `sample\_1`, we will map IBD between single individuals and keep track of the reference individual IDs in the output file by specifying the reference individuals IDs in the `GID` fields. For `sample\_4` and `sample\_8`, we clustered some of the reference individuals into two groups and we specify the ID of each group in the `GID` field. 


---

### Unrelated individuals
For the model to work correctly, we also need to provide a set of individuals unrelated to the focal individual, so that our model will compare the probability of sharing IBD with the reference individuals to the probability of sharing IBD with unrelated individuals. The set of unrelated individuals must be unrelated to all the focal individuals listed in the group file.
Assuming that we have a file listing our set of unrelated individual, with one ID per line, we can build this file using:


<div class="code-example" markdown="1">
```bash
CHR=20
IN=data_chr20.bcf
OUT=data_chr20.unrelated_samples.bcf
SAMP=unrelated_samples.txt

bcftools view -S ${SAMP} -Ob -o ${OUT} ${IN} && bcftools index ${OUT} 
```
</div>

---

### Input data
The input data must be indexed `.vcf.gz` or `.bcf` format containing at least all the focal individuals and reference individuals listed in your `.group` file. To speed-up computation of large datasets, you can subset your cohort file to include only those individuals using:

<div class="code-example" markdown="1">
```bash
CHR=20
IN=data_chr20.bcf
OUT=data_chr20.focal_and_reference_samples.bcf
SAMP=focal_and_reference_samples.txt

bcftools view -S ${SAMP} -Ob -o ${OUT} ${IN} && bcftools index ${OUT}         
```
</div>


---
### Outputs
The THORIN tool allow for different types of outputs.

#### IBD per variant site
This is the basic output of the THORIN tool, also present in v1.0.0. It reports the probability of sharing IBD per variant site on each of the focal individual haplotypes with (i) each of the reference individuals or group listed in the `.group` file, and (ii) the set of unrelated individuals. 

An example of output file is presented below:


| #CHROM | POS | IDX | CM  | sample\_1_sample\_2_0 | sample\_1_sample\_3_0 | sample\_1_HOLE_0 | sample\_1_sample\_2_1 | sample\_1_sample\_3_1 | sample\_1_HOLE_1 | sample\_4_GID\_1_0 | sample\_4_GID\_2_0 | sample\_4_HOLE_0 | sample\_4_GID\_1_1 | sample\_4_GID\_2_1 | sample\_4_HOLE_1|
|:-------|:----|:----|:----|:----------------------|:----------------------|:-----------------|:----------------------|:----------------------|:-----------------|:-------------------|:-------------------|:-----------------|:-------------------|:-------------------|:----------------|
|:-------|: 1  | 0   | 0.0 | 0.0                   | 1.0                   | 0.0              | 0.0                   | 0.0                   | 1.0              | 0.98               | 0.0                | 0.02             | 0.0                | 1.0                | 0.0             |
|:-------|: 12 | 1   | 0.1 | 0.0                   | 0.99                  | 0.01             | 0.0                   | 0.0                   | 1.0              | 1.0                | 0.0                | 0.0              | 0.0                | 0.99               | 0.01            |



In the output, columns correspond to:
- CHROM : chromosome ID
- POS : genomic position, as indicated in the input `.vcf.gz` or `.bcf` file
- IDX : index of the variant, from 0 to N - 1.
- CM : centimorgan
- additional columns: IBD probability based on the `.group` file.


---

### Command line options

#### Basic options

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-\-help             | NA      | NA       | Produces help message |
| \-\-seed             | INT     | 15052011 | Seed of the random number generator  |
| \-T \[ \-\-thread \] | INT     | 1        | Number of thread used|

#### Input files

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-I \[\-\-input \]   | STRING  | NA       | Genotypes to be phased in VCF/BCF/XCF format |
| \-H \[\-\-reference \]| STRING  | NA       | Reference panel of haplotypes in VCF/BCF/XCF format  |
| \-S \[\-\-scaffold \]| STRING  | NA       | Scaffold of haplotypes in VCF/BCF/XCF format  |
| \-M \[\-\-map \]     | STRING  | NA       | Genetic map  |
| \-\-pedigree         | STRING  | NA       | Pedigree information (offspring father mother triplets) |
| \-\-haploids         | STRING  | NA       | List of samples that are haploids (e.g. males for chrX) |
| \-R \[\-\-region \]  | STRING  | NA       | Target region  |


#### Filter parameters

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-\-filter-snp       | NA      | NA       | If specified, the program only consider SNPs |
| \-\-filter-maf       | FLOAT   | 0        | \[Expert option\] Only consider variants with MAF above the specifed value. It requires AC/AN tags in VCF/BCF file. |


#### MCMC parameters [Expert]

| Option name 	      | Argument| Default              | Description |
|:--------------------|:--------|:---------------------|:-------------------------------------|
| \-\-mcmc-iterations | STRING  | 5b,1p,1b,1p,1b,1p,5m | Iteration scheme of the MCMC (burnin=b, pruning=p, main=m) |
| \-\-mcmc-prune      | FLOAT   | 0.999                | Pruning threshold for genotype graphs (internal memory structures)  |
| \-\-mcmc-noinit     | NA      | NA                   | If specified, phasing initialization by PBWT sweep is disabled |

#### PBWT parameters [Expert]

| Option name 	      | Argument|  Default  | Description |
|:--------------------|:--------|:----------|:-------------------------------------|
| \-\-pbwt-modulo     | FLOAT   | 0.1       | Storage frequency of PBWT indexes in cM |
| \-\-pbwt-depth      | INT     | 4         | Depth of PBWT indexes to condition on  |
| \-\-pbwt-mac        | INT     | 5         | Minimal Minor Allele Count at which PBWT is evaluated |
| \-\-pbwt-mdr        | FLOAT   | 0.1       | Maximal Missing Data Rate at which PBWT is evaluated |
| \-\-pbwt-window     | INT     | 4         | Run PBWT selection in windows of this size |

#### HMM parameters [Expert]

| Option name 	      | Argument|  Default  | Description |
|:--------------------|:--------|:----------|:-------------------------------------|
| \-\-hmm-window      | INT     | 4         | Minimal size of the phasing window in cM |
| \-\-hmm-ne          | INT     | 15000     | Effective size of the population |

#### Output files

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-O \[\-\-output \]  | STRING  | NA       | Phased haplotypes in VCF/BCF/XCF format |
| \-\-output-format     | STRING  | bcf       | File format for the output ([bcf] standard VCF/BCF format / [graph] graph format that intergrates phasing uncertainty / [bh] XCF binary format for fast loading in Impute5)  |
| \-\-log              | STRING  | NA       | Log file  |
