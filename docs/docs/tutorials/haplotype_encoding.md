---
layout: default
title: 7. Parental haplotype encoding
nav_order: 7
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
In this section, we will start with our haploid-based imputed data, that we will recode into paternal haplotype, maternal haplotype and differential haplotype using of PofO predictions.

---


## Pipeline
This step rely on a single script that you can found in folder `pipeline/step5_encode_parental_haplotypes/src/encode.py`.

It takes as arguments: 

| Option name 	       | Argument| Default  | Description |
|:---------------------|:--------|:---------|:-------------------------------------|
| \-i \[\-\-input\_vcf \]  | STRING  | NA       | Haploid-based imputed data in vcf.gz format |
| \-p \-\-prob\_file     | STRING  | NA       | Parent-of-origin assignment probability file (see step 5 of the tutorial |
| \-o \-\-output  | STRING  | NA       | Output file prefix |


This script uses the PofO assignment probability produced in step 5 of the tutorial (specifically, it uses the file `PofO_probability.txt` as input) to assign allele to maternal and paternal haplotypes.
Additionally, it encode the differential haplotype. For this, it uses only heterozygous sites, and encode the data as 1 if the allele in paternally inherited, and 0 if maternally inherited. All homozygous are set as missing.



To do this, you can simply run the following command and adapt the input files to your data:

<div class="code-example" markdown="1">
```bash 
bash step0_encode_haplotypes.sh
```
</div>



The above command will produce three output files in `.vcf.gz``format:
- `${output}.paternal\_haplotype.vcf.gz`
- `${output}.maternal\_haplotype.vcf.gz`
- `${output}.differential\_haplotype.vcf.gz`



---















