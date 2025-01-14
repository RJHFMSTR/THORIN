---
layout: default
title: Toy example
nav_order: 2
parent: Tutorials
---
# Toy example
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---



## Rationale
We provided test files so that you understand the structure of command line options and verify that your code is running, especially if you compiled THORIN yourself on a different platform.


## Toy example
Simply navigate in the folder `test/` and run the following code (also provided in the bash script `test/run.sh`).

<div class="code-example" markdown="1">
```bash
IN=related.chr20.vcf.gz
UNR=unrelated.chr20.vcf.gz
MAP=chr20.b37.gmap.gz
GRP=Trios.txt
OUT=test.prob

BIN=../bin/thorin_v1.2_static # change with your own version

./${BIN} -I ${IN} -H ${UNR} -M ${MAP} -R 20 -G ${GRP} -O ${OUT} --ibd ${OUT}.ibd --scaffold ${OUT}.ibd.scaffold.bcf
```
</div>


This code gives you 3 outputs:
- test.prob : the HMM IBD probability per variant site
- test.prob.ibd : the IBD segments and corresponding class (see documentation).
- test.prob.ibd.scaffold.bcf : the IBD-based scaffold in VCF/BCF format that you can use as input in the phasing software [SHAPEIT5](https://odelaneau.github.io/shapeit5/) to perform inter-chromosomal phasing.


---

If you want the probability per variant site in VCF/BCF format, you can simply add the suffix to the output file name, such as:

<div class="code-example" markdown="1">
```bash
IN=related.chr20.vcf.gz
UNR=unrelated.chr20.vcf.gz
MAP=chr20.b37.gmap.gz
GRP=Trios.txt
OUT=test.prob.vcf.gz

BIN=../bin/thorin_v1.2_static # change with your own version

./${BIN} -I ${IN} -H ${UNR} -M ${MAP} -R 20 -G ${GRP} -O ${OUT} --ibd ${OUT}.ibd --scaffold ${OUT}.ibd.scaffold.bcf
```
</div>

---




