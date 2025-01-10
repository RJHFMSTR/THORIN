---
layout: default
title: Usage
nav_order: 1
parent: Documentation
---

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---

### Usage


---

#### A simple run

To simply run the THORIN tool on your cohort, we assume you already have the following files ready:
- `data.group` : a group file listing focal individuals and reference individuals
- `data_chr20.bcf` : a genotype file containing at least your focal and reference individuals
- `data_chr20.unrelated_samples.bcf` : a genotype file containing unrelated individuals
- `chr20.b38.gmap.gz` : a genetic map file in the correct genome built (can be found [here](https://github.com/RJHFMSTR/THORIN/tree/main/maps))


<div class="code-example" markdown="1">
```bash
CHR=20
GRP=data.group
IN=data_chr20.bcf
UNR=data_chr20.unrelated_samples.bcf
MAP=chr20.b38.gmap.gz
OUT=data_chr20.thorin.prob

./thorin_v1.2 -I ${IN} -H ${UNR} -M ${MAP} -R chr${CHR} -G ${GRP} -O ${OUT}

```
</div>

---

#### The Variant Call Format output
For the output to be in `.vcf.gz` or `.bcf` format, simply add the desired extension to the output file name, such that:

<div class="code-example" markdown="1">
```bash
CHR=20
GRP=data.group
IN=data_chr20.bcf
UNR=data_chr20.unrelated_samples.bcf
MAP=chr20.b38.gmap.gz
OUT=data_chr20.thorin.prob.bcf

./thorin_v1.2 -I ${IN} -H ${UNR} -M ${MAP} -R chr${CHR} -G ${GRP} -O ${OUT}

```
</div>


---

#### IBD segments
To output IBD segments, use in addition the option `--ibd` :

<div class="code-example" markdown="1">
```bash
CHR=20
GRP=data.group
IN=data_chr20.bcf
UNR=data_chr20.unrelated_samples.bcf
MAP=chr20.b38.gmap.gz
OUT=data_chr20.thorin.prob.bcf
OUT_SEG=data_chr20.thorin_segments.prob

./thorin_v1.2 -I ${IN} -H ${UNR} -M ${MAP} -R chr${CHR} -G ${GRP} -O ${OUT} --ibd ${OUT_SEG}

```
</div>

---

#### IBD scaffold
To output the scaffold based on IBD probabilities, use the option `--phasing`:

<div class="code-example" markdown="1">
```bash
CHR=20
GRP=data.group
IN=data_chr20.bcf
UNR=data_chr20.unrelated_samples.bcf
MAP=chr20.b38.gmap.gz
OUT=data_chr20.thorin.prob.bcf
OUT_SCAF=data_chr20.thorin_scaffold.vcf.gz

./thorin_v1.2 -I ${IN} -H ${UNR} -M ${MAP} -R chr${CHR} -G ${GRP} -O ${OUT} --phasing ${OUT_SCAF}

```
</div>














