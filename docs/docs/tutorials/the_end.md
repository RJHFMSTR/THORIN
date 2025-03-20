---
layout: default
title: 8. The End
nav_order: 8
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

That's it, you are at the end of our parent-of-origin inference tutorial.

The previous step allowed you to produce three files:
- `${output}.paternal\_haplotype.vcf.gz`
- `${output}.maternal\_haplotype.vcf.gz`
- `${output}.differential\_haplotype.vcf.gz`


If you wish to run association testing, using for example [REGENIE](https://rgcgithub.github.io/regenie/), you can simply convert these file to `.bgen` or `.pgen` format using plink2 for example.


In our parent-of-origin effects analysis, we typically used the "differential" haplotype to obain what we refer to as a POE p-value and identify parent-of-origin effects. In addition, we also typically run a GWAS for the paternal haplotype and for the maternal haplotype, which allows us to classify the POEs into maternal effect, paternal effect, bi-polar effect, or a more complex POE pattern.

---


For any question on this tutorial, you can contact Robin Hofmeister at robin.j.hofmeister@gmail.com or robin.hofmeister@unil.ch

---















