---
layout: default
title: 5. Parent-or-origin determination genome-wide
nav_order: 5
parent: Parent-of-origin inference pipeline tutorial
has_children: true

---
# Parent-or-origin determination genome-wide using inter-chromosomally phased data 

{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---



## Rationale
In this section, the aim is to assign the relatives of a focal individual to a parental side (maternal or paternal) using all the pre-computed PofO predictors (i.e chrX, sib-score and mtDNA).

---


## Pipeline

This step rely on a single script that you can found in folder `pipeline/step3_parental_side_determination/step3_combined_predictors/`.
As explain in our manuscript, when a focal individual as several available predictors, we keep the parental side indicated by the predictor having the highest accuracy. If two predictors yielded the same accuracy, we verify that they indicate the same parental side.

To perform this step, simply run the following command:


<div class="code-example" markdown="1">
```bash
Rscript combined_probabilities.R

```
</div>



This script will give you the main output `PofO_probability.txt`, where:
- col1 = focal individual ID
- col2 = maternal group (G1 or G2)
- col3 = probability that the group indicated in col2 is maternal (between 0.5 and 1)
- col4 = predictor(s) that were used to assign the parental side
- col5 = probability that the group G1 was maternal


In addition, the script also outputs a two panels probability plot in folder `Plots/`. Panel a show the number of individuals (y-axis) per probability bin. Panel b shows the number of individuals (y-axis) per PofO predictor (x-axis) for individuals within the probability bin 0.99-1.

---



