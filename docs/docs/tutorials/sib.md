---
layout: default
title: 3. Cross-over inference
nav_order: 3
parent: 4. Parental side inference
has_children: true

---
# Parental side inference from cross-overs inference in siblings

{: .no_toc }

## Table of contents
{: .no_toc .text-delta }

1. TOC
{:toc}

---



## Rationale

In this step, we aim to assign the parent-of-origin of individual using an IBD-based crossover inference approach. The underlying idea is to overlap inferred crossover event with sex-specific genetic maps to derive the likelihood of a crossover occuring in a female and in a male individuals, pointing toward maternal or maternal inheritance of the crossover, respectively.

For this, we used only siblings, which are the only relationship that can confidently be used to ensure that the crossover event identified occured in the parent, and not in more distant relatives.

In this approach, we favor inter-chromosomally phased data, since it allow identifying the entire set of crossover inherited from the same parent across the 22 autosomes. Conversly, when using intra-chromosomally phased data only, we can obtain the set of crossovers inherited from the same parent within the same chromosome only, which drastically reduces the accuracy of this approach (see our manuscript for the benchmark of both approaches).

For this step, we will first use THORIN to map IBD between each sibling pair of our data. Then, we will use custom R scripts to derive a metric that we termed "sibling score", which directly gives us a probabilistic estimation of the parent-of-origin of haplotypes.

All scripts for this step are in folder `pipeline/step3_parental_side_determination/step1_sibling_score/`.

**Note** : When running this approach on the 1,000 Genome Project data set (KGP), we have only a few individuals with both a sibling and available parental genomes. This validation cohort, that we use to derive the probabilistic assignment, is not large enough. For the example on the KGP data, we increase this set by adding random individuals, just for the code to run smoothly. This is indicated in the R script.

---

## Pipeline

### IBD mapping

We use THORIN to map IBD between sibling pairs that were identified in a previous step (see section cluster relatives). In this step, we will first reformat sibling pairs to fit THORIN input files, and then run THORIN to map IBD between siblings.



<div class="code-example" markdown="1">
```bash

mkdir -p data/

# format sibling groups
Rscript src/format_sib_group.R


# Run THORIN per chromosome
BIN=THORIN/bin/thorin_v1.2_static # adapt path
GRP=data/sibs_with_phasing.group
ODIR=data/THORIN
mkdir -p ${ODIR}
threads=1

for CHR in {1..22}; do # you can also use xargs to parallelize, or sbatch depending on your cluster.

	IN=../../step2_interchromosomal_phasing/data/SHAPEIT5/KGP.chr${CHR}.thorin_shapeit5.bcf # here, we use the inter-chromosomally phased data
	UNR=../../step1_surrogate_parents/data/relatedness/unrelated/KGP.chr${CHR}.gsa.unrelated.bcf
	MAP=/home/rhofmeis/Dropbox/Ressources/Git_repository/THORIN/maps/chr${CHR}.b38.gmap.gz

	OUT=${ODIR}/KGP.chr${CHR}.thorin.prob

	if [ ! -f "${OUT}.ibd" ]; then
		${BIN} -I ${IN} -H ${UNR} -M ${MAP} -R chr${CHR} -G ${GRP} -O ${OUT} --ibd ${OUT}.ibd -T ${threads}
	fi
done

```
</div>




---

### IBD-based crossover sibling score

In this step, we will use the IBD segment produced in the previous step to infer crossover event and overlap them with sex-specific genetic maps to compute the sibling score. Make sure to first download sex-specific genetic map from the [official github page](https://github.com/cbherer/Bherer_etal_SexualDimorphismRecombination) from Bherer et al.

The main scripts for this step are `pipeline/step3_parental_side_determination/step1_sibling_score/src/sib_score.R` and `pipeline/step3_parental_side_determination/step1_sibling_score/src/utils.R`.

To run this step, simply navigate in folder `pipeline/step3_parental_side_determination/step1_sibling_score/` and run:


<div class="code-example" markdown="1">
```bash
Rscript src/sib_score.R
```
</div>


As notified above, when using the KGP data, we dont have enough individuals in the validation cohort. In the R script `sib_score.R`, we use in addition random individual from our data, just to have enough data for the code to run to make prediction. These prediction will obviously be biologically incorrect.



---

When running the above code, it will also produce the two following figures in folder `pipeline/step3_parental_side_determination/step1_sibling_score/Plots/`, which are the equivalent of the figure from our manuscript.


**The benchmark figure**: sib\_score.benchmark.jpg

![](https://github.com/RJHFMSTR/THORIN/blob/main/pipeline/step3_parental_side_determination/step1_sibling_score/Plots/sib_score.benchmark.jpg?raw=true)

**The prediction figure**: sib\_score.call.jpg

![](https://github.com/RJHFMSTR/THORIN/blob/main/pipeline/step3_parental_side_determination/step1_sibling_score/Plots/sib_score.call.jpg?raw=true)






















