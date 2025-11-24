---
layout: default
title: B. Using inter-chromosomally phased data
nav_order: 2
parent: Genetic Assortative Mating pipeline tutorial
---
# Genetic Assortative Mating pipeline tutorial
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }


---



## Rationale

Estimating GAM from inter-chromosomally phased data is conceptually very similar to using mate-pairs data; the only difference in the input data.

Instead of using two individual forming a mate-pair, we use the two haplotype of the same individual as proxy for paternal and maternal genomes. The mate-pair is then formed by the haplotype 1 and the haplotype 2 of the same individual.
Practically, it means that, for each individual, we will use haplotype 1 and haplotype 2 as if these were two different individuals:
- we compute a PGS for each haplotype
- we correlate the haplotype 1's PGS and haplotype 2's PGS (instead of mate 1 and mate 2) to estimate GAM.

---


## Inter-chromosomal phasing and parental haplotype segregation

The inter-chromosomal phasing pipeline is available as part of our [parent-of-origin inference tutorial](https://rjhfmstr.github.io/THORIN/docs/tutorials/interchrs.html). Follow this tutorial up to [step 6; parental haplotypes imputation](https://rjhfmstr.github.io/THORIN/docs/tutorials/genotype_imputation.html).

Once you have imputed you parental haplotypes, we simply need to split the genotype file in .vcf.gz format into two separate files, one for each parental haplotypes.
This can be done using the python script [provided on github](https://github.com/RJHFMSTR/THORIN/blob/main/alternative_projects/GAM/split_parental_haplotypes.py) and pasted below:




<div class="code-example" markdown="1">
```bash
import argparse
import gzip
import os
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input_vcf')
parser.add_argument('-h1', '--output_haplotype1')
parser.add_argument('-h2', '--output_haplotype2')
args = parser.parse_args()

vcf_in=args.input_bcf
out_h1=args.output_haplotype1
out_h2=args.output_haplotype2

outfile1=open(out_h1,'w')
outfile2=open(out_h2,'w')
with gzip.open(vcf_in,'rb') as f:
    for line1 in f:
        line=line1.decode()
        tmp=line.split()
        if line.find('#')==0:
            outfile1.write(line)
            outfile2.write(line)

        else:
            w1=tmp[:8]; w1.append('GT:DS:GP')
            w2=tmp[:8]; w2.append('GT:DS:GP')
            for fields in tmp[9:]:

                GT=fields.split(':')[0]
                GT1=GT.split('|')[0]
                GT2=GT.split('|')[1]

                AP=fields.split(':')[2]
                AP1=AP.split(',')[0]
                AP2=AP.split(',')[1]

                GP1=str((1-float(AP1)))+','+AP1+',0'
                GP2=str((1-float(AP2)))+','+AP2+',0'

                field_h1=GT1+'|0:'+AP1+':'+GP1
                field_h2=GT2+'|0:'+AP2+':'+GP2

                w1.append(field_h1)
                w2.append(field_h2)

            outfile1.write('\t'.join(w1)+'\n')
            outfile2.write('\t'.join(w2)+'\n')

outfile1.close()
outfile2.close()

os.system('bgzip -f '+out_h1)
os.system('bgzip -f '+out_h2)
os.system('bcftools index -f '+out_h1+'.gz')
os.system('bcftools index -f '+out_h2+'.gz')

```
</div>





## Haploid-based PGS
The previous step produces two vcf.gz file per chromosome, one for each parental haplotype.
Using the score files available on our [zenodo repository](https://doi.org/10.5281/zenodo.17600467), you can compute the PGS for each parental haplotype of each individual:

<div class="code-example" markdown="1">
```bash

PHECODE=21001
SCORE=${PHECODE}.betas.tsv.gz

# Compute PGS for each haplotype
OUT=${PHECODE}_haploid_based_PGS.txt
for HAP in 1 2; do
	VCF=haplotype${HAP}.vcf.gz
	pgs-calc apply --ref ${SCORE} ${VCF} --out PGS_haplotype${HAP}.txt --threads 1
done

# Combine output into a single file with col1=ID, col2=PGS_h1, col3=PGS_h2
awk -F',' 'NR==FNR {if(FNR>1){gsub(/"/,""); hap1[$1]=$2} next}
           FNR==1 {next}
           {gsub(/"/,"");
            if(!p){print "sample,hap1_score,hap2_score"; p=1}
            print $1","hap1[$1]","$2}' \
  <(cat PGS_haplotype1.txt) \
  <(cat PGS_haplotype2.txt) \
> ${OUT}

```
</div>






# Haploid-based GAM estimates
Now, we can simply use the PGS derived from each haplotype to estimate GAM, as implemented in the R script provided [on github](https://github.com/RJHFMSTR/THORIN/blob/main/alternative_projects/GAM/compute_haploid_gam.R) and pasted below:

<div class="code-example" markdown="1">
```bash
library(data.table)
library(parallel)

N_CORE <- 32

## Input files
p_file <- "Phenotype_manifest.txt"    # should contain at least columns: PGS, PHENO (or similar)
c_file <- "Covariates.txt.gz"         # columns: ID, PC1–PC10, etc.

## Read phenotype manifest
pheno <- fread(p_file)
## adapt these if your column names are different
pgs_list <- pheno$PGS

## Read covariates
cov <- fread(c_file)  # assumes a column "ID" and PCs PC1–PC10
pc_vars <- paste0("PC", 1:10)

run <- function(pgs) {

  ## Read haplotype-based PGS file for this trait
  ## Adapt the path/name if needed (e.g. include chr or other components)
  infile <- paste0( pgs, "_haploid_based_PGS.txt")
  scores <- fread(infile)
  setnames(scores, 1:3, c("ID", "PGS_h1", "PGS_h2"))  # sample, hap1_score, hap2_score

  ## hap1_score and hap2_score must be genome-wide scores. If you run the PGS per chromosome, make sure to sum up the scores for the 22 autosomes into a single genome-wide score.
  ## if you have missing chromosome for some individuals, replace the score for the missing chromosome using the average score for that chromosome computed from the reste of the data.

  ## Merge with covariates
  d <- merge(scores, cov[, c("ID", pc_vars), with = FALSE],
             by = "ID", all.x = TRUE)

  ## Keep complete cases for hap1, hap2, and PCs
  d_cc <- d[complete.cases(PGS_h1, PGS_h2, d[, ..pc_vars])]

  N_ind <- nrow(d_cc)

  ## Raw correlation between haplotypes
  cor_res <- cor.test(d_cc$PGS_h1, d_cc$PGS_h2)

  ## PC-adjusted model: PGS_h2 ~ PGS_h1 + PC1–PC10
  formula_pc <- as.formula(
    paste("PGS_h2 ~ PGS_h1 +", paste(pc_vars, collapse = " + "))
  )
  fit_pc <- lm(formula_pc, data = d_cc)
  cf_pc  <- coef(summary(fit_pc))
  ss_pc  <- cf_pc["PGS_h1", ]

  out <- data.frame(
    PGS    = pgs,
    PHENO  = pheno$PHENO[match(pgs, pheno$PGS)],   # adapt PHENO column name if needed
    N      = N_ind,
    R      = unname(cor_res$estimate),
    RP     = cor_res$p.value,
    beta_pc = unname(ss_pc["Estimate"]),
    se_pc   = unname(ss_pc["Std. Error"]),
    p_pc    = unname(ss_pc["Pr(>|t|)"])
  )

  return(out)
}

## Parallel run over all PGS
results_list <- mclapply(pgs_list, run, mc.cores = N_CORE)
results <- rbindlist(results_list)
fwrite(results, "Haplotype_GAM.txt", sep = "\t")

```
</div>
















