---
layout: default
title: A. Using mate-pairs
nav_order: 1
parent: Genetic Assortative Mating pipeline tutorial
---
# Genetic Assortative Mating pipeline tutorial
{: .no_toc }

## Table of contents
{: .no_toc .text-delta }


---



## Rationale

When computing genetic assortative mating (GAM) using our haplotype-based method, we likely also want to compare those estimates to a set of established mate-pair.

In this section, we will see how to infer mate-pairs in the UK Biobank cohort, and how to use these pairs to estimate GAM.

For the Estonian Biobank, we simply use the large (>10k) number of mate-pairs that can directly be inferred from genetically determined parent-offspring trios (see [here](https://rjhfmstr.github.io/THORIN/docs/tutorials/cluster_relatives.html) for the tutorial).

---


## Mate pair inference

Following [Yengo et al., 2018](https://www.nature.com/articles/s41562-018-0476-3), we infer mate-pairs using a combination of phenotypes from the UK Biobank:

- home location east 1km : code 22702
- townsend_deprivation index : code 22189
- home location north 1km : code 22704
- average total household income before tax : code 738
- smokers_in_household : code 1259
- number_in_household : code 709
- lengt of time at current address : code 699
- how are people in household related to participant : code 6141


Assuming you already have :
- a phenotype matrix containing these phenotypes
- a file listing relate dindividuals (typically the output of the KING software, provided by the UKBB; or you can also compute it following [this tutorial](https://rjhfmstr.github.io/THORIN/docs/tutorials/cluster_relatives.html#infer-relatedness))
- a file with the participants' age (extracted from year-of-birth is more accurate than "age at recruitment")
- a file with the participants' genetically inferred sex


You can infer mate-pairs using the script [provided on github](https://github.com/RJHFMSTR/GAM/blob/main/src/cluster_ukbb_mate_pairs.R) and included below: 


<div class="code-example" markdown="1">
```bash
# Create a data with all phenotypes used to determine mate pairs. Use the following fields
#Phenotype_Name	Phenotype_Code
#home_location_east_1km 22702
#townsend_deprivation_index 22189
#home_location_north_1km 22704
#Average_total_household_income_before_tax 738
#smokers_in_household 1259
#number_in_household 709
#lengt_of_time_at_current_address 699
#How_are_people_in_household_related_to_participant 6141

phenotype_file="Out/Mate_phenotypes.txt"
relatedness_file="/data/FAC/FBM/DBC/zkutalik/default_sensitive/uk_biobank/useful_data/ukb_relatedness.up_to_degree_4.txt.gz" # the KING file provided by UKBB
age_file="/data/FAC/FBM/DBC/zkutalik/default_sensitive/uk_biobank/useful_data/ukb_age.txt.gz" # col1=IID; col2=age
wb_file="/data/FAC/FBM/DBC/zkutalik/default_sensitive/uk_biobank/useful_data/ukb_caucasians_for_regenie.txt"

d<-as.data.frame(data.table::fread(phenotype_file, hea=T))

# select only indiv reporting living with mate.
d<-d[d$living_with_husband_or_wife_or_partner==1,]

# remove indiv with missing data.
d<-d[!is.na(d$Average_total_household_income_before_tax),]
d<-d[!is.na(d$home_location_east_1km),]
d<-d[!is.na(d$townsend_deprivation_index),]
d<-d[!is.na(d$home_location_north_1km),]
d<-d[!is.na(d$Average_total_household_income_before_tax),]
d<-d[!is.na(d$smokers_in_household),]
d<-d[!is.na(d$number_in_household),]
d<-d[!is.na(d$lengt_of_time_at_current_address),]
d<-d[!is.na(d$living_with_husband_or_wife_or_partner),]
d<-d[!is.na(d$genetic_sex),]

# remove non-respondant
d<-d[d$Average_total_household_income_before_tax!=1 & d$Average_total_household_income_before_tax!=-3,]
d<-d[d$smokers_in_household!=-3,]
d<-d[d$number_in_household>0,]
d<-d[d$lengt_of_time_at_current_address!=-1 & d$lengt_of_time_at_current_address!=-3,]

# create unique identifier for all phenotypes
d$pheno_comb<-paste0(d$home_location_east_1km,'_',d$home_location_north_1km,'_',d$townsend_deprivation_index,'_',d$smokers_in_household,'_',d$number_in_household,'_',d$lengt_of_time_at_current_address,'_',d$Average_total_household_income_before_tax)

# get individual with same identifier
cluster_pairs<-function(pheno) {
  individuals <- d$IID[d$pheno_comb == pheno]
  if (length(individuals) > 1) {
    pairs <- t(combn(individuals, 2))
    return(data.frame(ID1 = pairs[,1], ID2 = pairs[,2], pheno_comb = pheno, stringsAsFactors = FALSE))
  }
}
library(parallel); library(dplyr)
x<-mclapply(unique(d$pheno_comb), cluster_pairs, mc.cores=36)
result<-rbind(bind_rows(x))

# add genetic sex
result$sex_id1<-d$genetic_sex[match(result$ID1, d$IID)]
result$sex_id2<-d$genetic_sex[match(result$ID2, d$IID)]

# remove non-opposite mate pairs
res<-result[result$sex_id1 != result$sex_id2,]

# add genetic relatedness
rel<-as.data.frame(data.table::fread(relatedness_file, hea=T))
rel_pairs<-unique(c(paste0(rel$ID1,'_',rel$ID2), paste0(rel$ID2,'_',rel$ID1)))

# remove related pairs
res$pair_id<-paste0(res$ID1, '_', res$ID2)
res<-res[!(res$pair_id %in% rel_pairs),]

# add age
age<-as.data.frame(data.table::fread(age_file, hea=T))
res$age_id1<-age$age[match(res$ID1, age$IID)]
res$age_id2<-age$age[match(res$ID2, age$IID)]
res<-res[complete.cases(res),]

# remove pairs with >=10 years diff
res$diff_age<-abs(res$age_id1-res$age_id2)
res<-res[res$diff_age<=10,]

# remove non white British (field 22006)
cau<-as.data.frame(data.table::fread(wb_file, hea=F))
res<-res[res$ID1 %in% cau$V1,]
res<-res[res$ID2 %in% cau$V1,]

# remove pairs for which one of the mate appear in >1 pair
t<-table(c(res$ID1, res$ID2)); n<-names(t)[unname(t)>1]
res<-res[!(res$ID1 %in% n) & !(res$ID1 %in% n),]

write.table(res, 'Mate_pairs.txt', quote=F, col.names=T, row.names=F, sep='\t')

```
</div>






## Mate-pair polygenic score (PGS) calculation

Once you have identify the mate-pairs, you have two alternative to estimate assortative mating:
- use the true phenotypic value. For this, you'll need to have access to the biobank phenotypes of interest.
- use polygenic scores. This is what we will do in this tutorial

As part of our [preprint](https://www.biorxiv.org/content/10.1101/2025.09.24.678243v1), we provide the PGS file that were derived from a subset of ~ 150,000 UK Biobank individuals, not overlapping mate-pair individuals. These score files can be accessed on our [zenodo repository](https://doi.org/10.5281/zenodo.17600467).


To compute PGS for the mate-pairs individuals, you'll need:
- the score files (argument ${SCORE}
- the genotype file for the mate-pairs (argument ${VCF})
- the [pgs-calc software](https://pgsc-calc.readthedocs.io/en/latest/)


### PGS computation:

<div class="code-example" markdown="1">
```bash
PHECODE=21001
SCORE=${PHECODE}.betas.tsv.gz
OUT=${PHECODE}.mate_pairs.txt
pgs-calc apply --ref ${SCORE} ${VCF} --out ${OUT} --threads 1

```
</div>




## Genetic assortative mating estimates

The last step to estimate GAM from mate-pair data consists in computing the correlation between mate-pairs PGS.
For this, you will need:
- the mate-pair inferred in previous step
- the per individual PGS, one score per individual, genome-wide. If you computed PGS for each chromosome separately, you can sum up the 22 scores to get a final whole-genome score per individual.
- a file with Principal Components (PCs)



Using these files, you can estimate GAM from mate-pairs using the script [provided on github](https://github.com/RJHFMSTR/GAM/blob/main/src/compute_gam_mate_pairs.R) and included below:

<div class="code-example" markdown="1">
```bash
library(data.table)
library(parallel)

N_CORE <- 32

## Input files
p_file <- "Phenotype_manifest.txt"    # should contain at least columns: PGS, PHENO (or similar)
m_file <- "Mate_pairs.txt"            # two columns: IDs of partners
c_file <- "Covariates.txt.gz"         # columns: ID, PC1–PC10, etc.

## Read phenotype manifest
pheno <- fread(p_file)
## adapt these if your column names are different
## e.g. pheno$PGS and pheno$PHENO
pgs_list <- pheno$PGS

## Read mate pairs
couples <- fread(m_file)
setnames(couples, 1:2, c("ID1", "ID2"))

## Read covariates
cov <- fread(c_file)  # assumes a column "ID" and PCs PC1–PC10

pc_vars <- paste0("PC", 1:10)

run <- function(pgs) {

  ## Read mate-pair PGS file for this trait
  infile <- paste0(pgs, ".mate_pairs.txt") # one score per individual for the 22 autosomes.
  scores <- fread(infile)
  setnames(scores, 1:2, c("ID", "SCORE"))

  ## Merge PGS for both partners
  d <- merge(couples, scores, by.x = "ID1", by.y = "ID", all.x = TRUE)
  setnames(d, "SCORE", "score_id1")
  d <- merge(d, scores, by.x = "ID2", by.y = "ID", all.x = TRUE)
  setnames(d, "SCORE", "score_id2")

  ## Add PCs for partner 1
  d <- merge(d, cov[, c("ID", pc_vars), with = FALSE],
             by.x = "ID1", by.y = "ID", all.x = TRUE)

  ## Keep complete cases for scores + PCs
  d_cc <- d[complete.cases(score_id1, score_id2, d[, ..pc_vars])]

  N_pairs <- nrow(d_cc)

  ## Raw correlation between partners' PGS
  cor_res <- cor.test(d_cc$score_id1, d_cc$score_id2)

  ## PC-adjusted model: score_id2 ~ score_id1 + PC1–PC10
  formula_pc <- as.formula(
    paste("score_id2 ~ score_id1 +", paste(pc_vars, collapse = " + "))
  )
  fit_pc <- lm(formula_pc, data = d_cc)
  cf_pc  <- coef(summary(fit_pc))
  ss_pc  <- cf_pc["score_id1", ]

  out <- data.frame(
    PGS   = pgs,
    PHENO = pheno$PHENO[match(pgs, pheno$PGS)],   # adapt PHENO column name if needed
    N     = N_pairs,
    R     = unname(cor_res$estimate),
    RP    = cor_res$p.value,
    beta_pc = unname(ss_pc["Estimate"]),
    se_pc   = unname(ss_pc["Std. Error"]),
    p_pc    = unname(ss_pc["Pr(>|t|)"])
  )

  return(out)
}

## Parallel run over all PGS
results_list <- mclapply(pgs_list, run, mc.cores = N_CORE)
results <- rbindlist(results_list)
fwrite(results, "Mate_pairs_GAM.txt", sep = "\t")
```
</div>

































