#!bin/bash



PFX="KGP"
REL="data/relatedness/KGP.king_relatedness.kin0"
SEX="../step0_download_genotype/data/1kGP.3202_samples.pedigree_info.txt"
AGE="../step0_download_genotype/data/1kGP.3202_samples.pedigree_info.txt" # for the KGP example, we don't have the age of individual. Here, we provided the same file as for the sex, only for the example.

Rscript src/grouping.R ${PFX} ${REL} ${SEX} ${AGE}

