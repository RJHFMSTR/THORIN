#!bin/bash


CHR=20

mkdir -p data/Haplotypes/

VCF=../step4_impute_haplotypes/data/IMPUTE5/test.KGP.chr${CHR}.impute5.info08.maf01.vcf.gz # imputed, scaffold-phased
OUT=data/Haplotypes/KGP.chr${CHR}.impute5.info08.maf01 # output file prefix
PROB=../step3_parental_side_determination/step3_combined_predictors/PofO_probability.txt # PofO probability assignment

python3 src/encode.py -i ${VCF} -o ${OUT} -p ${PROB}

bcftools index -f ${OUT}.paternal_haplotype.vcf.gz
bcftools index -f ${OUT}.maternal_haplotype.vcf.gz
bcftools index -f ${OUT}.differential_haplotype.vcf.gz

