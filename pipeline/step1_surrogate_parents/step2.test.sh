#!bin/bash

threads=8

# get 100 random unrelated
UNR=data/relatedness/unrelated/random_unrelated_samples.txt



# subset each genotype file to keep only unrelated samples
for CHR in {1..22}; do
	IN=../step0_download_genotype/data/genotype/vcf/KGP.chr${CHR}.gsa.vcf.gz
	OUT=data/relatedness/unrelated/KGP.chr${CHR}.gsa.unrelated.bcf
	bcftools view -S ${UNR} -Ob -o ${OUT} ${IN} --threads ${threads}
	bcftools index ${OUT} --threads ${threads}
done

for CHR in X; do
	IN=../step0_download_genotype/data/genotype/vcf/KGP.chr${CHR}.gsa.diploidize.bcf
	OUT=data/relatedness/unrelated/KGP.chr${CHR}.gsa.unrelated.bcf
	bcftools view -S ${UNR} -Ob -o ${OUT} ${IN} --threads ${threads}
	bcftools index ${OUT} --threads ${threads}
done

