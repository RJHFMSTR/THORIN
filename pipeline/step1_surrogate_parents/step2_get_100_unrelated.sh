#!bin/bash

mkdir -p data/relatedness/unrelated
threads=8


# list related individuals.
REL=data/relatedness/KGP.king_relatedness.kin0
awk '$14 != "UN" {print $1 "\n" $3}' ${REL} | grep -v FID | sort | uniq > data/relatedness/related.samples



# list all samples
IN=../step0_download_genotype/data/genotype/vcf/KGP.chr22.gsa.vcf.gz
bcftools query -l ${IN} > data/all.samples



# get 100 random unrelated
UNR=data/relatedness/unrelated/random_unrelated_samples.txt
comm -23 <(sort data/all.samples) <(sort data/relatedness/related.samples) | shuf -n 100 > ${UNR}



# subset each genotype file to keep only unrelated samples
for CHR in {1..22}; do
	IN=../step0_download_genotype/data/genotype/vcf/KGP.chr${CHR}.gsa.vcf.gz
	OUT=data/relatedness/unrelated/KGP.chr${CHR}.gsa.unrelated.bcf
	bcftools view -S ${UNR} -Ob -o ${OUT} ${IN} --threads ${threads}
	bcftools index ${OUT} --threads ${threads}
done

for CHR in X; do
	IN=../step0_download_genotype/data/genotype/vcf/KGP.chr${CHR}.gsa.diploidized.bcf
	OUT=data/relatedness/unrelated/KGP.chr${CHR}.gsa.unrelated.bcf
	bcftools view -S ${UNR} -Ob -o ${OUT} ${IN} --threads ${threads}
	bcftools index ${OUT} --threads ${threads}
done
