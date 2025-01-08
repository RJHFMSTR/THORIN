#!bin/bash


BIN=diploidize_static # available from https://github.com/odelaneau/otools
SEX=data/1kGP.3202_samples.pedigree_info.txt
VCF=data/genotype/vcf/KGP.chrX.gsa.vcf.gz
threads=8


##
## Get male individuals
##
awk '$4 == 2 {print $1}' data/1kGP.3202_samples.pedigree_info.txt > data/KGP.males.txt



##
## QC male individuals --> transform mixed diploid and haploid call into fully diploid calls.
##

IN=data/genotype/vcf/KGP.chrX.gsa.vcf.gz

OUT_M=data/genotype/vcf/KGP.chrX.gsa.males.bcf
bcftools view -S data/KGP.males.txt -Ob -o ${OUT_M} ${IN} --threads ${threads}
bcftools index ${OUT_M} --threads ${threads}


OUT_F=data/genotype/vcf/KGP.chrX.gsa.females.bcf
bcftools view -S ^data/KGP.males.txt -Ob -o ${OUT_F} ${IN} --threads ${threads}
bcftools index ${OUT_F} --threads ${threads}

OUT_D=data/genotype/vcf/KGP.chrX.gsa.males.diploidized.bcf
./${BIN} --input ${OUT_M} --output ${OUT_D}
bcftools index ${OUT_D} --threads ${threads}

OUT=data/genotype/vcf/KGP.chrX.gsa.diploidized.bcf
bcftools merge ${OUT_D} ${OUT_F} -Ob -o ${OUT}
bcftools index ${OUT} --threads ${threads}

rm ${OUT_D}* && rm ${OUT_M}* && rm ${OUT_F}*



