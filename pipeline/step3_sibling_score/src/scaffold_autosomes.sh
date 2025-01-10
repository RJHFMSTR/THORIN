#!bin/bash

CHR=$1

IN=../step0_download_genotype/data/genotype/vcf/KGP.chr${CHR}.gsa.vcf.gz
UNR=../step1_surrogate_parents/data/relatedness/unrelated/KGP.chr${CHR}.gsa.unrelated.bcf
MAP=../../maps/chr${CHR}.b38.gmap.gz
GRP=data/Relatives.group
BIN=../../bin/thorin_v1.2_static


ODIR=data/THORIN/autosomes/from_relatives
mkdir -p ${ODIR}
OUT=${ODIR}/KGP.chr${CHR}.thorin.prob


./${BIN} -I ${IN} -H ${UNR} -M ${MAP} -R chr${CHR} -G ${GRP} -O ${OUT} --ibd ${OUT}.ibd --phasing ${OUT}.scaffold.bcf

bcftools index -f ${OUT}.scaffold.bcf







