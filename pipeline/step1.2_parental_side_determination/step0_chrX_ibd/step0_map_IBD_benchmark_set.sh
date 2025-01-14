#!bin/bash

CHR=X

IN=../step0_download_genotype/data/genotype/vcf/KGP.chr${CHR}.gsa.diploidized.bcf
UNR=../step1_surrogate_parents/data/relatedness/unrelated/KGP.chr${CHR}.gsa.unrelated.bcf
MAP=../../maps/chrX.b38.gmap.gz
GRP=../step1_surrogate_parents/data/benchmark/Relatives.benchmark.group

BIN=../../bin/thorin_v1.2_static


ODIR=data/THORIN/benchmark
mkdir -p ${ODIR}
OUT=${ODIR}/KGP.chr${CHR}.benchmark.thorin.prob


./${BIN} -I ${IN} -H ${UNR} -M ${MAP} -R chr${CHR} -G ${GRP} -O ${OUT} --ibd ${OUT}.ibd









