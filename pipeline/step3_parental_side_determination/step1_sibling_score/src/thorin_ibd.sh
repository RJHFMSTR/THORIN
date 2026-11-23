#!/bin/bash

CHR=$1
threads=1

IN=../../step2_interchromosomal_phasing/data/SHAPEIT5/KGP.chr${CHR}.thorin_shapeit5.bcf # here, we use the inter-chromosomally phased data
UNR=../../step1_surrogate_parents/data/relatedness/unrelated/KGP.chr${CHR}.gsa.unrelated.bcf
MAP=../../../maps/chr${CHR}.b38.gmap.gz # adapt path
GRP=data/sibs_with_phasing.group
BIN=../../../bin/thorin_v1.3_static # adapt path. Sibling labels (A/B/C/D) must be read with the current ibd2co (src/utils.R)

ODIR=data/THORIN
mkdir -p ${ODIR}
OUT=${ODIR}/KGP.chr${CHR}.thorin.prob

if [ ! -f "${OUT}.ibd" ]; then
	${BIN} -I ${IN} -H ${UNR} -M ${MAP} -R chr${CHR} -G ${GRP} -O ${OUT} --ibd ${OUT}.ibd -T ${threads}
fi


