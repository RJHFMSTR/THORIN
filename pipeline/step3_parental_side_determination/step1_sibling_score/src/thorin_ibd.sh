#!bin/bash

CHR=$1
threads=1

IN=../../step2_interchromosomal_phasing/data/SHAPEIT5/KGP.chr${CHR}.thorin_shapeit5.bcf # here, we use the inter-chromosomally phased data
UNR=../../step1_surrogate_parents/data/relatedness/unrelated/KGP.chr${CHR}.gsa.unrelated.bcf
MAP=/home/rhofmeis/Dropbox/Ressources/Git_repository/THORIN/maps/chr${CHR}.b38.gmap.gz
GRP=data/sibs_with_phasing.group
BIN=/home/rhofmeis/Dropbox/Ressources/Git_repository/THORIN/bin/thorin_v1.2_static

ODIR=data/THORIN
mkdir -p ${ODIR}
OUT=${ODIR}/KGP.chr${CHR}.thorin.prob

if [ ! -f "${OUT}.ibd" ]; then
	${BIN} -I ${IN} -H ${UNR} -M ${MAP} -R chr${CHR} -G ${GRP} -O ${OUT} --ibd ${OUT}.ibd -T ${threads}
fi


