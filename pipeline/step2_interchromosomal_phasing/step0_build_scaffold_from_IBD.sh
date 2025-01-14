#!bin/bash

CHR=20

IN=KGP.chr${CHR}.gsa.vcf.gz
UNR=KGP.chr${CHR}.gsa.unrelated.bcf
MAP=chr${CHR}.b38.gmap.gz
GRP=Relatives.benchmark.group


BIN=./../../../../bin/thorin_v1.2_static


OUT=test.thorin.prob


./${BIN} -I ${IN} -H ${UNR} -M ${MAP} -R chr${CHR} -G ${GRP} -O ${OUT} --ibd ${OUT}.ibd --scaffold ${OUT}.ibd.scaffold.bcf




