#!bin/bash


IN=related.chr20.vcf.gz
UNR=unrelated.chr20.vcf.gz
MAP=chr20.b37.gmap.gz
GRP=Trios.txt
OUT=test.prob

BIN=../bin/thorin_v1.2_static # change with your own version

./${BIN} -I ${IN} -H ${UNR} -M ${MAP} -R 20 -G ${GRP} -O ${OUT} --ibd ${OUT}.ibd --scaffold ${OUT}.ibd.scaffold.bcf

