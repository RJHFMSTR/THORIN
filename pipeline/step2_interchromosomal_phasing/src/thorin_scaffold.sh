#!bin/bash

CHR=$1
threads=1

IN=../step0_download_genotype/data/genotype/vcf/KGP.chr${CHR}.gsa.vcf.gz
UNR=../step1_surrogate_parents/data/relatedness/unrelated/KGP.chr${CHR}.gsa.unrelated.bcf
MAP=/home/rhofmeis/Dropbox/Ressources/Git_repository/THORIN/maps/chr${CHR}.b38.gmap.gz
GRP=../step1_surrogate_parents/data/Relatives.group

BIN=/home/rhofmeis/Dropbox/Ressources/Git_repository/THORIN/bin/thorin_v1.2_static

ODIR=data/THORIN
mkdir -p ${ODIR}
OUT=${ODIR}/KGP.chr${CHR}.thorin.prob

if [ ! -f "${OUT}.ibd.scaffold.bcf.csi" ]; then
	${BIN} -I ${IN} -H ${UNR} -M ${MAP} -R chr${CHR} -G ${GRP} -O ${OUT} --ibd ${OUT}.ibd --scaffold ${OUT}.ibd.scaffold.bcf -T ${threads}
	bcftools index ${OUT}.ibd.scaffold.bcf --threads ${threads}
fi



# Get scaffolded samples

SAMP=${OUT}.ibd.scaffolded_samples.txt
awk -F'\t' 'NR > 1 && $5 >= 3 && ($4 == "A" || $4 == "B") { print $6 }' "${OUT}.ibd" | sort | uniq > ${SAMP}

