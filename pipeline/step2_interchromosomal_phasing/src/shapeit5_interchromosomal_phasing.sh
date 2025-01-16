#!bin/bash

CHR=$1
threads=1

IN=../step0_download_genotype/data/genotype/vcf/KGP.chr${CHR}.gsa.vcf.gz
MAP=/home/rhofmeis/Dropbox/Ressources/Git_repository/THORIN/maps/chr${CHR}.b38.gmap.gz
SCAF=data/THORIN/KGP.chr${CHR}.thorin.prob.ibd.scaffold.bcf

BIN=/home/rhofmeis/Dropbox/Ressources/Git_repository/shapeit5/phase_common/bin/phase_common_static

ODIR=data/SHAPEIT5
mkdir -p ${ODIR}
OUT=${ODIR}/KGP.chr${CHR}.thorin_shapeit5.vcf.gz

if [ ! -f "${OUT}.csi" ]; then
	${BIN} -I ${IN} -M ${MAP} -O ${OUT} -R ${CHR} -T ${threads} -S ${SCAF}
	bcftools index ${OUT} --threads ${threads}
fi


