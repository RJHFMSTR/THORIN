#!bin/bash

CHR=$1
threads=2



# You may need to re-compute AC and AN in you data, as required by shapeit5.
IN=../step0_download_genotype/data/genotype/vcf/KGP.chr${CHR}.gsa.vcf.gz
OUT=../step0_download_genotype/data/genotype/vcf/KGP.chr${CHR}.gsa.tags.bcf
#bcftools +fill-tags ${IN} -Ob -o ${OUT} -- -t AC,AN
#bcftools index ${OUT}

IN=data/THORIN/KGP.chr${CHR}.thorin.prob.ibd.scaffold.bcf
OUT=data/THORIN/KGP.chr${CHR}.thorin.prob.ibd.scaffold.tags.bcf
#bcftools +fill-tags ${IN} -Ob -o ${OUT} -- -t AC,AN
#bcftools index ${OUT}
#



IN=../step0_download_genotype/data/genotype/vcf/KGP.chr${CHR}.gsa.tags.bcf
MAP=/home/rhofmeis/Dropbox/Ressources/Git_repository/THORIN/maps/chr${CHR}.b38.gmap.gz
SCAF=data/THORIN/KGP.chr${CHR}.thorin.prob.ibd.scaffold.tags.bcf

BIN=/home/rhofmeis/Dropbox/Ressources/Git_repository/shapeit5/phase_common/bin/phase_common_static

ODIR=data/SHAPEIT5
mkdir -p ${ODIR}
OUT=${ODIR}/KGP.chr${CHR}.thorin_shapeit5.bcf

if [ ! -f "${OUT}.csi" ]; then
	${BIN} -I ${IN} -M ${MAP} -O ${OUT} -R chr${CHR} -T ${threads} -S ${SCAF}
	bcftools index ${OUT} --threads ${threads}
fi


