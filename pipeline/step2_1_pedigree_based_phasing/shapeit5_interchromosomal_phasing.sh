#!bin/bash

threads=2


PED=data/pedigree.ped
SAMP=data/pedigree.kids.samples
rm -f ${PED}
rm -f ${SAMP}
cat ../step1_surrogate_parents/data/Trios.ped >> ${PED}
cat ../step1_surrogate_parents/data/Duos.ped >> ${PED}
cat ../step1_surrogate_parents/data/Trios.ped | cut -f 1 >> ${SAMP}
cat ../step1_surrogate_parents/data/Duos.ped | cut -f 1 >> ${SAMP}




BIN=/home/rhofmeis/Dropbox/Ressources/Git_repository/shapeit5/phase_common/bin/phase_common_static
ODIR=data/SHAPEIT5
mkdir -p ${ODIR}

for CHR in 20; do

	# You may need to re-compute AC and AN in you data, as required by shapeit5.
	IN=../step0_download_genotype/data/genotype/vcf/KGP.chr${CHR}.gsa.vcf.gz
	OUT=../step0_download_genotype/data/genotype/vcf/KGP.chr${CHR}.gsa.tags.bcf
	#bcftools +fill-tags ${IN} -Ob -o ${OUT} -- -t AC,AN
	#bcftools index ${OUT}

	IN=../step0_download_genotype/data/genotype/vcf/KGP.chr${CHR}.gsa.tags.bcf
	MAP=/home/rhofmeis/Dropbox/Ressources/Git_repository/THORIN/maps/chr${CHR}.b38.gmap.gz
	OUT=${ODIR}/KGP.chr${CHR}.pedigree_shapeit5.bcf

	if [ ! -f "${OUT}.csi" ]; then
		${BIN} -I ${IN} -M ${MAP} -O ${OUT} -R chr${CHR} -T ${threads} --pedigree ${PED}
		bcftools index ${OUT} --threads ${threads}
	
		bcftools view -S ${SAMP} -Oz -o ${ODIR}/KGP.chr${CHR}.pedigree_shapeit5.kids.vcf.gz ${ODIR}/KGP.chr${CHR}.pedigree_shapeit5.bcf
		bcftools index ${ODIR}/KGP.chr${CHR}.pedigree_shapeit5.kids.vcf.gz
	fi
done

