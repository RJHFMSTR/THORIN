#!bin/bash

CHR=$1

threads=2
DIR="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"
PRE="CCDG_14151_B01_GRM_WGS_2020-08-05"
SUF="filtered.shapeit2-duohmm-phased.vcf.gz"
ODIR=data/genotype
mkdir -p ${ODIR}/vcf



if [ ${CHR} == "X" ]; then
        SUF="filtered.eagle2-phased.v2.vcf.gz"
fi



if [ ! -f "${ODIR}/vcf/${PRE}_chr${CHR}.gsa.${SUF}.csi" ]; then
	
	# download chromosome-wide vcf.gz file
        wget ${DIR}/${PRE}_chr${CHR}.${SUF} -O ${ODIR}/vcf/${PRE}_chr${CHR}.${SUF}
        wget ${DIR}/${PRE}_chr${CHR}.${SUF}.tbi -O ${ODIR}/vcf/${PRE}_chr${CHR}.${SUF}.tbi

	# subset to GSA sites
	bcftools view -T data/GSA.variant_positions.txt ${ODIR}/vcf/${PRE}_chr${CHR}.${SUF} -Oz -o ${ODIR}/vcf/${PRE}_chr${CHR}.gsa.${SUF} --threads ${threads}
	bcftools index -f ${ODIR}/vcf/${PRE}_chr${CHR}.gsa.${SUF} --threads ${threads}
	rm ${ODIR}/vcf/${PRE}_chr${CHR}.${SUF}*

fi



if [ ! -f "${ODIR}/plink/KGP.chr${CHR}.gsa.bed" ]; then

	# convert to plink file
	mkdir -p ${ODIR}/plink
	plink2 --vcf ${ODIR}/vcf/${PRE}_chr${CHR}.gsa.${SUF} --make-bed --out ${ODIR}/plink/KGP.chr${CHR}.gsa

	mv ${ODIR}/vcf/${PRE}_chr${CHR}.gsa.${SUF} ${ODIR}/vcf/KGP.chr${CHR}.gsa.vcf.gz
        mv ${ODIR}/vcf/${PRE}_chr${CHR}.gsa.${SUF}.csi ${ODIR}/vcf/KGP.chr${CHR}.gsa.vcf.gz.csi
fi


