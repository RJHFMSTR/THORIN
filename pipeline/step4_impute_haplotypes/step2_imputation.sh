#!bin/bash



###
# IMPUTE5 documentation: https://jmarchini.org/software/#impute-5
# IMPUTE 5 is freely available for academic use only
###


# Note: for the example on the KDP project, we will run imputation on the entire chromosome. When dealing with larger dataset, you may want tto run the imputation per chunks speed up the process by parallelizing. When using chunks, specify the region you want to impute (here, IRG), and also specify a buffer region (here, BRG). The buffer region is typically a 250kb window aroung the specified input region.


CHR=20
threads=4


ODIR=data/IMPUTE5
mkdir -p ${ODIR}
IN=${ODIR}/KGP.chr${CHR}.targets.scaffold_phased_shapeit5.pofo_merged.bcf
REF=data/reference_panel/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz
MAP=/home/rhofmeis/Dropbox/Ressources/Git_repository/THORIN/maps/chr${CHR}.b38.gmap.gz
OUT=${ODIR}/KGP.chr${CHR}.impute5.bcf
OUT_F=${ODIR}/KGP.chr${CHR}.impute5.info08.maf01.vcf.gz


if [ ! -e "${OUT_F}.csi" ]; then
	STA=60137 # start of the input region to impute
	END=242195 # end of the input region to impute
	BSTA=50137 # start of the buffer region
	BEND=342195 # end of the buffer region
	IRG=chr${CHR}:${STA}-${END}	
	BRG=chr${CHR}:${BSTA}-${BEND}
	./impute5_static_v1.2.1 --g ${IN} --h ${REF} --m ${MAP} --r ${IRG} --buffer-region ${BRG} --o ${OUT} --l ${OUT}.log --threads ${threads} --out-ap-field
	bcftools filter -e 'INFO/INFO<0.8' --threads ${threads} ${OUT} -Ou | bcftools view --threads ${threads} -q 0.01:minor -Oz -o ${OUT_F}
	bcftools index ${OUT_F} --threads ${threads}
fi




