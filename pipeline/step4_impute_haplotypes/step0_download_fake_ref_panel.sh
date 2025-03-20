#!bin/bash



# For the example, we started with the KGP WGS data, that we subsetted to GSA sites, simulating a SNP-array data.
# Here, we will simulate imputation.
# Since we need publicly available data for this example, we will use the KGP WGS data as reference panel.

# Note 1: we will do this only on chromosome 20.

# Note 2: we do this only for the example, but this is completely wrong to do on your real data. Please use a proper reference panel for imputation, depending on your input data. You can obviously not use the same data as input and reference panel to perform imputation. The goal here is just to provide example script that are running.


CHR=20

threads=2
DIR="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"
PRE="CCDG_14151_B01_GRM_WGS_2020-08-05"
SUF="filtered.shapeit2-duohmm-phased.vcf.gz"
ODIR=data/reference_panel
mkdir -p ${ODIR}




if [ ! -f "${ODIR}/${PRE}_chr${CHR}.${SUF}.csi" ]; then
	
	# download chromosome-wide vcf.gz file
        wget ${DIR}/${PRE}_chr${CHR}.${SUF} -O ${ODIR}/${PRE}_chr${CHR}.${SUF}
        wget ${DIR}/${PRE}_chr${CHR}.${SUF}.tbi -O ${ODIR}/${PRE}_chr${CHR}.${SUF}.tbi

fi



