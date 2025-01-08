#!bin/bash


threads=16
ODIR=data
mkdir -p ${ODIR}


###
##### 1. download and pre-filter genotype data
###



##
## GSA variant, use to subset the WGS 30X files to SNP-array sites to simulate large genotyped biobanks (e.g. UK Biobank, Estonian Biobank, Finngen)
##

# for hg38
wget https://emea.support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/global-screening-array-24/infinium-global-screening-array-24-v1-0-c2-annotated.zip -O ${ODIR}/gsa.zip
unzip data/gsa.zip -d ${ODIR}/ && rm ${ODIR}/gsa.zip
awk 'NR > 1 {print "chr"$2"\t"$3}' ${ODIR}/GSA-24v1-0_C2.hg38.annotated.txt > ${ODIR}/GSA.variant_positions.txt
rm ${ODIR}/GSA-24v1-0_C2.hg38.annotated.txt


##
## Download and filter genotype files
##
mkdir -p ${ODIR}/genotype
printf "%s\n" {1..22} X > list_chrs.txt
cat list_chrs.txt | xargs -P ${threads} -n 1 bash src/download_chr.sh
rm list_chrs.txt


##
## merge all chromomsome into a single file for KING
##

printf "%s\n" ${ODIR}/genotype/plink/KGP.chr{1..22}.gsa > merge_list.txt
plink1.9 --merge-list merge_list.txt --make-bed --out data/genotype/plink/KGP.merged_chromosomes 
rm ${ODIR}/genotype/plink/KGP.chr*.gsa.*
rm merge_list.txt

##
## modify fam file so that they don't have all the same family ID == 0.
##
awk '{print $2" "$2" "$3" "$4" "$5" "$6}' data/genotype/plink/KGP.merged_chromosomes.fam > data/genotype/plink/KGP.merged_chromosomes.v2.fam
mv data/genotype/plink/KGP.merged_chromosomes.v2.fam data/genotype/plink/KGP.merged_chromosomes.fam


##
## get age and sex (here, age is not available).
##
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/1kGP.3202_samples.pedigree_info.txt -O data/1kGP.3202_samples.pedigree_info.txt

