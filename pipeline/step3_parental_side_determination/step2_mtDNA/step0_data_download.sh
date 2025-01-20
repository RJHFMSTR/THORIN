#!bin/bash


## Here for the example on the KGP data, we will directly start with the VCF files.
## However, in the our original study on the UK Biobank whole-genome sequencing data, we re-called mtDNA variant using the software MitoHPC.


# Download VCF files
mkdir -p data/
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz
mv ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz data/



