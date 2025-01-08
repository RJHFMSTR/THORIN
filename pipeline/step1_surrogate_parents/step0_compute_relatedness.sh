#!bin/bash

wget https://www.kingrelatedness.com/Linux-king.tar.gz
tar -xvzf Linux-king.tar.gz
chmod +x king

BED=../step0_download_genotype/data/genotype/plink/KGP.merged_chromosomes.bed
PFX=data/relatedness/KGP.king_relatedness
mkdir -p data/relatedness


./king -b ${BED} --related --degree 5 --cpus 4 --prefix ${PFX}



