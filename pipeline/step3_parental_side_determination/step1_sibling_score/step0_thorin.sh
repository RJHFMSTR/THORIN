#!bin/bash

mkdir -p data/

# format sibling groups
Rscript src/format_sib_group.R

# run thorin
threads=2
printf "%s\n" {1..22} > list_chrs.txt
cat list_chrs.txt | xargs -P ${threads} -n 1 bash src/thorin_ibd.sh
rm list_chrs.txt


