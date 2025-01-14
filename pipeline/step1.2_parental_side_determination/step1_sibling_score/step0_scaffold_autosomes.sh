#!bin/bash


mkdir -p data/


# Keep only the group file for individual having a sibling
awk 'NR==FNR {sibs[$1]; next} $1 in sibs' ../step1_surrogate_parents/data/Sibs.list ../step1_surrogate_parents/data/Relatives.group > data/Relatives.group

printf "%s\n" {1..22} X > list_chrs.txt
cat list_chrs.txt | xargs -P ${threads} -n 1 bash src/scaffold_autosomes.sh
rm list_chrs.txt














