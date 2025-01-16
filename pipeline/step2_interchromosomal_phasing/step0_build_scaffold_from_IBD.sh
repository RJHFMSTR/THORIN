#!bin/bash


threads=2

printf "%s\n" {1..22} > list_chrs.txt
cat list_chrs.txt | xargs -P ${threads} -n 1 bash src/thorin_scaffold.sh
rm list_chrs.txt

