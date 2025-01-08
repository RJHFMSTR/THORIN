
# I added the --ouput-bcf option
# This file verify that the --output-bcf option matches
# the ouptut of the --output

OUT=output_test.txt
OUT_BCF=output_test.bcf
TRIO=../test/Trios.txt

../bin/thorin_v1.2 -I ../test/related.chr20.vcf.gz -H ../test/unrelated.chr20.vcf.gz -G $TRIO -M ../maps/chr20.b37.gmap.gz -R 20 -O $OUT 
../bin/thorin_v1.2 -I ../test/related.chr20.vcf.gz -H ../test/unrelated.chr20.vcf.gz -G $TRIO -M ../maps/chr20.b37.gmap.gz -R 20 -O $OUT_BCF

# check number of records
echo "Check number of records (numbers should be equal)"
grep -v "#" $OUT | wc -l
bcftools view $OUT_BCF -H | wc -l

# check CHROM POS
echo "Check CHROM and POS per row (should return nothing)"
awk 'NR>1{print $1, $2}' $OUT > tmp.txt
bcftools view $OUT_BCF -H | awk -F '\t' '{print $1, $2}' > tmp.bcf
diff tmp.txt tmp.bcf
rm tmp.txt tmp.bcf

# check values for each sample
echo "Check probabilities per row (should return nothing)"
awk -F '\t' 'NR>1{for(i=5;i<=NF;i++) printf "%s%s", $i*100/100, (i==NF ? "\n" : " ")}' $OUT  > tmp.txt # I divide by 100 and multiply by 100 so that I get "1" instead of "1.00"
bcftools view $OUT_BCF -H | awk -F '\t' '{for(i=10;i<=NF;i++) printf "%s%s", $i, (i==NF ? "\n" : " ")}' > tmp.bcf
diff tmp.txt tmp.bcf
rm tmp.txt tmp.bcf

# check samples names
echo "Check samples name (should return nothing)"
head -n 1 $OUT | awk -F '\t' '{for(i=5;i<=NF;i++) printf "%s%s", $i, (i==NF ? "\n" : " ")}'  > tmp.txt 
bcftools view $OUT_BCF -h | tail -1  | awk -F '\t' '{for(i=10;i<=NF;i++) printf "%s%s", $i, (i==NF ? "\n" : " ")}' > tmp.bcf
diff tmp.txt tmp.bcf
rm tmp.txt tmp.bcf
 
# rm files
rm $OUT $OUT_BCF
