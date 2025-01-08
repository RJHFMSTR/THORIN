
# Run THORIN
OUT=output_test.txt
OUT_IBD=output_test.ibd
TRIO=../test/Trios.txt

../bin/thorin_v1.2 -I ../test/related.chr20.vcf.gz -H ../test/unrelated.chr20.vcf.gz -G $TRIO -M ../maps/chr20.b37.gmap.gz -R 20 -O $OUT --ibd $OUT_IBD 

# Run v1 IBD
PROB_FILE="output_test.txt"
GROUP_FILE=$TRIO
PROB_OUT="v1_ibd.bed"
CHR="20"
N_CORE="1"
Rscript Scripts/compute_pofo_probabilities.R $PROB_FILE $GROUP_FILE $PROB_OUT $CHR $N_CORE

## Run comparison
python3 Scripts/compare_ibd_out.py $OUT_IBD $PROB_OUT
rm $OUT_IBD $PROB_OUT $OUT
