
REL_VCF=../test/related.chr20.vcf.gz
UNREL_VCF=../test/unrelated.chr20.vcf.gz
PHASE_VCF=phased_vcf.bcf
IBD=ibd.txt
TRIO=../test/Trios.txt

#------------------------------------------------------------------------------#
# RUN thorin
#------------------------------------------------------------------------------#
../bin/thorin_v1.2 -I $REL_VCF -H $UNREL_VCF -G $TRIO -M ../maps/chr20.b37.gmap.gz -R 20 -O output_test.bcf --ibd $IBD --scaffold $PHASE_VCF --scaffold-cM 0

#------------------------------------------------------------------------------#
# VALIDATION
#------------------------------------------------------------------------------#
sample=NA20900

# Filter non biallelic sites
filtered_rel_vcf=filtered.related.chr20.bcf
bcftools view -m2 -M2  $REL_VCF -Ob -o $filtered_rel_vcf && bcftools index -f $filtered_rel_vcf

# Check that original and new have the same alleles
echo "------------------------------------------------------------------------------
Check if the genotypes are the same (should not print anything)
------------------------------------------------------------------------------"
bcftools view $filtered_rel_vcf -s $sample | bcftools query -f "[%GT]\n" | sed "s/|/ /g" | awk '{print $1+$2}'> ori_genotypes.txt
bcftools view $PHASE_VCF -s $sample | bcftools query -f "[%GT]\n" | sed "s/|/ /g" | sed "s/\// /g" | awk '{print $1+$2}'> new_genotypes.txt
diff ori_genotypes.txt new_genotypes.txt
rm ori_genotypes.txt new_genotypes.txt

# Check that positions with "B" prob have indeed a changed phase
echo "------------------------------------------------------------------------------
Check where we had wrong phasing i.e. where whe had Prob=B we should have
inversed genotypes (should not print anything)
------------------------------------------------------------------------------"
bcftools view $filtered_rel_vcf -s $sample | bcftools query -f "%POS [%GT]\n" | sed "s/|/ /g" | awk '{print $1,$2,$3}' > ori_phase.txt
bcftools view $PHASE_VCF -s $sample | bcftools query -f "%POS [%GT]\n" | sed "s/|/ /g" | sed "s/\// /g" | awk '{print $1,$2,$3}' > new_phase.txt
paste ori_phase.txt new_phase.txt | awk '{ if ($2!=$5) print $1 }' > changed_phase.txt
python Scripts/check_phasing.py
rm ori_phase.txt new_phase.txt changed_phase.txt

# Check phasing regions in the new file
echo "------------------------------------------------------------------------------
Check that Prob=C and Prob=D positions are not phased '/' + Check that Prob=A 
and Prob=B positions are phased '|'
------------------------------------------------------------------------------"
bcftools view $PHASE_VCF -s $sample | bcftools query -f "%POS [%GT]\n" | awk '{if($2 ~ /\//) {print $1,"/"} else {print $1,"|"}}' > unphased_vs_phased.txt
python Scripts/check_unphased_vs_phased.py
rm unphased_vs_phased.txt

#------------------------------------------------------------------------------#
# RUN thorin but with scaffold != 0
#------------------------------------------------------------------------------#
../bin/thorin_v1.2 -I $REL_VCF -H $UNREL_VCF -G $TRIO -M ../maps/chr20.b37.gmap.gz -R 20 -O output_test.bcf --ibd $IBD --scaffold $PHASE_VCF --scaffold-cM 3

#------------------------------------------------------------------------------#
# Validation
#------------------------------------------------------------------------------#

# Check that positions with "B" prob > 3cM have indeed a changed phase
echo "------------------------------------------------------------------------------
Check where we had wrong phasing i.e. where whe had Prob=B we should have
inversed genotypes (should not print anything)
------------------------------------------------------------------------------"
bcftools view $filtered_rel_vcf -s $sample | bcftools query -f "%POS [%GT]\n" | sed "s/|/ /g" | awk '{print $1,$2,$3}' > ori_phase.txt
bcftools view $PHASE_VCF -s $sample | bcftools query -f "%POS [%GT]\n" | sed "s/|/ /g" | sed "s/\// /g" | awk '{print $1,$2,$3}' > new_phase.txt
paste ori_phase.txt new_phase.txt | awk '{ if ($2!=$5) print $1 }' > changed_phase.3cM.txt
python Scripts/check_phasing.3cM.py
rm ori_phase.txt new_phase.txt changed_phase.3cM.txt

# Check phasing regions in the new file
echo "------------------------------------------------------------------------------
Check that Prob=C and Prob=D positions are not phased '/' + Check that Prob=A 
and Prob=B positions are phased '|'
------------------------------------------------------------------------------"
bcftools view $PHASE_VCF -s $sample | bcftools query -f "%POS [%GT]\n" | awk '{if($2 ~ /\//) {print $1,"/"} else {print $1,"|"}}' > unphased_vs_phased.3cM.txt
python Scripts/check_unphased_vs_phased.3cM.py
rm unphased_vs_phased.3cM.txt


