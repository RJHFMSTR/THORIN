#!bin/bash




CHR=20
threads=4

mkdir -p data/support/

# remove kids from surrogate parents callset. Keep only scaffolded samples
IN=../step2_interchromosomal_phasing/data/SHAPEIT5/KGP.chr${CHR}.thorin_shapeit5.bcf
SAMP=../step2_1_pedigree_based_phasing/data/pedigree.kids.samples
SCAF=../step2_interchromosomal_phasing/data/THORIN/KGP.chr${CHR}.thorin.prob.ibd.scaffolded_samples.txt
OUT=data/support/KGP.chr${CHR}.thorin_shapeit5.wo_kids.bcf
bcftools view --threads ${threads} --force-samples -S ^${SAMP} -Ou ${IN} | bcftools view --force-samples -S ${SCAF} -Ob -o ${OUT} && bcftools index ${OUT} --threads ${threads}


# add ped-phased kids.
IN1=../step2_1_pedigree_based_phasing/data/SHAPEIT5/KGP.chr${CHR}.pedigree_shapeit5.kids.vcf.gz
IN2=data/support/KGP.chr${CHR}.thorin_shapeit5.wo_kids.bcf
ODIR=data/IMPUTE5/
mkdir -p ${ODIR}
OUT=${ODIR}/KGP.chr${CHR}.targets.scaffold_phased_shapeit5.pofo_merged.bcf
bcftools merge --threads ${threads} ${IN1} ${IN2} -Ob -o ${OUT} && bcftools index ${OUT} --threads ${threads}
rm ${IN2}*



