
REL="../../step1_surrogate_parents/data/Relatives.group"

group<-as.data.frame(data.table::fread(REL, hea=F))

# format relative pairs.
res<-data.frame()
for (i in 1:dim(group)[1]){
  tmp<-group[i,]
  target=tmp$V1
  g1<-unlist(lapply(gsub('G1=','',tmp$V2), FUN=function(x){unlist(strsplit(x,';'))}))
  for (g in g1){res<-rbind(res, data.frame(target=target, relative=g, group='G1'))}
  if (tmp$V3!=''){
    g2<-unlist(lapply(gsub('G2=','',tmp$V3), FUN=function(x){unlist(strsplit(x,';'))}))
    for (g in g2){res<-rbind(res, data.frame(target=target, relative=g, group='G2'))}
    
  }
}
group<-res



# For the KGP example, filter target and relatives to keep only those with available data.
t_vcf<-as.data.frame(data.table::fread('data/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz', hea=T))
group<-group[group$target %in% colnames(t_vcf) & group$relative %in% colnames(t_vcf),]


# function to get genotype from the vcf.
get_geno<-function(IID){
  #t_vcf<-as.data.frame(data.table::fread('data/ALL.chrMT.phase1_samtools_si.20101123.snps.low_coverage.genotypes.vcf.gz', hea=T))
  t_vcf<-as.data.frame(data.table::fread('data/ALL.chrMT.phase3_callmom-v0_4.20130502.genotypes.vcf.gz', hea=T))
  t_vcf<-cbind(t_vcf[1:5], t_vcf[colnames(t_vcf)==IID])
  colnames(t_vcf)<-c('CHR','POS','ID','REF','ALT','IID')
  t_vcf$ID<-paste0(t_vcf$CHR,'_', t_vcf$POS,'_', t_vcf$REF,'_', t_vcf$ALT)
  
  t_vcf$IID<-unlist(lapply(as.character(t_vcf$IID), FUN=function(x){unlist(strsplit(x,':'))[1]}))
  t_vcf$IID1<-unlist(lapply(as.character(t_vcf$IID), FUN=function(x){as.numeric(unlist(strsplit(x,'/'))[1])}))
  t_vcf$IID2<-unlist(lapply(as.character(t_vcf$IID), FUN=function(x){as.numeric(unlist(strsplit(x,'/'))[2])}))
  t_vcf$GENO<-unlist(apply(t_vcf, 1, FUN=function(x){sum(c(as.numeric(x[7]),as.numeric(x[8])),na.rm=T)}))
  return(t_vcf)
}




# process in turn each target-relative pair.
res<-data.frame()
for (i in 1:dim(group)[1]){
  print(i)
  tmp<-group[i,]
  target=tmp$target; relative=tmp$relative

  target_vcf<-get_geno(target)
  relative_vcf<-get_geno(relative)
  
  comp<-data.frame(ID=unique(c(target_vcf$ID, relative_vcf$ID)))
  comp$target<-target_vcf$GENO[match(comp$ID, target_vcf$ID)]
  comp$relative<-relative_vcf$GENO[match(comp$ID, relative_vcf$ID)]
  comp[is.na(comp)]<-0
  comp<-comp[comp$target!=0 | comp$relative!=0,]
  
  if (dim(comp)[1]>2){	
    c<-cor.test(comp$target, comp$relative)
    tmp$R=c$estimate; tmp$P=c$p.value; tmp$N=dim(comp)[1]; tmp$N_diff=dim(comp[comp$target!=comp$relative,])[1]
  }
  if (dim(comp)[1]<=2){
    tmp$R='U'; tmp$P='U'; tmp$N=dim(comp)[1]; tmp$N_diff=dim(comp[comp$target!=comp$relative,])[1]
  }
  
  res<-rbind(res, tmp)
}


# Add relatedness degree
kin<-as.data.frame(data.table::fread('../../step1_surrogate_parents/data/relatedness/KGP.king_relatedness.kin0', hea=T))
kin$pid1<-paste0(kin$ID1,kin$ID2)
kin$pid2<-paste0(kin$ID2,kin$ID1)

res$id<-paste0(res$target,res$relative)
res1<-res[res$id %in% kin$pid1,]; res1$degree<-kin$InfType[match(res1$id, kin$pid1)]
res2<-res[res$id %in% kin$pid2,]; res2$degree<-kin$InfType[match(res2$id, kin$pid2)]
res<-rbind(res1, res2)

# add parental side for the validation cohort subset
REL="../../step1_surrogate_parents/data/benchmark/Relatives.benchmark.side"
rel<-read.table(REL, hea=T)
res$side<-rel$side[match(paste0(res$target,res$relative), paste0(rel$target, rel$relative))]
res$mtdna<-1-(res$N_diff/res$N)
res<-res[,c(1:3,9:11)]
write.table(res, 'data/mtDNA_MVS.txt', quote=F, col.names=T, row.names=F, sep='\t')



