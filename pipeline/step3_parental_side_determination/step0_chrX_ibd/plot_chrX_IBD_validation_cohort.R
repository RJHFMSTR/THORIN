
# Read the file with relatives and assigned parental side
rel<-read.table('../../step1_surrogate_parents/data/benchmark/Relatives.benchmark.side', hea=T)
rel<-rel[rel$target_sex=='Male',]
rel<-rel[,c(1,3,4)]; rel<-rel[!duplicated(rel),]
rel2<-data.frame()
for (t in unique(rel$target)){
  print(t)
  tmp<-rel[rel$target==t,]
  if (dim(tmp)[1]==1){
    
    if (tmp$group=='G1' & tmp$side=='maternal'){
      rel2<-rbind(rel2, data.frame(target=tmp$target, group='G2', side='paternal'))
    }
    if (tmp$group=='G1' & tmp$side=='paternal'){
      rel2<-rbind(rel2, data.frame(target=tmp$target, group='G2', side='maternal'))
    }
    if (tmp$group=='G2' & tmp$side=='maternal'){
      rel2<-rbind(rel2, data.frame(target=tmp$target, group='G1', side='paternal'))
    }
    if (tmp$group=='G2' & tmp$side=='paternal'){
      rel2<-rbind(rel2, data.frame(target=tmp$target, group='G1', side='maternal'))
    }
  }
}
rel<-rbind(rel, rel2)


# Assign parental side to IBD tracks
dt<-read.table('data/THORIN/benchmark/KGP.chrX.benchmark.thorin.prob.ibd', hea=T)
dt<-dt[dt$target %in% rel$target & dt$Prob!='C',]
dt$copy_group[dt$Prob=='A']<-'G1'; dt$copy_group[dt$Prob=='B']<-'G2'
dt$side<-rel$side[match(paste0(dt$target,':',dt$copy_group),paste0(rel$target,':',rel$group))]
dt<-dt[!is.na(dt$side),]


# Plot IBD sharing per relative parental side
png('chrX_IBD_validation_cohort.png')
par(mfrow=c(1,2))
boxplot(dt$length_CM~dt$side, xlab='relative parental side', ylab='IBD segment length')
boxplot(dt$length_CM~dt$side, xlab='relative parental side', ylab='IBD segment length', ylim=c(0,10))
dev.off()
