
# get sib pairs
sibs<-read.table('../../step1_surrogate_parents/data/Sibs.list', hea=F)
sibs <- as.data.frame(t(apply(sibs, 1, sort)))
sibs<-unique(sibs)
sibs<-data.frame(V1=c(sibs$V1, sibs$V2), V2=c(sibs$V2, sibs$V1))

# get inter-chromosomally scaffolded individuals
d<-data.frame()
for (CHR in 1:22){
	tmp<-as.data.frame(data.table::fread(paste0('../../step2_interchromosomal_phasing/data/THORIN/KGP.chr',CHR,'.thorin.prob.ibd.scaffolded_samples.txt'), hea=F))
	d<-rbind(d, tmp)
}
targets<-unique(d$V1)
sibs1<-sibs[!(sibs$V1 %in% targets),]
sibs2<-sibs[sibs$V1 %in% targets,]

# format for THORIN
#sibs1$V2<-paste0(sibs1$V2,'=', sibs1$V2)
sibs1$V2<-paste0('G1=', sibs1$V2)
write.table(sibs1,paste0('data/sibs_wo_phasing.group'), quote=F, col.names=F, row.names=F)
#sibs2$V2<-paste0(sibs2$V2,'=', sibs2$V2)
sibs2$V2<-paste0('G1=', sibs2$V2)
write.table(sibs2,paste0('data/sibs_with_phasing.group'), quote=F, col.names=F, row.names=F)

