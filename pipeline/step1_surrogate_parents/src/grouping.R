#install.packages(c('dplyr','data.table','igraph','parallel','argparse'), dependencies=T)
library(dplyr)
library(igraph)
library(parallel)
library(argparse)
`%ni%` = Negate(`%in%`)
N_CORE=8
args=commandArgs(trailingOnly=T)

#source('src/grouping_functions.R')
source('src/grouping_functions.KGP.R')




#**********************************
# THIS FIRST PART MUST BE MODIFY TO FIT YOUR DATAS.
#**********************************

prefix=args[1]

## 0.0 Input files
relatedness_file=args[2]
sex_file=args[3]
age_file=args[4]





# for the KGP example:
prefix='KGP'
relatedness_file="data/relatedness/KGP.king_relatedness.kin0"
sex_file="../step0_download_genotype/data/1kGP.3202_samples.pedigree_info.txt"
age_file="../step0_download_genotype/data/1kGP.3202_samples.pedigree_info.txt"




## 0.1 Output files
odir='data/'; system(paste0('mkdir -p ',odir))
out_ped=paste0(odir,prefix,'.ped') # .ped file can be used as input for pedigree phasing using the SHAPEIT5 software
out_trios=paste0(odir,'Trios.ped') # .ped file can be used as input for pedigree phasing using the SHAPEIT5 software
out_duos=paste0(odir,'Duos.ped') # .ped file can be used as input for pedigree phasing using the SHAPEIT5 software
out_relatives=paste0(odir,'Relatives.group') # .group file can be used as input to map IBD using the THORIN software
out_sibs=paste0(odir,'Sibs.list') # # .list file lists all relatives of the specified degree for each target individuals. col1=target; col2=relative list; col3=target sex; col4= number of relatives
out_relatives_males=paste0(odir,'Relatives.male_targets.group') # .group file can be used as input to map IBD using the THORIN software




## 1. Read file and rename columns

# Relatedness file (KING output)
d<-as.data.frame(data.table::fread(relatedness_file, hea=T))

# if needed, rename the columns
#colnames(d)[1]<-'ID1'
#colnames(d)[2]<-'ID2'
#colnames(d)[5]<-'Kinship'
#colnames(d)[4]<-'IBS0'
#colnames(d)[5]<-'InfType'





# Age
age<-as.data.frame(data.table::fread(age_file, hea=T)) # May need to change to hea=F if the header is included in the file. In any case, rename the colnames to "ID1" and "age".Note: For the example on KGP files, we don't have the age information. The file used for the example contain incorrect age information and can not be reliably used. It is provided only as example. 
colnames(age)[1]<-'ID1'
colnames(age)[2]<-'age'
age<-age[complete.cases(age),]





# Sex
sex<-as.data.frame(data.table::fread(sex_file, hea=T)) # Same here, may need to change to hea=F.
colnames(sex)[1]<-'ID1'
colnames(sex)[4]<-'sex'
sex$sex[sex$sex==1]<-'Male'
sex$sex[sex$sex==2]<-'Female' # change this value to 0 is coded according to plink1.9
sex<-sex[complete.cases(sex),]









#**********************************
# IN THEORY, NOTHING SHOULD BE MODIFIED IN THE FOLLOWING PART.
#**********************************

# 1. set threshold used by KING. If available, use the 'InfType' column in the relatedness file.
po_down=1/(2**(5/2))
po_up=1/(2**(3/2))
rel2_down=1/(2**(7/2))
rel2_up=1/(2**(5/2))
rel3_down=1/(2**(9/2))
rel3_up=1/(2**(7/2))
rel4_down=NA
rel4_up=1/(2**(9/2))

if (!('InfType' %in% colnames(d))){
  d$InfType<-'UN'
  d$InfType[d$Kinship>po_up]<-'Dup/MZ'
  d$InfType[d$Kinship>=po_down & d$Kinship<po_up & d$IBS0<0.0012]<-'PO'
  d$InfType[d$Kinship>=po_down & d$Kinship<po_up & d$IBS0>0.0012]<-'FS'
  d$InfType[d$Kinship>=rel2_down & d$Kinship<rel2_up]<-'2nd'
  d$InfType[d$Kinship>=rel3_down & d$Kinship<rel3_up]<-'3rd'
  d$InfType[d$Kinship>=rel4_down & d$Kinship<rel4_up]<-'4th'
}



## 2. Keep only samples for which we have at least one relative, and sex and age are available
d<-d[(d$ID1 %in% sex$ID1) & (d$ID1 %in% age$ID1) & (d$ID2 %in% sex$ID1) & (d$ID2 %in% age$ID1),]
d<-d[d$InfType!='UN',]





## 3. Cluster relatives.

# In this part, we cluster relative into surrogate parent groups. We consider different cases for (i) parent-offspring duos and trios, (ii) siblings, and (iii) >=2 degree relatives.












### I. TRIOS AND DUOS
mz<-unique(c(d$ID1[d$InfType=='Dup/MZ'], d$ID2[d$InfType=='Dup/MZ']))
rel<-d[d$InfType=='PO',]
samples<-unique(c(rel$ID1, rel$ID2))
print(paste0('Total number of individuals to consider for the clustering of trios-duos: ',length(samples)))

x<-mclapply(samples, grouping_ped, mz=mz, mc.cores=N_CORE)
xx<-cbind(bind_rows(x))
PED<-xx[!is.na(xx$target),]




### For the KGP example, since we don't have the age of individuals, we don't know who is the proband and who is the parent.
# We will therefore used previously reported trios
ped<-as.data.frame(data.table::fread("../step0_download_genotype/data/1kGP.3202_samples.pedigree_info.txt", hea=T)) # May need to change to hea=F if the header is included in the file. In any case, rename the colnames to "ID1" and "age".Note: For the example on KGP files, we don't have the age information. The file used for the example contain incorrect age information and can not be reliably used. It is provided only as example. 
ped<-ped[ped$fatherID!=0 | ped$motherID!=0,]
ped<-ped[,1:3]
ped[ped==0]<-NA
PED<-ped
colnames(PED)<-c('target','father','mother')




### II. SIBS
rel<-d[d$InfType=='FS',]
samples<-unique(c(rel$ID1, rel$ID2))
print(paste0('Total number of individuals to consider for the clustering of sibs: ',length(samples)))

x<-mclapply(samples, grouping_sibs, mc.cores=N_CORE)
SIBS<-cbind(bind_rows(x)); colnames(SIBS)<-c('target','sib')










### 2nd 3rd 4th degree relatives
rel<-d[d$InfType=='2nd' | d$InfType=='3rd' | d$InfType=='4th',]
tt<-c(rel$ID1, rel$ID2); t<-table(tt)
samples<-unique(names(t)[unname(t)<50]) # remove indiv with abnormal number of relatives
print(paste0('Total number of individuals to consider for the clustering of relatives: ',length(samples)))

x<-mclapply(samples, grouping_rel, mc.cores=N_CORE)
OUT<-cbind(bind_rows(x))







## Write final files
SIBS[SIBS=='NA']<-NA
OUT[OUT=='NA']<-NA
PED[PED=='NA']<-NA

TRIOS<-PED[complete.cases(PED),]
DUOS<-PED[PED$target %ni% TRIOS$target,]

TRIOS<-TRIOS[!is.na(TRIOS$target),]
DUOS<-DUOS[!is.na(DUOS$target),]
OUT<-OUT[!is.na(OUT$target),]
SIBS<-SIBS[!is.na(SIBS$target),]

write.table(OUT, out_relatives, quote=F, row.names=F, col.names=F, sep='\t') # .group file can be used as input to map IBD using the THORIN software
write.table(TRIOS,out_trios, quote=F, row.names=F, col.names=F, sep='\t') # .ped file can be used as input for pedigree phasing using the SHAPEIT5 software
write.table(DUOS, out_duos, quote=F, row.names=F, col.names=F, sep='\t')
write.table(SIBS, out_sibs, quote=F, row.names=F, col.names=F, sep='\t') # .list file lists all relatives of the specified degree for each target individuals. col1=target; col2=relative list; col3=target sex; col4= number of relatives
write.table(PED, out_ped, quote=F, row.names=F, col.names=F, sep='\t')

OUT<-OUT[OUT$target %in% sex$ID1[sex$sex=='Male'],]
write.table(OUT, out_relatives_males, quote=F, row.names=F, col.names=F, sep='\t') # .group file can be used as input to map IBD using the THORIN software












# *************** Build benchmark set: indiviuals with both parental and close relative groups. Assign close relative to parental side using the relatedness of the close relative with the parents.



rel<-as.data.frame(data.table::fread(out_relatives, hea=F))
ped<-as.data.frame(data.table::fread(out_ped, hea=F))
ped<-ped[complete.cases(ped),] # keep only trios



###
# UNTESTED: Remove trios when parental are related
###
remove_related_parents <- function(t){
	fa<-ped$V2[ped$V1==t]
	mo<-ped$V3[ped$V1==t]
	sub1<-d[d$ID1==fa & d$ID2==mo,]
	sub2<-d[d$ID2==fa & d$ID1==mo,]

	if ((dim(sub1)[1]==0) & (dim(sub2)[1]==0)){
		return(ped[ped$V1==t,])
	}
}


x<-mclapply(unique(ped$V1), remove_related_parents, mc.cores=N_CORE)
ped<-cbind(bind_rows(x))

###
#
###


rel<-rel[rel$V1 %in% ped$V1,]
samples<-rel$V1

x<-mclapply(samples, get_side, mc.cores=N_CORE)
SIDE<-cbind(bind_rows(x))
SIDE<-SIDE[!duplicated(SIDE),]
dim(SIDE)

# verify that the same group is not attributed to bother maternal and paternal side
res<-data.frame()
for (t in unique(SIDE$target)){
  sub<-SIDE[SIDE$target==t,]
  g1_side<-unique(as.character(sub$side[sub$group=='G1']))
  g2_side<-unique(as.character(sub$side[sub$group=='G2']))
  
  if (length(g2_side)==0 & length(g1_side)==1){
    res<-rbind(res, sub)
  }
  
  if (length(g1_side)==0 & length(g2_side)==1){
    res<-rbind(res, sub)
  }
  
  else {
    if (length(g1_side)==1 & length(g2_side)==1){
      if (g1_side!=g2_side){
        res<-rbind(res, sub)
      }
    }
  }
  
}
dim(res)



# check cases filtered out
sub<-SIDE[!(SIDE$target %in% res$target),]

res$target_sex<-sex$sex[match(res$target, sex$ID1)]
res$relative_sex<-sex$sex[match(res$relative, sex$ID1)]


d1<-d; d1$id<-paste0(d1$ID1,'_',d1$ID2)
d2<-d; d2$id<-paste0(d2$ID2,'_',d2$ID1)
d<-rbind(d1, d2)
res$id<-paste0(res$target,'_',res$relative)
res$degree<-d$InfType[match(res$id, d$id)]
res2<-res[,c(1,2,3,4,5,6,8)]

odir='data/benchmark/'
system(paste0('mkdir -p ',odir))
write.table(res2, paste0(odir,'Relatives.benchmark.side'), quote=F, col.names=T, row.names=F, sep='\t')

rel<-rel[rel$V1 %in% res$target,]
write.table(rel, paste0(odir,'Relatives.benchmark.group'), quote=F, col.names=F, row.names=F, sep='\t')



