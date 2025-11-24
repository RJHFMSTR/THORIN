# Create a data with all phenotypes used to determine mate pairs. Use the following fields
#Phenotype_Name	Phenotype_Code
#home_location_east_1km 22702
#townsend_deprivation_index 22189
#home_location_north_1km 22704
#Average_total_household_income_before_tax 738
#smokers_in_household 1259
#number_in_household 709
#lengt_of_time_at_current_address 699
#How_are_people_in_household_related_to_participant 6141

phenotype_file="Out/Mate_phenotypes.txt"
relatedness_file="/data/FAC/FBM/DBC/zkutalik/default_sensitive/uk_biobank/useful_data/ukb_relatedness.up_to_degree_4.txt.gz" # the KING file provided by UKBB
age_file="/data/FAC/FBM/DBC/zkutalik/default_sensitive/uk_biobank/useful_data/ukb_age.txt.gz" # col1=IID; col2=age
wb_file="/data/FAC/FBM/DBC/zkutalik/default_sensitive/uk_biobank/useful_data/ukb_caucasians_for_regenie.txt"

d<-as.data.frame(data.table::fread(phenotype_file, hea=T))

# select only indiv reporting living with mate.
d<-d[d$living_with_husband_or_wife_or_partner==1,]

# remove indiv with missing data.
d<-d[!is.na(d$Average_total_household_income_before_tax),]
d<-d[!is.na(d$home_location_east_1km),]
d<-d[!is.na(d$townsend_deprivation_index),]
d<-d[!is.na(d$home_location_north_1km),]
d<-d[!is.na(d$Average_total_household_income_before_tax),]
d<-d[!is.na(d$smokers_in_household),]
d<-d[!is.na(d$number_in_household),]
d<-d[!is.na(d$lengt_of_time_at_current_address),]
d<-d[!is.na(d$living_with_husband_or_wife_or_partner),]
d<-d[!is.na(d$genetic_sex),]

# remove non-respondant
d<-d[d$Average_total_household_income_before_tax!=1 & d$Average_total_household_income_before_tax!=-3,]
d<-d[d$smokers_in_household!=-3,]
d<-d[d$number_in_household>0,]
d<-d[d$lengt_of_time_at_current_address!=-1 & d$lengt_of_time_at_current_address!=-3,]

# create unique identifier for all phenotypes
d$pheno_comb<-paste0(d$home_location_east_1km,'_',d$home_location_north_1km,'_',d$townsend_deprivation_index,'_',d$smokers_in_household,'_',d$number_in_household,'_',d$lengt_of_time_at_current_address,'_',d$Average_total_household_income_before_tax)

# get individual with same identifier
cluster_pairs<-function(pheno) {
  individuals <- d$IID[d$pheno_comb == pheno]
  if (length(individuals) > 1) {
    pairs <- t(combn(individuals, 2))
    return(data.frame(ID1 = pairs[,1], ID2 = pairs[,2], pheno_comb = pheno, stringsAsFactors = FALSE))
  }
}
library(parallel); library(dplyr)
x<-mclapply(unique(d$pheno_comb), cluster_pairs, mc.cores=36)
result<-rbind(bind_rows(x))

# add genetic sex
result$sex_id1<-d$genetic_sex[match(result$ID1, d$IID)]
result$sex_id2<-d$genetic_sex[match(result$ID2, d$IID)]

# remove non-opposite mate pairs
res<-result[result$sex_id1 != result$sex_id2,]

# add genetic relatedness
rel<-as.data.frame(data.table::fread(relatedness_file, hea=T))
rel_pairs<-unique(c(paste0(rel$ID1,'_',rel$ID2), paste0(rel$ID2,'_',rel$ID1)))

# remove related pairs
res$pair_id<-paste0(res$ID1, '_', res$ID2)
res<-res[!(res$pair_id %in% rel_pairs),]

# add age
age<-as.data.frame(data.table::fread(age_file, hea=T))
res$age_id1<-age$age[match(res$ID1, age$IID)]
res$age_id2<-age$age[match(res$ID2, age$IID)]
res<-res[complete.cases(res),]

# remove pairs with >=10 years diff
res$diff_age<-abs(res$age_id1-res$age_id2)
res<-res[res$diff_age<=10,]

# remove non white British (field 22006)
cau<-as.data.frame(data.table::fread(wb_file, hea=F))
res<-res[res$ID1 %in% cau$V1,]
res<-res[res$ID2 %in% cau$V1,]

# remove pairs for which one of the mate appear in >1 pair
t<-table(c(res$ID1, res$ID2)); n<-names(t)[unname(t)>1]
res<-res[!(res$ID1 %in% n) & !(res$ID1 %in% n),]

write.table(res, 'Mate_pairs.txt', quote=F, col.names=T, row.names=F, sep='\t')

