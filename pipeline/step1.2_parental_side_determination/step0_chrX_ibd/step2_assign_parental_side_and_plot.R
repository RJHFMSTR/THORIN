setwd('~/Dropbox/Ressources/Git_repository/THORIN/pipeline/step2_chrX_ibd/')
'%ni%'=Negate('%in%')
library(dplyr)
library(parallel)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(ggdist)
library(ggh4x)
library(scales)
library(ggpubr)

axis_text_size=6
axis_title_size=8
title_size=10
facet_text_size=6
cols<-c("#6699CC","#CC6677","#88CCEE", "#CC6677", "#DDCC77", "#44AA99", "#AA4499", 
        "#999933", "#882255", "#661100", "#888888","#117733", "#332288")







#### detect extreme outliers ####
find_outliers<-function(data, value){
  q1=quantile(data,0.25)
  q3=quantile(data,0.75)
  IQR=q3-q1
  UB=q3+(value*IQR)
  LB=q1-(value*IQR)
  return(data[data>UB | data<LB])
}



# format benchmark groups
rel<-read.table('../step1_surrogate_parents/data/benchmark/Relatives.benchmark.side', hea=T)
rel<-rel[rel$target_sex=='Male',]
rel<-rel[,c(1,3,4)]; rel<-rel[!duplicated(rel),]
rel2<-data.frame()
for (t in unique(rel$target)){
  tmp<-rel[i,]
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
data='data/THORIN/benchmark/KGP.chrX.benchmark.thorin.prob.ibd'
dt<-as.data.frame(data.table::fread(data, hea=T))
dt$side<-rel$side[match(paste0(dt$target,':',dt$UPD),paste0(rel$target,':',rel$group))]
dt<-dt[!is.na(dt$side),]

# add target with paternal or relatives, but no IBD sharing at all
rel<-read.table('../step1_surrogate_parents/data/benchmark/Relatives.benchmark.side', hea=T)
for (i in 1:dim(rel)[1]){
  tmp<-rel[i,]
  t=tmp$target
  if (!(t %in% dt$target)){
    dt<-rbind(dt, data.frame(CHR='chrx', start=0, end=0, Prob='D', length_CM=0, target=t, UPD=tmp$group, side=tmp$side))
  }
}

write.table(dt, 'data/chrX.summary_benchmark.txt', quote=F, col.names=T, row.names=F)



#########################################################################################################
#### Accuracy approach - benchmark #######################################################################
#########################################################################################################


run_benchmark<-function(i){
  
  data='data/chrX.summary_benchmark.txt'
  dt<-as.data.frame(data.table::fread(data, hea=T))
  test<-dt[i,]
  train<-dt[-i,]
  
  dtm<-dt[dt$side=='maternal',]; dtm<-dtm[!is.na(dtm$target),]
  dtp<-dt[dt$side=='paternal',]; dtp<-dtp[!is.na(dtp$target),]
  
  #remove extreme outliers
  olm<-find_outliers(dtm$length_CM, 5)  
  olp<-find_outliers(dtp$length_CM, 5)  
  dtm<-dtm[!(dtm$length_CM %in% olm),]
  dtp<-dtp[!(dtp$length_CM %in% olp),]
  
  # idxm<-sample(1:dim(dtm)[1], dim(dtm)[1]/2)
  # idxp<-sample(1:dim(dtp)[1], dim(dtp)[1]/2)
  
  # testm<-dtm[idxm, ]; dtm<-dtm[-idxm, ]
  # testp<-dtp[idxp, ]; dtp<-dtp[-idxp, ]
  
  dtm$weight<-1/dim(dtm)[1]
  dtp$weight<-1/dim(dtp)[1]
  
  train<-rbind(dtm, dtp)
  # dt<-rbind(testm, testp)
  
  
  
  value2accuracyM<-function(value, data){
    tpm=sum(data$weight[data$side=='maternal' & data$length_CM>=value])
    fpm=sum(data$weight[data$side=='paternal' & data$length_CM>=value])
    am<-tpm/(tpm+fpm)
    if (value>max(data$length_CM)){
      value=max(data$length_CM)
      tpm=sum(data$weight[data$side=='maternal' & data$length_CM>=value])
      fpm=sum(data$weight[data$side=='paternal' & data$length_CM>=value])
      am<-tpm/(tpm+fpm)
    }
    return(am)
  }
  
  value2accuracyP<-function(value, data){
    tpp=sum(data$weight[data$side=='paternal' & data$length_CM<=value])
    fpp=sum(data$weight[data$side=='maternal' & data$length_CM<=value])
    ap<-tpp/(tpp+fpp)
    if (value<min(data$length_CM)){
      value=min(data$length_CM)
      tpp=sum(data$weight[data$side=='paternal' & data$length_CM<=value])
      fpp=sum(data$weight[data$side=='maternal' & data$length_CM<=value])
      ap<-tpp/(tpp+fpp)
    }
    return(ap)
  }
  
  
  test$accuracyM<-unlist(lapply(test$length_CM, value2accuracyM, data=train))
  test$accuracyP<-unlist(lapply(test$length_CM, value2accuracyP, data=train))
  
  
  # keep the maximum accuracy as parental assignment.
  test <- as.data.frame( test %>% rowwise() %>% mutate(accuracy = max(accuracyM, accuracyP), accuracy_type = ifelse(accuracyM >= accuracyP, "accuracyM", "accuracyP")) %>% ungroup())
  
  
  
  # Reverse probability when paternal side prediction is the best choice
  setM<-test[test$accuracy_type=='accuracyM',]
  setP<-test[test$accuracy_type=='accuracyP',]
  
  if (dim(setM)[1]!=0){
    setM$maternal_probability<-setM$accuracy
    setM$maternal_group<-setM$UPD
  }
  if (dim(setP)[1]!=0){
    setP$maternal_probability<-setP$accuracy
    setP$maternal_group<-NA; setP$maternal_group[setP$UPD=='G1']<-'G2'; setP$maternal_group[setP$UPD=='G2']<-'G1'
  }
  dt<-rbind(setM, setP)
  return(dt)
}



# run the benchmark using a Leave-One-Individual-Out approach.
data='data/chrX.summary_benchmark.txt'
dt<-as.data.frame(data.table::fread(data, hea=T))

res<-data.frame()
for (i in 1:dim(dt)[1]){
  print(i)
  res<-rbind(res, run_benchmark(i))
}
res<-res[!duplicated(res),]
res<-res[!is.na(res$target),]
write.table(res, 'data/chrx_accuracy_derived_prediction.benchmark.txt', quote=F, col.names=T, row.names=F, sep='\t')



#########################################################################################################
#### Accuracy approach - call #######################################################################
#########################################################################################################




# training on benchmark set
data='data/chrX.summary_benchmark.txt'
dt<-as.data.frame(data.table::fread(data, hea=T))

dtm<-dt[dt$side=='maternal',]; dtm<-dtm[!is.na(dtm$target),]
dtp<-dt[dt$side=='paternal',]; dtp<-dtp[!is.na(dtp$target),]

olm<-find_outliers(dtm$length_CM, 5)  
olp<-find_outliers(dtp$length_CM, 5)  
dtm<-dtm[!(dtm$length_CM %in% olm),]
dtp<-dtp[!(dtp$length_CM %in% olp),]

dtm$weight<-1/dim(dtm)[1]
dtp$weight<-1/dim(dtp)[1]

train<-rbind(dtm, dtp)


# call on full sample set
data='data/THORIN/call/KGP.chrX.thorin.prob.ibd'
dt<-as.data.frame(data.table::fread(data, hea=T))
dt<-dt[dt$UPD!='U',]
dt<-dt[!(dt$target %in% train$target),]

value2accuracyM<-function(value, data){
  tpm=sum(data$weight[data$side=='maternal' & data$length_CM>=value])
  fpm=sum(data$weight[data$side=='paternal' & data$length_CM>=value])
  am<-tpm/(tpm+fpm)
  if (value>max(data$length_CM)){
    value=max(data$length_CM)
    tpm=sum(data$weight[data$side=='maternal' & data$length_CM>=value])
    fpm=sum(data$weight[data$side=='paternal' & data$length_CM>=value])
    am<-tpm/(tpm+fpm)
  }
  return(am)
}

value2accuracyP<-function(value, data){
  tpp=sum(data$weight[data$side=='paternal' & data$length_CM<=value])
  fpp=sum(data$weight[data$side=='maternal' & data$length_CM<=value])
  ap<-tpp/(tpp+fpp)
  if (value<min(data$length_CM)){
    value=min(data$length_CM)
    tpp=sum(data$weight[data$side=='paternal' & data$length_CM<=value])
    fpp=sum(data$weight[data$side=='maternal' & data$length_CM<=value])
    ap<-tpp/(tpp+fpp)
  }
  return(ap)
}


dt$accuracyM<-unlist(lapply(dt$length_CM, value2accuracyM, data=train))
dt$accuracyP<-unlist(lapply(dt$length_CM, value2accuracyP, data=train))


# keep the maximum accuracy as parental assignment.
dt <- as.data.frame( dt %>% rowwise() %>% mutate(accuracy = max(accuracyM, accuracyP), accuracy_type = ifelse(accuracyM >= accuracyP, "accuracyM", "accuracyP")) %>% ungroup())



# Reverse probability when paternal side prediction is the best choice
setM<-dt[dt$accuracy_type=='accuracyM',]
setP<-dt[dt$accuracy_type=='accuracyP',]

setM$maternal_probability<-setM$accuracy
setM$maternal_group<-setM$UPD

setP$maternal_probability<-setP$accuracy
setP$maternal_group<-NA; setP$maternal_group[setP$UPD=='G1']<-'G2'; setP$maternal_group[setP$UPD=='G2']<-'G1'


dt<-rbind(setM, setP)

write.table(dt, 'data/chrx_accuracy_derived_prediction.txt', quote=F, col.names=T, row.names=F, sep='\t')













#########################################################################################################
#### Plots #######################################################################
#########################################################################################################

### ***
### *** Benchmark
### ***
data='data/chrX.summary_benchmark.txt'
dt<-as.data.frame(data.table::fread(data, hea=T))
colnames(dt)[5]<-'chrx'
subm<-dt[!is.na(dt$side) & dt$side=='maternal', ]
subp<-dt[!is.na(dt$side) & dt$side=='paternal', ]
labm<-paste0('Maternal relative groups (N=',dim(subm)[1],')')
labp<-paste0('Paternal relative groups (N=',dim(subp)[1],')')


gm<-ggplot() +
  geom_boxplot(data=subm, aes(x=chrx, y=-0.15), width=0.2, fill=cols[2], alpha=.3) +
  geom_jitter(data=subm, aes(x=chrx, y=-0.15), alpha=.3, height=0.05) +
  geom_density(data=subm, aes(x=chrx, ..scaled.., fill=labm), alpha=.5) +
  geom_density(data=data.frame(x=0, y=0), aes(x=x, fill=labp), alpha=.5) + # this is only to add paternal legend.
  theme_classic() + theme(legend.position='none') +
  labs(title=paste0(paste0('a. Length of chromosome X IBD with close relative groups - validation cohort (N=',dim(dt[!is.na(dt$side),]))[1],')'),
       x='Chromosome X IBD length [cM]', y='Density', col='', fill='') +
  theme( axis.title.y = element_text(margin = margin(r = 3)),
         axis.title.x = element_blank(), axis.text = element_text(size=axis_text_size),
         axis.text.x = element_blank(), axis.ticks.x = element_blank(),
         axis.title = element_text(size=axis_title_size),
         plot.title = element_text(size=title_size),
         plot.title.position = "plot", legend.position = 'none') +
  scale_x_continuous(limits=c(0, max(dt$chrx))) + scale_fill_manual(values=c(cols[2], cols[3])) + scale_y_continuous(breaks=seq(0,1,0.25), labels=seq(0,1,0.25), limits=c(-0.25,1))

gp<-ggplot(subp, aes(x=chrx)) +
  geom_boxplot(alpha=.3, aes(x=chrx, y=-0.15), fill=cols[3], width=0.2) +
  geom_jitter(alpha=.3, aes(x=chrx, y=-0.15), height=0.05) +
  geom_density(aes(x=chrx,..scaled.., fill=labp), alpha=.5) +
  geom_density(data=data.frame(x=0, y=0), aes(x=x, fill=labm), alpha=.5) + # this is only to add maternal legend.
  theme_classic() + theme(legend.position='none') +
  labs(x='Chromosome X IBD length [cM]', y='Density', fill='', col='') +
  theme( axis.title.y = element_text(margin = margin(r = 3)),
         axis.text = element_text(size=axis_text_size), axis.title = element_text(size=axis_title_size),
         plot.title = element_text(size=title_size),
         plot.title.position = "plot",
         legend.key.size = unit(0.4, 'cm'),
         legend.text = element_text(size=axis_text_size),
         legend.key.spacing.y = unit(0.1, 'cm'),
         legend.position = c(0.8,0.6)) +
  scale_color_manual(values=c(cols[2], cols[3])) +
  scale_x_continuous(limits=c(0, max(dt$chrx))) + scale_fill_manual(values=c(cols[2], cols[3])) + scale_y_continuous(breaks=seq(0,1,0.25), labels=seq(0,1,0.25), limits=c(-0.25,1))

a1<-plot_grid(gm, NULL, gp, ncol=1, rel_heights = c(1,-0.05,1))







### ***
### *** Call probabilities
### ***



dt<-read.table("data/chrx_accuracy_derived_prediction.txt", hea=T)
colnames(dt)[5]<-'chrx'
dt<-dt[!is.na(dt$target),]
data<-data.frame(chrx=rep(dt$chrx, 2),
                 accuracy=c(dt$accuracyM, dt$accuracyP), 
                 type=c(rep('maternal', length(dt$accuracyM)), rep('paternal', length(dt$accuracyP))),
                 accuracy_type=c(dt$accuracy_type, dt$accuracy_type))
data$linecol<-paste0(data$type, data$accuracy_type)
data$linecol<-factor(data$linecol, levels=c("maternalaccuracyP", "paternalaccuracyM","maternalaccuracyM", "paternalaccuracyP"))


a2_1<-ggplot(data, aes(x=chrx, y=accuracy, col=linecol, linetype=linecol, linewidth=linecol)) + geom_line() +
  lims(x=c(0,8)) +
  theme_classic() +
  theme(axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size=axis_text_size),
        axis.title = element_text(size=axis_title_size),
        plot.title = element_text(size=title_size),
        plot.title.position = "plot", legend.position = c(0.6, 0.25),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=axis_text_size),
        legend.key.spacing.y = unit(0.1, 'cm'),
        legend.key.spacing.x = unit(0.1, 'cm'),
        legend.title=element_text(size=axis_text_size)) +
  scale_color_manual(values=c(cols[2], cols[3],cols[2], cols[3]), labels=c('/','/','(un)-selected maternal probability','(un)-selected paternal probability')) +
  scale_linetype_manual(values=c('dotted','dotted','solid','solid'), labels=c('/','/','(un)-selected maternal probability','(un)-selected paternal probability')) +
  scale_linewidth_manual(values=c(0.8,0.8,1.5,1.5), labels=c('/','/','(un)-selected maternal probability','(un)-selected paternal probability')) +
  guides(color=guide_legend(nrow=2)) +
  labs(title='b. Parental-side assignment probability', x='Chromosome X IBD length [cM]', y='Probability',
       col='', linetype='', linewidth='')








### ***
### *** Call probabilities accuracy
### ***



dt<-read.table("data/chrx_accuracy_derived_prediction.benchmark.txt", hea=T)
colnames(dt)[5]<-'chrx'
data<-dt[!is.na(dt$side),]
tmp<-data.frame()
step=0.005
for (i in c(seq(min(data$maternal_probability)+step,1,step),1)){
  tpM<-dim(data[data$maternal_probability>=(i-step) & data$side=='maternal' & data$accuracy_type=='accuracyM',])[1]
  tpP<-dim(data[data$maternal_probability>=(i-step) & data$side=='paternal' & data$accuracy_type=='accuracyP',])[1]
  fpM<-dim(data[data$maternal_probability>=(i-step) & data$side=='maternal' & data$accuracy_type=='accuracyP',])[1]
  fpP<-dim(data[data$maternal_probability>=(i-step) & data$side=='paternal' & data$accuracy_type=='accuracyM',])[1]
  tp=tpM+tpP
  fp=fpM+fpP
  tmp<-rbind(tmp, data.frame(cut=i, accuracy=tp/(tp+fp)))
}

a3_1<-ggplot(tmp, aes(x=cut, y=accuracy)) + geom_line(alpha=.8, , linewidth=1) + geom_point() +
  theme_classic() +
  theme(axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size=axis_text_size),
        axis.title = element_text(size=axis_title_size),
        plot.title = element_text(size=title_size),
        plot.title.position = "plot", legend.position = c(0.75, 0.65),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=axis_text_size),
        legend.key.spacing.y = unit(0.1, 'cm'),
        legend.title=element_text(axis_title_size)) +
  labs(title='c. Probability accuracy', x='Probability', y='Accuracy',
       col='') + lims(x=c(0.5, 1), y=c(0.5, 1))





### ***
### *** Call probabilities frequency
### ***
dt<-read.table("chrx_accuracy_derived_prediction.txt", hea=T)
colnames(dt)[5]<-'chrx'
dt<-dt[!is.na(dt$target),]
step=0.0025
d_data<-data.frame(bin=c(0.5,0.5, min(dt$accuracy)-step, min(dt$accuracy)-step), n=c(0,0), type=c('Paternal-side', 'Maternal-side'))
for (bin in seq(min(dt$accuracy), max(dt$accuracy)+step, step)){
  nm<-dim(dt[dt$accuracy>bin & dt$accuracy<(bin+step) & dt$accuracy_type=='accuracyM',])[1]
  np<-dim(dt[dt$accuracy>bin & dt$accuracy<(bin+step) & dt$accuracy_type=='accuracyP',])[1]
  d_data<-rbind(d_data, data.frame(bin=bin, n=nm, type='Maternal-side'))
  d_data<-rbind(d_data, data.frame(bin=bin, n=np, type='Paternal-side'))
}

a4<-ggplot(data=d_data, aes(x=bin, y=n)) + geom_area(aes(fill=type), col='black', alpha=.5) + facet_grid(type~.) +
  theme_classic() +
  labs(title=paste0('d. Parental side predictions from chromosome X IBD for N=',comma_format()(dim(dt)[1]),' male individuals'), x='Probability', y='Number of individuals', fill='We predict') +
  theme(axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size=axis_text_size),
        axis.title = element_text(size=axis_title_size),
        plot.title = element_text(size=title_size),
        plot.title.position = "plot", legend.position = c(0.85, 0.3),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=axis_text_size),
        legend.title = element_text(size=axis_text_size),
        legend.key.spacing.y = unit(0.15, 'cm')) +
  scale_fill_manual(values=c(cols[2], cols[3]), labels=c(paste0('Maternal side (N=',comma_format()(dim(dt[dt$accuracy_type=='accuracyM',])[1]),')'),
                                                         paste0('Paternal side (N=',comma_format()(dim(dt[dt$accuracy_type=='accuracyP',])[1]),')')))







a<-plot_grid(a1, NULL, ggarrange(a2_1, a3_1, ncol=2, nrow=1, common.legend = TRUE, legend="bottom"),NULL, a4, ncol=1, rel_heights = c(1,0.1,0.85,0.1,1))
ggsave('data/chrx_accuracy_derived_prediction.jpg', a,width=8, height=10)



