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




#########################################################################################################
#### Accuracy approach - benchmark #######################################################################
#########################################################################################################





run_benchmark<-function(i, degree, data){
  
  dt<-data
  dt<-dt[dt$degree==degree, ]
  test<-dt[i,]
  train<-dt[-i,]
  
  dtm<-dt[dt$side=='maternal',]; dtm<-dtm[!is.na(dtm$target),]
  dtp<-dt[dt$side=='paternal',]; dtp<-dtp[!is.na(dtp$target),]
  
  #remove extreme outliers
  olm<-find_outliers(dtm$mtdna, 5)  
  olp<-find_outliers(dtp$mtdna, 5)  
  dtm<-dtm[!(dtm$mtdna %in% olm),]
  dtp<-dtp[!(dtp$mtdna %in% olp),]
  
  dtm$weight<-1/dim(dtm)[1]
  dtp$weight<-1/dim(dtp)[1]
  
  train<-rbind(dtm, dtp)

  value2accuracyM<-function(value, data){
    tpm=sum(data$weight[data$side=='maternal' & data$mtdna>=value])
    fpm=sum(data$weight[data$side=='paternal' & data$mtdna>=value])
    am<-tpm/(tpm+fpm)
    if (value>max(data$mtdna)){
      value=max(data$mtdna)
      tpm=sum(data$weight[data$side=='maternal' & data$mtdna>=value])
      fpm=sum(data$weight[data$side=='paternal' & data$mtdna>=value])
      am<-tpm/(tpm+fpm)
    }
    return(am)
  }
  
  value2accuracyP<-function(value, data){
    tpp=sum(data$weight[data$side=='paternal' & data$mtdna<=value])
    fpp=sum(data$weight[data$side=='maternal' & data$mtdna<=value])
    ap<-tpp/(tpp+fpp)
    if (value<min(data$mtdna)){
      value=min(data$mtdna)
      tpp=sum(data$weight[data$side=='paternal' & data$mtdna<=value])
      fpp=sum(data$weight[data$side=='maternal' & data$mtdna<=value])
      ap<-tpp/(tpp+fpp)
    }
    return(ap)
  }
  
  
  test$accuracyM<-unlist(lapply(test$mtdna, value2accuracyM, data=train))
  test$accuracyP<-unlist(lapply(test$mtdna, value2accuracyP, data=train))
  
  
  # keep the maximum accuracy as parental assignment.
  test <- as.data.frame( test %>% rowwise() %>% mutate(accuracy = max(accuracyM, accuracyP), accuracy_type = ifelse(accuracyM >= accuracyP, "accuracyM", "accuracyP")) %>% ungroup())
  
  
  
  # Reverse probability when paternal side prediction is the best choice
  setM<-test[test$accuracy_type=='accuracyM',]
  setP<-test[test$accuracy_type=='accuracyP',]
  
  if (dim(setM)[1]!=0){
    setM$maternal_probability<-setM$accuracy
    setM$maternal_group<-setM$group
  }
  if (dim(setP)[1]!=0){
    setP$maternal_probability<-setP$accuracy
    setP$maternal_group<-NA; setP$maternal_group[setP$group=='G1']<-'G2'; setP$maternal_group[setP$group=='G2']<-'G1'
  }
  dt<-rbind(setM, setP)
  return(dt)
}





###
### run the benchmark using a Leave-One-Individual-Out approach.
###

data='data/mtDNA_MVS.txt'
dt<-as.data.frame(data.table::fread(data, hea=T))
dt<-dt[!is.na(dt$side),]
dt_benchmark<-dt


# !! note: for the example on the KGP data, we don't have enough individuals from each group since the sample size of the entire cohort is significantly smaller than biobanks.
## !! Hence, we also use random individual that we add in the validation cohort just for the example.
data='data/mtDNA_MVS.txt'
dt<-as.data.frame(data.table::fread(data, hea=T))
dt1<-dt[!is.na(dt$side),] # correct validation cohort
dt_sim2<-dt[dt$degree=='2nd',]; dt_sim2<-dt_sim2[1:10,] # random individuals of the 2nd degree
dt_sim3<-dt[dt$degree=='3rd',]; dt_sim3<-dt_sim3[1:10,] # random individuals of the 3rd degree
dt_sim4<-dt[dt$degree=='4th',]; dt_sim4<-dt_sim4[1:10,] # random individuals of the 4th degree
dt_benchmark<-rbind(dt1, dt_sim2, dt_sim3, dt_sim4)
dt_benchmark$side[dt_benchmark$mtdna>0.98]<-'maternal'
dt_benchmark$side[dt_benchmark$mtdna<=0.98]<-sample(c('paternal', 'maternal'), dim(dt_benchmark[dt_benchmark$mtdna<=0.98,])[1], replace=T)
## !! End of modification




# 2nd degree
dt2<-data.frame()
for (i in 1:dim(dt_benchmark[dt_benchmark$degree=='2nd',])[1]){
  print(i)
  dt2<-rbind(dt2, run_benchmark(i, '2nd', data=dt_benchmark))
}


# 3rd degree
dt3<-data.frame()
for (i in 1:dim(dt_benchmark[dt_benchmark$degree=='3rd',])[1]){
  print(i)
  dt3<-rbind(dt3, run_benchmark(i, '3rd', data=dt_benchmark))
}


# 4th degree
dt4<-data.frame()
for (i in 1:dim(dt_benchmark[dt_benchmark$degree=='4th',])[1]){
  print(i)
  dt4<-rbind(dt4, run_benchmark(i, '4th', data=dt_benchmark))
}


write.table(rbind(dt2, dt3, dt4), 'data/mtdna_accuracy_derived_prediction.benchmark_not_uniq.txt', quote=F, col.names=T, row.names=F, sep='\t')

















#########################################################################################################
#### Accuracy approach - call #######################################################################
#########################################################################################################







derive_probs<-function(data, degree){
  
  dt<-data
  dt<-dt[dt$degree==degree, ]
  train<-dt[!is.na(dt$side),]
  
  # !! note: for the example on the KGP data, we don't have enough individuals from each group since the sample size of the entire cohort is significantly smaller than biobanks.
  ## !! Hence, we also use random individual that we add in the validation cohort just for the example.
  dt1<-dt[!is.na(dt$side),] # correct validation cohort
  dt_sim<-dt[dt$degree==degree,]; dt_sim<-dt_sim[1:10,] # random individuals 
  train<-rbind(dt1, dt_sim)
  train$side[train$mtdna>0.98]<-'maternal'
  train$side[train$mtdna<=0.98]<-sample(c('paternal', 'maternal'), dim(train[train$mtdna<=0.98,])[1], replace=T)
  ## !! End of modification
  
  
  
  test<-dt[!(dt$target %in% train$target),]
  
  dtm<-train[train$side=='maternal',]; dtm<-dtm[!is.na(dtm$target),]
  dtp<-train[train$side=='paternal',]; dtp<-dtp[!is.na(dtp$target),]
  
  #remove extreme outliers
  olm<-find_outliers(dtm$mtdna, 5)  
  olp<-find_outliers(dtp$mtdna, 5)  
  dtm<-dtm[!(dtm$mtdna %in% olm),]
  dtp<-dtp[!(dtp$mtdna %in% olp),]
  
  dtm$weight<-1/dim(dtm)[1]
  dtp$weight<-1/dim(dtp)[1]
  
  train<-rbind(dtm, dtp)
  
  value2accuracyM<-function(value, data){
    tpm=sum(data$weight[data$side=='maternal' & data$mtdna>=value])
    fpm=sum(data$weight[data$side=='paternal' & data$mtdna>=value])
    am<-tpm/(tpm+fpm)
    if (value>max(data$mtdna)){
      value=max(data$mtdna)
      tpm=sum(data$weight[data$side=='maternal' & data$mtdna>=value])
      fpm=sum(data$weight[data$side=='paternal' & data$mtdna>=value])
      am<-tpm/(tpm+fpm)
    }
    return(am)
  }
  
  value2accuracyP<-function(value, data){
    tpp=sum(data$weight[data$side=='paternal' & data$mtdna<=value])
    fpp=sum(data$weight[data$side=='maternal' & data$mtdna<=value])
    ap<-tpp/(tpp+fpp)
    if (value<min(data$mtdna)){
      value=min(data$mtdna)
      tpp=sum(data$weight[data$side=='paternal' & data$mtdna<=value])
      fpp=sum(data$weight[data$side=='maternal' & data$mtdna<=value])
      ap<-tpp/(tpp+fpp)
    }
    return(ap)
  }
  
  
  test$accuracyM<-unlist(lapply(test$mtdna, value2accuracyM, data=train))
  test$accuracyP<-unlist(lapply(test$mtdna, value2accuracyP, data=train))
  
  
  # keep the maximum accuracy as parental assignment.
  test <- as.data.frame( test %>% rowwise() %>% mutate(accuracy = max(accuracyM, accuracyP), accuracy_type = ifelse(accuracyM >= accuracyP, "accuracyM", "accuracyP")) %>% ungroup())
  
  
  
  # Reverse probability when paternal side prediction is the best choice
  setM<-test[test$accuracy_type=='accuracyM',]
  setP<-test[test$accuracy_type=='accuracyP',]
  
  if (dim(setM)[1]!=0){
    setM$maternal_probability<-setM$accuracy
    setM$maternal_group<-setM$group
  }
  if (dim(setP)[1]!=0){
    setP$maternal_probability<-setP$accuracy
    setP$maternal_group<-NA; setP$maternal_group[setP$group=='G1']<-'G2'; setP$maternal_group[setP$group=='G2']<-'G1'
  }
  dt<-rbind(setM, setP)
  return(dt)
}


data<-as.data.frame(data.table::fread('data/mtDNA_MVS.txt', hea=T))
dt2<-derive_probs(data, '2nd')
dt3<-derive_probs(data, '3rd')
dt4<-derive_probs(data, '4th')

write.table(rbind(dt2, dt3, dt4), 'data/mtdna_accuracy_derived_prediction.not_uniq.txt', quote=F, col.names=T, row.names=F, sep='\t')










#########################################################################################################
#### Plots #######################################################################
#########################################################################################################


plot_data<-function(degree){
  
  
  ### ***
  ### *** Benchmark
  ### ***
  
  
  data='data/mtDNA_MVS.txt'
  dt<-as.data.frame(data.table::fread(data, hea=T))
  dt<-dt[dt$degree==degree,]
  
  
  
  # !! note: for the example on the KGP data, we don't have enough individuals from each group since the sample size of the entire cohort is significantly smaller than biobanks.
  ## !! Hence, we also use random individual that we add in the validation cohort just for the example.
  dt1<-dt[!is.na(dt$side),] # correct validation cohort
  dt_sim<-dt[dt$degree==degree,]; dt_sim<-dt_sim[1:10,] # random individuals 
  train<-rbind(dt1, dt_sim)
  train$side[train$mtdna>0.98]<-'maternal'
  train$side[train$mtdna<=0.98]<-sample(c('paternal', 'maternal'), dim(train[train$mtdna<=0.98,])[1], replace=T)
  dt<-train
  ## !! End of modification
  
  
  
  
  subm<-dt[!is.na(dt$side) & dt$side=='maternal' & !is.na(dt$mtdna), ]
  subp<-dt[!is.na(dt$side) & dt$side=='paternal' & !is.na(dt$mtdna), ]
  
  #remove extreme outliers
  olm<-find_outliers(subm$mtdna,5)  
  olp<-find_outliers(subp$mtdna, 5)  
  subm<-subm[!(subm$mtdna %in% olm),]
  subp<-subp[!(subp$mtdna %in% olp),]
  
  
  labm<-paste0('Maternal relative pairs (N=',comma_format()(dim(subm)[1]),')')
  labp<-paste0('Paternal relative pairs (N=',(dim(subp)[1]),')')
  
  
  gm<-ggplot() +
    geom_boxplot(data=subm, aes(x=mtdna, y=-0.15), width=0.2, fill=cols[2], alpha=.3) +
    geom_jitter(data=subm, aes(x=mtdna, y=-0.15), alpha=.3, height = 0.05) +
    geom_density(data=subm, aes(x=mtdna, ..scaled.., fill=labm), alpha=.5) +
    geom_density(data=data.frame(x=0, y=0), aes(x=x, fill=labp), alpha=.5) + # this is only to add paternal legend.
    theme_classic() + theme(legend.position='none') +
    labs(title=paste0(paste0('a. mtDNA variant sharing proportion with ', degree,' degree relatives - validation cohort (N=',comma_format()(dim(dt[!is.na(dt$side),]))[1]),')'),
         x='mtDNA variant sharing proportion', y='Density', col='', fill='') +
    theme( axis.title.y = element_text(margin = margin(r = 3)),
           axis.title.x = element_blank(), axis.text = element_text(size=axis_text_size),
           axis.text.x = element_blank(), axis.ticks.x = element_blank(),
           axis.title = element_text(size=axis_title_size),
           plot.title = element_text(size=title_size),
           plot.title.position = "plot", legend.position = 'none') +
    scale_x_continuous(limits=c(0, 1)) + scale_fill_manual(values=c(cols[2], cols[3])) + scale_y_continuous(breaks=seq(0,1,0.25), labels=seq(0,1,0.25), limits=c(-0.25,1))
  
  gp<-ggplot(subp, aes(x=mtdna)) +
    geom_boxplot(alpha=.3, aes(x=mtdna, y=-0.15), fill=cols[3], width=0.2) +
    geom_jitter(alpha=.3, aes(x=mtdna, y=-0.15), height=0.05) +
    geom_density(aes(x=mtdna,..scaled.., fill=labp), alpha=.5) +
    geom_density(data=data.frame(x=0, y=0), aes(x=x, fill=labm), alpha=.5) + # this is only to add maternal legend.
    theme_classic() + theme(legend.position='none') +
    labs(x='mtDNA variant sharing proportion', y='Density', fill='', col='') +
    theme( axis.title.y = element_text(margin = margin(r = 3)),
           axis.text = element_text(size=axis_text_size), axis.title = element_text(size=axis_title_size),
           plot.title = element_text(size=title_size),
           plot.title.position = "plot",
           legend.key.size = unit(0.4, 'cm'),
           legend.text = element_text(size=axis_text_size),
           legend.key.spacing.y = unit(0.1, 'cm'),
           legend.position = c(0.8,0.6)) +
    scale_color_manual(values=c(cols[2], cols[3])) +
    scale_x_continuous(limits=c(0, 1)) + scale_fill_manual(values=c(cols[2], cols[3])) + scale_y_continuous(breaks=seq(0,1,0.25), labels=seq(0,1,0.25), limits=c(-0.25,1))
  
  a1<-plot_grid(gm, NULL, gp, ncol=1, rel_heights = c(1,-0.05,1))
  
  
  
  
  
  
  
  ### ***
  ### *** Call probabilities
  ### ***
  
  
  
  
  
  
  dt<-read.table("data/mtdna_accuracy_derived_prediction.not_uniq.txt", hea=T)
  dt<-dt[dt$degree==degree,]
  dt<-dt[!is.na(dt$target),]
  data<-data.frame(mtdna=rep(dt$mtdna, 2),
                   accuracy=c(dt$accuracyM, dt$accuracyP), 
                   type=c(rep('maternal', length(dt$accuracyM)), rep('paternal', length(dt$accuracyP))),
                   accuracy_type=c(dt$accuracy_type, dt$accuracy_type))
  data$linecol<-paste0(data$type, data$accuracy_type)
  data$linecol<-factor(data$linecol, levels=c("maternalaccuracyP", "paternalaccuracyM","maternalaccuracyM", "paternalaccuracyP"))
  
  
  a2_1<-ggplot(data, aes(x=mtdna, y=accuracy, col=linecol, linetype=linecol, linewidth=linecol)) + geom_line() +
    lims(x=c(0,1)) +
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
    labs(title='b. Parental-side assignment probability', x='mtDNA variant sharing proportion', y='Probability',
         col='', linetype='', linewidth='')
  
  
  
  
  
  
  
  
  ### ***
  ### *** Call probabilities accuracy
  ### ***
  
  
  
  
  dt<-read.table("data/mtdna_accuracy_derived_prediction.benchmark_not_uniq.txt", hea=T)
  dt<-dt[dt$degree==degree,]; dt<-dt[!is.na(dt$target),]
  data<-dt[!is.na(dt$side),]
  tmp<-data.frame()
  step=0.001
  for (i in seq(min(data$maternal_probability),max(data$maternal_probability),step)){
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
  dt<-read.table("data/mtdna_accuracy_derived_prediction.not_uniq.txt", hea=T)
  dt<-dt[dt$degree==degree,]; dt<-dt[!is.na(dt$target),]
  step=0.0025
  d_data<-data.frame(bin=c(0.5,0.5, min(dt$accuracy)-0.005, min(dt$accuracy)-0.005), n=c(0,0), type=c('Paternal-side', 'Maternal-side'))
  for (bin in seq(min(dt$accuracy), max(dt$accuracy)+step, step)){
    nm<-dim(dt[dt$accuracy>=bin & dt$accuracy<(bin+step) & dt$accuracy_type=='accuracyM',])[1]
    np<-dim(dt[dt$accuracy>=bin & dt$accuracy<(bin+step) & dt$accuracy_type=='accuracyP',])[1]
    d_data<-rbind(d_data, data.frame(bin=bin, n=nm, type='Maternal-side'))
    d_data<-rbind(d_data, data.frame(bin=bin, n=np, type='Paternal-side'))
  }
  
  a4<-ggplot(data=d_data, aes(x=bin, y=n)) + geom_area(aes(fill=type), col='black', alpha=.5) + facet_grid(type~.) +
    theme_classic() +
    labs(title=paste0('d. Parental side predictions from mtDNA for N=',comma_format()(dim(dt)[1]),' target-',degree,' degree relative pairs'), x='Probability', y='Number of individuals', fill='We predict') +
    theme(axis.title.y = element_text(margin = margin(r = 10)),
          axis.text = element_text(size=axis_text_size),
          axis.title = element_text(size=axis_title_size),
          plot.title = element_text(size=title_size),
          plot.title.position = "plot", legend.position = c(0.8, 0.3),
          legend.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          legend.key.size = unit(0.4, 'cm'),
          legend.text = element_text(size=axis_text_size),
          legend.title = element_text(size=axis_text_size),
          legend.key.spacing.y = unit(0.15, 'cm')) +
    scale_fill_manual(values=c(cols[2], cols[3]), labels=c(paste0('Maternal side (N=',comma_format()(dim(dt[dt$accuracy_type=='accuracyM',])[1]),' pairs)'),
                                                           paste0('Paternal side (N=',comma_format()(dim(dt[dt$accuracy_type=='accuracyP',])[1]),' pairs)')))
  
  
  
  
  plot_grid(a1, NULL, ggarrange(a2_1, a3_1, ncol=2, nrow=1, common.legend = TRUE, legend="bottom"),NULL, a4, ncol=1, rel_heights = c(1,0.1,0.85,0.1,1))
  
  system('mkdir -p Plots/')
  ggsave(paste0('Plots/mtdna_accuracy_derived_prediction.degree_',degree,'.jpg'), width=8, height=10)
  
  
}



plot_data('2nd')
plot_data('3rd')
plot_data('4th')








#########################################################################################################
#### Uniq: keep only the best prediction per target per degree #######################################################################
#########################################################################################################



data2best<-function(target, data, degree){
  tmp<-data[data$target==target & data$degree==degree,]; tmp<-tmp[!is.na(tmp$target),]
  accuracy<-tmp$accuracy
  accuracy_type<-tmp$accuracy_type
  idx<-which(accuracy==max(accuracy, na.rm=T))
  t<-table(accuracy_type[idx])
  accuracy_type<-names(t)[unname(t)==max(t)]
  if (length(accuracy_type)>1){accuracy_type<-"accuracyM"} # we do better at defining maternal side.
  maternal_group<-tmp$maternal_group[tmp$accuracy_type==accuracy_type & tmp$accuracy==max(accuracy, na.rm=T)]
  t<-table(maternal_group)
  group<-names(t)[unname(t)==max(t)]
  if (length(group)==1){
    tmp<-tmp[tmp$maternal_group==group,]
    return(data.frame(target=target, degree=degree, maternal_group=group, n_rel=length(tmp$maternal_group[tmp$accuracy_type==accuracy_type & tmp$accuracy==max(accuracy, na.rm=T)]),
                      mtdna_mean=mean(tmp$mtdna[tmp$accuracy_type==accuracy_type & tmp$accuracy==max(accuracy, na.rm=T)]),
                      accuracyM=mean(tmp$accuracyM[tmp$accuracy_type==accuracy_type & tmp$accuracy==max(accuracy, na.rm=T)]),
                      accuracyP=mean(tmp$accuracyP[tmp$accuracy_type==accuracy_type & tmp$accuracy==max(accuracy, na.rm=T)]),
                      accuracy=max(accuracy, na.rm=T), accuracy_type=accuracy_type,
                      maternal_probability=max(tmp$maternal_probability[tmp$accuracy_type==accuracy_type & tmp$accuracy==max(accuracy, na.rm=T)])))
  }
}


x<-mclapply(unique(dt2$target), data2best, data=dt2, degree='2nd', mc.cores=15)
xx2<-rbind(bind_rows(x))
xx2<-xx2[!duplicated(xx2),]
dim(dt2)[1]; dim(xx2); length(unique(dt2$target)); length(unique(xx2$target))




x<-mclapply(unique(dt3$target), data2best, data=dt3, degree='3rd', mc.cores=15)
xx3<-rbind(bind_rows(x))
dim(dt3)[1]; dim(xx3); length(unique(dt3$target)); length(unique(xx3$target))


x<-mclapply(unique(dt4$target), data2best, data=dt4, degree='4th', mc.cores=15)
xx4<-rbind(bind_rows(x))
dim(dt4)[1]; dim(xx4); length(unique(dt4$target)); length(unique(xx4$target))



dt<-rbind(xx2, xx3, xx4)
dt<-dt[,c(1,2,3,10)]
write.table(dt, 'data/mtdna_accuracy_derived_prediction.uniq_per_degree.txt', quote=F, col.names=T, row.names=F, sep='\t')





######################################################
######### PLOT #######################################
######################################################

## Probability frequency binned and colored by selection
dt<-read.table('data/mtdna_accuracy_derived_prediction.uniq_per_degree.txt', hea=T)
dt$bin<-NA
step=0.01
for (bin in seq(0.5, 1, step)){
  dt$bin[dt$maternal_prob>bin & dt$maternal_prob<=(bin+step)]<-paste0(bin,'-',bin+step)
}

dt<-dt[!is.na(dt$bin),]

dt$degree<-factor(dt$degree, levels=c('2nd','3rd','4th'), labels=c(paste0('2nd degree (N=',comma_format()(dim(dt[dt$degree=='2nd',])[1]),')'),
                                                                   paste0('3rd degree (N=',comma_format()(dim(dt[dt$degree=='3rd',])[1]),')'),
                                                                   paste0('4th degree (N=',comma_format()(dim(dt[dt$degree=='4th',])[1]),')')))
g1<-ggplot(dt, aes(x=factor(bin), fill=degree)) + geom_bar() +
  theme_classic() +
  scale_x_discrete(guide = guide_axis(angle = 45)) +
  labs(title=paste0('a. mtDNA probability distribution for ',comma_format()(dim(dt)[1]),' target-relative pairs'),
       x='Probability bin', y='Number of target-relative pairs', fill='Relatedness degree: ') +
  theme(axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size=axis_text_size),
        axis.title = element_text(size=axis_title_size),
        plot.title = element_text(size=title_size),
        plot.title.position = "plot", legend.position = 'bottom',
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=axis_text_size),
        legend.title = element_text(size=axis_text_size),
        legend.key.spacing.y = unit(0.15, 'cm'),
        legend.key.spacing.x = unit(0.55, 'cm'),
        legend.direction="horizontal") +
  scale_fill_manual(values=c(cols[7], cols[6], cols[5]))






# add: number of relative supporting the group assignment.
t<-data.frame(t(table(dt$n_rel, dt$degree)))

t$Freq_adjusted <- t$Freq + 1
g2<-ggplot(t, aes(x = as.numeric(as.character(Var2)), y = Freq_adjusted, fill = Var1, group = Var1)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_y_log10(breaks=1+c(0, 1,10,10**2,10**3, 10**4, 10**5),
                labels=c('0','1','10','10e1','10e2','10e3','10e4')) +
  labs(title='b. Number of relatives supporting the same parental assignment',
       x = "Number of relatives supporting the parental assignment", y = "Number of target individual", fill = "Relatedness degree") +
  theme_classic() +
  scale_x_continuous(breaks=1:10) +
  scale_fill_manual(values=c(cols[7], cols[6], cols[5])) +
  theme(axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size=axis_text_size),
        axis.title = element_text(size=axis_title_size),
        plot.title = element_text(size=title_size),
        plot.title.position = "plot", legend.position = 'none',
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=axis_text_size),
        legend.title = element_text(size=axis_text_size),
        legend.key.spacing.y = unit(0.15, 'cm'),
        legend.key.spacing.x = unit(0.55, 'cm'),
        legend.direction="horizontal"); g2


plot_grid(g1, g2, ncol=1, rel_heights = c(1,0.8))

ggsave('Plots/KGP_mtdna_predictions.jpg', width=8,  height=8)





