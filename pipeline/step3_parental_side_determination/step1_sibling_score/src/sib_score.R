setwd('~/Dropbox/Ressources/Git_repository/THORIN/pipeline/step3_parental_side_determination/step1_sibling_score/')
library(ggplot2)
library(cowplot)
library(parallel)
library(dplyr)
library(data.table)
library(scales)

axis_text_size=6
axis_title_size=8
title_size=10
facet_text_size=6
cols<-c("#6699CC","#CC6677","#88CCEE", "#CC6677", "#DDCC77", "#44AA99", "#AA4499", 
        "#999933", "#882255", "#661100", "#888888","#117733", "#332288")




source('src/utils.R')

###### Define main function

process_target<-function(t, mismatch=1, ibd_length=3, chromosomes=20:22, window_size=1000, threshold_diff=0){
  print(t)
  T_CO<-data.frame()
  for (CHR in chromosomes){
    tmp<-read.table(paste0('data/THORIN/KGP.chr',CHR,'.thorin.prob.ibd'), hea=T)
    t_co<-ibd2co(data=tmp, t=t, ibd_length=ibd_length, CHR=CHR)
    if (is.data.frame(t_co)){ # solve case where no IBD is found on one chromosome.
      mm<-read_rec_map(CHR, 'male')
      fm<-read_rec_map(CHR, 'female')
      t_co$bp_window_down<-bp_window(bps=t_co$bp, window_size=window_size, map=mm)[[1]]
      t_co$bp_window_up<-bp_window(bps=t_co$bp, window_size=window_size, map=mm)[[2]]
      
      # remove CO positions that are outside the recombination map (likely at the end of haplotypes)
      t_co<-t_co[t_co$bp>min(mm$pos) & t_co$bp<max(mm$pos),]
      
      if (dim(t_co)[1]>0){
        # get the Morgan value for each window position
        t_co$morgan_male_window_up<-bp2morgan(t_co$bp_window_up, map=mm)
        t_co$morgan_male_window_down<-bp2morgan(t_co$bp_window_down, map=mm)
        t_co$morgan_female_window_up<-bp2morgan(t_co$bp_window_up, map=fm)
        t_co$morgan_female_window_down<-bp2morgan(t_co$bp_window_down, map=fm)
        
        # get the diff of Morgan within the window
        t_co$morgan_diff_m<-t_co$morgan_male_window_up-t_co$morgan_male_window_down
        t_co$morgan_diff_f<-t_co$morgan_female_window_up-t_co$morgan_female_window_down
        
        # get the expected number of CO on that chromosome (not sex-specific)
        t_co$CHR<-CHR
        t_co$exp_CO<-chr2expectedCO(CHR, mm, fm)
        
        T_CO<-rbind(T_CO, t_co)
      }
    }
  }
  
  
  if (dim(T_CO)[1]>0){
    
    # look at the Morgan difference between female recomb and male recomb.
    # negative -> larger Morgan distance in males --> more likely to recombine in father.
    # positive -> larger Morgan distance in females --> more likely to recombine in mother
    T_CO$diff<-T_CO$morgan_diff_f-T_CO$morgan_diff_m 
    
    T_CO<-T_CO[abs(T_CO$diff)>=threshold_diff,]
    
    # in log10:
    # negative --> larger male log10 --> smaller male length --> more likely to recomb in mother
    T_CO$diff<-log10(T_CO$morgan_diff_f)-log10(T_CO$morgan_diff_m) 
    
    
    # Aggregate delta for the two haplotypes, across all interchrs-phased chrs 
    S<-sum(T_CO$diff[T_CO$haplotype=='H0'])-sum(T_CO$diff[T_CO$haplotype=='H1'])
    
    
    
    # check expected number of CO given the chrs that could be used
    sum_morgan<-sum(unique(T_CO$exp_CO))
    
    
    res<-data.frame(target=t, score=S, n_chrs=length(chromosomes), n_chrs_used=length(unique(T_CO$CHR)),
                    chrs=paste0(unique(T_CO$CHR), collapse=','), n_co_used=dim(T_CO)[1], mismatch=mismatch,
                    ibd_length=ibd_length, sum_morgan=sum_morgan)
    
    return(res)
  }
}



####################################################################
####################################################################
### Benchmark
####################################################################
####################################################################

rel<-read.table('../../step1_surrogate_parents/data/benchmark/Relatives.benchmark.side', hea=T); rel<-rel[,c(1,3,4)]
rel$maternal_group<-"G1"
rel$maternal_group[rel$group=='G2' & rel$side=='maternal']<-"G2"
rel$maternal_group[rel$group=='G1' & rel$side=='paternal']<-"G2"
rel<-rel[,c(1,4)]



# get the benchmark targets
tmp=paste0('data/sibs_with_phasing.group'); tmp<-read.table(tmp, hea=F)
targets<-unique(tmp$V1)
targets<-targets[targets %in% rel$target]
# on the KGP example, we have only one individual here.
# It will not be enough to derive probabilities.
# We will still go through all the steps, but make sure you have a decent number of individuals in your benchmark cohort for this method.
# For the sample, let's use more benchmark individual and assign randomly their group parental side, just to have enough individuals.
tmp=paste0('data/sibs_with_phasing.group'); tmp<-read.table(tmp, hea=F)
targets<-unique(tmp$V1)
targets<-unique(c(targets[targets %in% rel$target], targets[1:9]))


### *** Compute the score for each chromosomes
mismatch=1; ibd_length=3; threshold_diff=0
score<-data.frame()
for (N in 1:22){
  xx<-mclapply(targets, process_target,  mismatch=mismatch, ibd_length=ibd_length, chromosomes=N:N, window_size=1000, threshold_diff=threshold_diff, mc.cores=15)
  RES<-rbind(bind_rows(xx))
  if (dim(RES)[1]>0){
    RES$type<-'sib_score'
    RES$threshold_diff<-threshold_diff
    score<-rbind(score, RES)
  }
}
score$maternal_group<-rel$maternal_group[match(score$target, rel$target)]

# here, since we used also random benchmark target for the KGP example,
# we randomly assign some relative groups.
score$maternal_group[is.na(score$maternal_group)]<-'G1'
write.table(score, 'data/benchmark_score.unique_chr.txt', quote=F, col.names=T, row.names=F, sep='\t')



### *** Compute the score for different nbr of chromosomes
mismatch=1; ibd_length=3; threshold_diff=0
score<-data.frame()
for (N in 1:22){
  xx<-mclapply(targets, process_target,  mismatch=mismatch, ibd_length=ibd_length, chromosomes=1:N, window_size=1000, threshold_diff=threshold_diff, mc.cores=15)
  RES<-rbind(bind_rows(xx))
  RES$type<-'sib_score'
  RES$threshold_diff<-threshold_diff
  score<-rbind(score, RES)
}

score$maternal_group<-rel$maternal_group[match(score$target, rel$target)]
# here, since we used also random benchmark target for the KGP example,
# we randomly assign some relative groups.
score$maternal_group[is.na(score$maternal_group)]<-'G1'
score<-score[score$n_chrs_used!=0,]
write.table(score, 'data/benchmark_score.n_chrs.txt', quote=F, col.names=T, row.names=F, sep='\t')








####################################################################
####################################################################
### Call
####################################################################
####################################################################

# get targets
tmp=paste0('data/sibs_with_phasing.group'); tmp<-read.table(tmp, hea=F)
targets<-unique(tmp$V1)

# remove individual use for the benchmark
dt<-read.table('data/benchmark_score.n_chrs.txt', hea=T)
targets<-targets[!(targets %in% dt$target)]

# compute score
xx<-mclapply(targets, process_target,  mismatch=1, ibd_length=3, chromosomes=1:22, window_size=1000, threshold_diff=0, mc.cores=4)
RES<-rbind(bind_rows(xx))
write.table(RES, 'data/score_interchrs_full.txt', quote=F, col.names=T, row.names=F, sep='\t')





####################################################################
####################################################################
### Derive a probabilistic assignment
####################################################################
####################################################################

###### 
###### Benchmark
###### 



dt<-read.table('data/benchmark_score.n_chrs.txt', hea=T)

derive_prob<-function(value, cm, train, window=2){
  data<-train[train$sum_morgan>=(cm-window) & train$sum_morgan<=(cm+window),]
  tp=dim(data[data$score>=abs(value) & data$maternal_group=='G1',])[1] + dim(data[data$score<= -abs(value) & data$maternal_group=='G2',])[1]
  fp=dim(data[data$score>=abs(value) & data$maternal_group=='G2',])[1] + dim(data[data$score<= -abs(value) & data$maternal_group=='G1',])[1]
  return(tp/(tp+fp))
}

get_prob<-function(subset, data=dt){
  subset<-as.data.frame((subset))
  train<-dt[dt$target!=subset$target, ]
  
  return(derive_prob(subset$score, subset$sum_morgan, train, window=3))
}

# Assign group and probability.
dt$assignment<-NA
dt$assignment[dt$score<0]<-'G2'
dt$assignment[dt$score>0]<-'G1'
dt$prob<-NA
for (i in 1:dim(dt)[1]){
  dt$prob[i]<-get_prob(dt[i,])
}
dt$prob[abs(dt$score)>=10]<-1 # for cases when the training set does not contain large value, otherwise assign NA probability.
write.table(dt, 'data/sib_score.benchmark.txt', quote=F, col.names=T, row.names = F, sep='\t')



# Accuracy
dt<-dt[!is.na(dt$prob),]
step=0.001
tmp<-data.frame()
for (i in seq(min(dt$prob), max(dt$prob), step)){
  tp=dim(dt[dt$prob>=i & dt$maternal_group==dt$assignment,])[1]
  fp=dim(dt[dt$prob>=i & dt$maternal_group!=dt$assignment,])[1]
  tmp<-rbind(tmp, data.frame(cut=i, tp=tp, fp=fp, accuracy=tp/(tp+fp)))  
}







###### 
###### Call
###### 
dt<-read.table('data/benchmark_score.n_chrs.txt', hea=T)

score<-read.table('data/score_interchrs_full.txt', hea=T);
score$assignment<-NA; score$assignment[score$score<0]<-'G2'; score$assignment[score$score>0]<-'G1'
score$prob<-NA
for (i in 1:dim(score)[1]){
  score$prob[i]<-get_prob(score[i,])
}

score$prob[is.na(score$prob) & abs(score$score)>10]<-1





score$outliers<-0
step=.5
res<-data.frame()
for (i in seq(min(score$sum_morgan), max(score$sum_morgan)+step, step)){
  tmp<-score[score$sum_morgan>=i & score$sum_morgan<(i+step),]
  ol<-find_outliers(tmp$n_co_used, 20)
  tmp$outliers[tmp$n_co_used %in% ol]<-1
  res<-rbind(res, tmp)
}

write.table(res, 'data/Sib_score.processed.txt', quote=F, col.names=T, row.names=F, sep='\t')















####################################################################
####################################################################
### Plots
####################################################################
####################################################################

system('mkdir -p Plots/')


### Intra- vs inter-chromosomal phasing
score<-read.table('data/benchmark_score.unique_chr.txt', hea=T)
score$morgan_label<-paste0('chr',score$chrs,':',round(score$sum_morgan,2),' M')
tmp<-data.frame(chr=score$chrs, label=score$morgan_label); tmp<-tmp[!duplicated(tmp),]; tmp<-tmp[order(tmp$chr),]
score$morgan_label<-factor(score$morgan_label, levels=tmp$label)
g1<-ggplot(score, aes(x=factor(morgan_label), y=n_co_used)) + geom_boxplot(fill='grey',col='black') +
  theme_classic() +
  labs(title="a. Crossover inference \n    - intra-chromosomal phased data",
       y='Number of crossovers', x='Length of genomic region used in Morgan') +
  theme(axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size=axis_text_size),
        axis.title = element_text(size=axis_title_size),
        plot.title = element_text(size=title_size),
        plot.title.position = "plot", legend.position = 'none',
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.25, 'cm'),
        legend.text = element_text(size=axis_text_size),
        legend.title = element_text(size=axis_text_size),
        legend.key.spacing.y = unit(0.15, 'cm'),axis.text.x = element_text(angle = 45, hjust = 1)); g1
score<-read.table('data/benchmark_score.n_chrs.txt', hea=T)
g2<-ggplot(score, aes(x=sum_morgan, y=n_co_used)) + geom_point(fill='grey',col='black', shape=21) +
  theme_classic() +
  labs(title='b. Crossover inference \n    - inter-chromosomal phased data',
       y='Number of crossovers', x='Length of genomic region used in Morgan') +
  theme(axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size=axis_text_size),
        axis.title = element_text(size=axis_title_size),
        plot.title = element_text(size=title_size),
        plot.title.position = "plot", legend.position = 'none',
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.25, 'cm'),
        legend.text = element_text(size=axis_text_size),
        legend.title = element_text(size=axis_text_size),
        legend.key.spacing.y = unit(0.15, 'cm')); g2




### intra-chrs phasing
score<-read.table('data/benchmark_score.unique_chr.txt', hea=T)
score$chrs<-factor(score$chrs, levels=1:22)
g3<-ggplot(score, aes(x=factor(chrs), y=score)) +
  geom_jitter(aes(col=maternal_group), alpha=.5, width=.1) + theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed", col='grey10') +
  geom_hline(yintercept = 2, linetype="dashed", col='grey50') +
  geom_hline(yintercept = -2, linetype="dashed", col='grey50') +
  theme(legend.position = 'bottom') +
  scale_color_manual(values=c(cols[2], cols[3]), labels=c('maternal','paternal')) +
  labs(title='c. Sibling score \n    - intra-chromosomal phased data', x='Chromosomes', y='Score',col='Haplotype 0\nparental origin:') +
  theme(axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.text = element_text(size=axis_text_size),
        axis.title = element_text(size=axis_title_size),
        plot.title = element_text(size=title_size),
        plot.title.position = "plot", legend.position = c(0.7, 0.1),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=axis_text_size,margin = margin(l = -0.01, unit = "cm")),
        legend.title = element_text(size=axis_text_size),
        legend.key.spacing.y = unit(0.15, 'cm'), legend.direction = 'horizontal') ; g3


### inter-chrs phasing
score<-read.table('data/benchmark_score.n_chrs.txt', hea=T)
score<-score[score$n_chrs_used>1,]
g4<-ggplot(score, aes(x=factor(n_chrs_used), y=score)) +
  geom_jitter(aes(col=maternal_group), alpha=.5, width=.1) + theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed", col='grey10') +
  geom_hline(yintercept = 2, linetype="dashed", col='grey50') +
  geom_hline(yintercept = -2, linetype="dashed", col='grey50') +
  theme(legend.position = 'none') +
  scale_color_manual(values=c(cols[2], cols[3]), labels=c('H1','H2')) +
  labs(title='d. Sibling score \n    - inter-chromosomal phased data', x='Number of chromosomes used', y='Score') +
  theme(axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size=axis_text_size),
        axis.title = element_text(size=axis_title_size),
        plot.title = element_text(size=title_size),
        plot.title.position = "plot", legend.position = 'none',
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=axis_text_size,margin = margin(l = -0.01, unit = "cm")),
        legend.title = element_text(size=axis_text_size),
        legend.key.spacing.y = unit(0.15, 'cm')) +
  scale_x_discrete(breaks=2:22, labels=paste0('N=',2:22), guide = guide_axis(angle = 45)); g4


legend <- get_legend(
  ggplot(score, aes(x=factor(chrs), y=score)) +
    geom_jitter(aes(col=maternal_group), alpha=.5, width=.1) + theme_classic() +
    geom_hline(yintercept = 0, linetype="dashed", col='grey10') +
    geom_hline(yintercept = 2, linetype="dashed", col='grey50') +
    geom_hline(yintercept = -2, linetype="dashed", col='grey50') +
    theme(legend.position = 'bottom') +
    scale_color_manual(values=c(cols[2], cols[3]), labels=c('maternal','paternal')) +
    labs(title='c. Sibling score - intra-chromosomal phased data', x='Chromosomes', y='Score',col='Haplotype 0\nparental origin:') +
    theme(axis.title.y = element_text(margin = margin(r = 10)),
          axis.text = element_text(size=axis_text_size),
          axis.title = element_text(size=axis_title_size),
          plot.title = element_text(size=title_size),
          plot.title.position = "plot", legend.position = 'left',
          legend.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          legend.key.size = unit(0.4, 'cm'),
          legend.text = element_text(size=axis_text_size,margin = margin(l = -0.01, unit = "cm")),
          legend.title = element_text(size=axis_text_size),
          legend.key.spacing.y = unit(0.15, 'cm')) 
)





res<-data.frame()
for (cut in seq(0,10,0.5)){
  #intra
  score<-read.table('data/benchmark_score.unique_chr.txt', hea=T)
  call<-dim(score[score$score< (-cut) | score$score>cut,])[1]/dim(score)[1]
  tmp<-score[score$score> (-cut) | score$score<cut,]
  fn=dim(score[score$maternal_group=='G1' & score$score<0,])[1] + dim(score[score$maternal_group=='G2' & score$score>0,])[1]
  score<-score[score$score< (-cut) | score$score>cut,]
  tp=dim(score[score$maternal_group=='G1' & score$score>0,])[1] + dim(score[score$maternal_group=='G2' & score$score<0,])[1]
  fp=dim(score[score$maternal_group=='G1' & score$score<0,])[1] + dim(score[score$maternal_group=='G2' & score$score>0,])[1]
  errors_rate=fp/(tp+fp)
  recall<-tp/(tp+fn)
  res<-rbind(res, data.frame(cut=cut, errors_rate=errors_rate, call_rate=call, recall=recall, type="intra"))
  
  #inter
  score<-read.table('data/benchmark_score.n_chrs.txt', hea=T)
  score<-score[score$n_chrs_used>1,]
  call=dim(score[score$score< (-cut) | score$score>cut,])[1]/dim(score)[1]
  tmp<-score[score$score> (-cut) | score$score<cut,]
  fn=dim(score[score$maternal_group=='G1' & score$score<0,])[1] + dim(score[score$maternal_group=='G2' & score$score>0,])[1]
  score<-score[score$score< (-cut) | score$score>cut,]
  tp=dim(score[score$maternal_group=='G1' & score$score>0,])[1] + dim(score[score$maternal_group=='G2' & score$score<0,])[1]
  fp=dim(score[score$maternal_group=='G1' & score$score<0,])[1] + dim(score[score$maternal_group=='G2' & score$score>0,])[1]
  errors_rate=fp/(tp+fp)
  recall<-tp/(tp+fn)
  res<-rbind(res, data.frame(cut=cut, errors_rate=errors_rate, call_rate=call, recall=recall, type="inter"))
}

g5<-ggplot() +
  geom_point(data=res, aes(x=cut, y=errors_rate*100, col=factor(type))) +
  geom_line(data=res, aes(x=cut, y=errors_rate*100, col=factor(type))) +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(title='e. Sibling score error rate', x='Threshold on score (abs value)', y='Error rate',  col='Phased data') +
  theme(axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size=axis_text_size),
        axis.title = element_text(size=axis_title_size),
        plot.title = element_text(size=title_size),
        plot.title.position = "plot", legend.position = c(0.75, 0.6),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=axis_text_size,margin = margin(l = -0.01, unit = "cm")),
        legend.title = element_text(size=axis_text_size),
        legend.key.spacing.y = unit(0.15, 'cm'))  +
  scale_color_manual(values=cols[5:6], labels=c('inter-chromosomal phasing', 'intra-chromosomal phasing')) +
  scale_y_continuous(breaks=c(0,1,2.5,5,10,15,20, 25), labels=paste0(c(0,1,2.5,5,10,15,20, 25),'%')) +
  scale_x_continuous(breaks=seq(0,10), labels=seq(0,10))

g6<-ggplot() +
  geom_point(data=res, aes(x=cut, y=call_rate*100, col=factor(type))) +
  geom_line(data=res, aes(x=cut, y=call_rate*100, col=factor(type))) +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(title='f. Sibling score call rate', x='Threshold on score (abs value)', y='Call rate',  col='') +
  theme(axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size=axis_text_size),
        axis.title = element_text(size=axis_title_size),
        plot.title = element_text(size=title_size),
        plot.title.position = "plot", legend.position = 'none',
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=axis_text_size,margin = margin(l = -0.01, unit = "cm")),
        legend.title = element_text(size=axis_text_size),
        legend.key.spacing.y = unit(0.15, 'cm')) + guides(col = guide_legend(title.position = "left")) +
  scale_color_manual(values=cols[5:6], labels=c('inter-chromosomal phasing', 'intra-chromosomal phasing')) +
  scale_y_continuous(breaks=seq(0,100,25), labels=paste0(seq(0,100,25),'%')) +
  scale_x_continuous(breaks=seq(0,10), labels=seq(0,10))

legend2 <- get_legend(
  ggplot() +
    geom_point(data=res, aes(x=cut, y=call_rate*100, col=factor(type))) +
    geom_line(data=res, aes(x=cut, y=call_rate*100, col=factor(type))) +
    theme_classic() +
    theme(legend.position = 'none') +
    labs(title='f. Sibling score call rate', x='Threshold on score (abs value)', y='Call rate',  col='Phased data') +
    theme(axis.title.y = element_text(margin = margin(r = 10)),
          axis.text = element_text(size=axis_text_size),
          axis.title = element_text(size=axis_title_size),
          plot.title = element_text(size=title_size),
          plot.title.position = "plot", legend.position = 'left',
          legend.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          legend.key.size = unit(0.4, 'cm'),
          legend.text = element_text(size=axis_text_size),
          legend.title = element_text(size=axis_text_size),
          legend.key.spacing.y = unit(0.15, 'cm')) +
    scale_color_manual(values=cols[5:6], labels=c('inter-chromosomal', 'intra-chromosomal')) +
    scale_y_continuous(breaks=seq(0,100,25), labels=paste0(seq(0,100,25),'%')) +
    scale_x_continuous(breaks=seq(0,10), labels=seq(0,10))
)




a1<-plot_grid(g1, NULL,plot_grid(g2,NULL, ncol=1, rel_heights=c(1,0.13)), rel_widths = c(1,0.1,1), ncol=3); a1
a2<-plot_grid(g3, NULL, g4, ncol=3, rel_widths = c(1,0.15,1)); a2
a3<-plot_grid(g5, NULL, g6, ncol=3, rel_widths = c(1,0.15,1)); a3
plot_grid(a1, NULL, a2, NULL, a3, ncol=1, rel_heights = c(1,0.1,1,0.1,1))
ggsave('Plots/sib_score.benchmark.jpg', width=8,  height=10)



















### 2. Probabilities
dt<-read.table('data/sib_score.benchmark.txt', hea=T)
dt$col<-0; dt$col[dt$maternal_group==dt$assignment]<-1
dt$best<-0
for (t in unique(dt$target)){
  tmp<-dt[dt$target==t,]
  dt$col[dt$target==t & dt$sum_morgan==max(tmp$sum_morgan)]<-2
  dt$best[dt$target==t & dt$sum_morgan==max(tmp$sum_morgan)]<-1
}

g1<-ggplot(dt, aes(x=score, y=prob, col=factor(col))) + geom_point()+
  theme_classic() +
  theme(legend.position = 'bottom') +
  scale_color_manual(values=c(cols[4], cols[6], cols[5]), labels=c('incorrect','correct','correct (all available chrs)')) +
  labs(title='a. Probabilities vs. sib-score - validation cohort', x='Score', y='Probability', col='Assignment :  ') +
  theme(axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size=axis_text_size),
        axis.title = element_text(size=axis_title_size),
        plot.title = element_text(size=title_size),
        plot.title.position = "plot", legend.position = c(0.75, 0.5),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=axis_text_size, margin = margin(l = -0.01, unit = "cm")),
        legend.title = element_text(size=axis_text_size),
        legend.key.spacing.y = unit(0.15, 'cm')); g1

g2<-ggplot(dt, aes(x=sum_morgan, y=prob, col=factor(col))) + geom_point() +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_color_manual(values=c(cols[4], cols[6], cols[5]), labels=c('incorrect','correct','correct - all available chromosomes')) +
  labs(title='b. Probabilities vs. Morgan length - validation cohort', x='Length of genomic region used in Morgan', y='Probability',  col='Assignment :  ') +
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
        legend.key.spacing.x = unit(0, 'cm')) + guides(col = guide_legend(title.position = "left")); g2


legend1 <- get_legend(
  ggplot(dt, aes(x=score, y=prob, col=factor(col))) + geom_point()+
    theme_classic() +
    theme(legend.position = 'bottom') +
    scale_color_manual(values=c(cols[4], cols[6], cols[5]), labels=c('incorrect','correct','correct\n(all available chrs)')) +
    labs(title='b. Probabilities vs. sib-score - validation cohort', x='Score', y='Probability', col='Assignment :  ') +
    theme(axis.title.y = element_text(margin = margin(r = 10)),
          axis.text = element_text(size=axis_text_size),
          axis.title = element_text(size=axis_title_size),
          plot.title = element_text(size=title_size),
          plot.title.position = "plot", legend.position = 'left',
          legend.background = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank(),
          legend.key.size = unit(0.4, 'cm'),
          legend.text = element_text(size=axis_text_size, margin = margin(l = -0.01, unit = "cm")),
          legend.title = element_text(size=axis_text_size),
          legend.key.spacing.y = unit(0.15, 'cm'))
)


a1<-plot_grid(g1, legend1, g2, ncol=3, rel_widths = c(1,0.25,1)); a1



# Call- Crossovers
dt<-read.table('data/Sib_score.processed.txt', hea=T)
g3<-ggplot(dt[dt$outliers==0,], aes(x=sum_morgan, y=n_co_used)) + geom_point(fill='grey',col='black', shape=21, alpha=.2) +
  theme_classic() +
  labs(title=paste0("c. Crossovers inference for N=26,635 individuals"),
       y='Number of crossovers inferred', x='Length of available genomic region in Morgan') +
  theme(axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size=axis_text_size),
        axis.title = element_text(size=axis_title_size),
        plot.title = element_text(size=title_size),
        plot.title.position = "plot", legend.position = 'none',
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.25, 'cm'),
        legend.text = element_text(size=axis_text_size),
        legend.title = element_text(size=axis_text_size),
        legend.key.spacing.y = unit(0.15, 'cm')); g3




# Call -
score<-read.table('data/Sib_score.processed.txt', hea=T); score$maternal_group<-NA
res<-score[score$outliers==0,]
sub<-res[!is.na(res$maternal_group) & res$outliers==0,]
score$assignment[score$score>=0]<-'Paternal'; score$assignment[score$score<0]<-'Maternal'
score$assignment<-factor(score$assignment, levels<-c('Paternal', 'Maternal'))

g5<-ggplot() +
  geom_jitter(data=res, aes(x=score, y=n_co_used, fill=assignment), alpha=.5, width=.1, shape=21, col='grey10') +
  theme_classic() +
  theme(legend.position = c(0.8,0.1)) +
  geom_vline(xintercept = 0, linetype="dashed", col='grey10') +
  
  labs(title=paste0('d. Inter-chromosomal sibling score for N=',comma_format()(dim(score)[1]),' individuals\n'),
       y='Number of crossovers used', x='Score', fill='Inferred origin for haplotype 0:') +
  theme(axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size=axis_text_size),
        axis.title = element_text(size=axis_title_size),
        plot.title = element_text(size=title_size),
        plot.title.position = "plot", legend.position = 'none',
        strip.background = element_blank(),
        legend.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=axis_text_size),
        legend.title = element_text(size=axis_text_size),
        legend.key.spacing.y = unit(0.05, 'cm')) +
  scale_fill_manual(values=c(cols[2],cols[3], 'black'), labels=c(paste0('Maternal (N=',comma_format()(dim(res[res$assignment=='G1',])[1]),')'),
                                                                 paste0('Paternal (N=',comma_format()(dim(res[res$assignment=='G2',])[1]),')'))) +
  guides(fill = guide_legend(title.position = "left")); g5

g52<-ggplot() +
  geom_jitter(data=res, aes(y=score, x=sum_morgan, fill=assignment), alpha=.5, width=.1, shape=21, col='grey10') +
  theme_classic() +
  geom_hline(yintercept = 0, linetype="dashed", col='grey10') +
  labs(title=paste0('c. Inter-chromosomal sibling score for N=',comma_format()(dim(score)[1]),' individuals\n'),
       x='Length of available genomic region in Morgan', y='Score', fill='Inferred origin\nfor haplotype 0:') +
  theme(axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size=axis_text_size),
        axis.title = element_text(size=axis_title_size),
        plot.title = element_text(size=title_size),
        plot.title.position = "plot", legend.position = c(0.15,0.8),
        strip.background = element_blank(),
        legend.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        legend.text = element_text(size=axis_text_size),
        legend.title = element_text(size=axis_text_size),
        legend.key.spacing.y = unit(0.05, 'cm')) +
  scale_fill_manual(values=c(cols[2],cols[3], 'black'), labels=c(paste0('Maternal (N=',comma_format()(dim(score[score$assignment=='Maternal',])[1]),')'),
                                                                 paste0('Paternal (N=',comma_format()(dim(score[score$assignment=='Paternal',])[1]),')'))) +
  guides(fill = guide_legend(title.position = "left")); g52


dt<-read.table('data/Sib_score.processed.txt', hea=T)
dt<-dt[dt$outliers==0,]
step=0.0025
d_data<-data.frame(bin=c(0.5,0.5, min(dt$prob)-0.005, min(dt$prob)-0.005), n=c(0,0), type=c('Paternal', 'Maternal'))
for (bin in seq(min(dt$prob), max(dt$prob)+step, step)){
  nm<-dim(dt[dt$prob>=bin & dt$prob<(bin+step) & dt$assignment=='G1',])[1]
  np<-dim(dt[dt$prob>=bin & dt$prob<(bin+step) & dt$assignment=='G2',])[1]
  d_data<-rbind(d_data, data.frame(bin=bin, n=nm, type='Maternal'))
  d_data<-rbind(d_data, data.frame(bin=bin, n=np, type='Paternal'))
}

g6<-ggplot(data=d_data, aes(x=bin, y=n)) + geom_area(aes(fill=type), col='black', alpha=.5) + facet_grid(type~.) +
  theme_classic() +
  labs(title=paste0("d. Haplotypes' parent-of-origin probabilities for N=",comma_format()(dim(score)[1])," individuals"),
       x='Probability', y='Number of individuals', fill='Inferred origin\nfor haplotype 0:') +
  theme(axis.title.y = element_text(margin = margin(r = 10)),
        axis.text = element_text(size=axis_text_size),
        axis.title = element_text(size=axis_title_size),
        plot.title = element_text(size=title_size),
        plot.title.position = "plot", legend.position = c(0.2,0.8),
        legend.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.key.size = unit(0.25, 'cm'),
        legend.text = element_text(size=axis_text_size),
        legend.title = element_text(size=axis_text_size),
        legend.key.spacing.y = unit(0.15, 'cm')) +
  scale_fill_manual(values=c(cols[2], cols[3]), labels=c(paste0('Maternal (N=',comma_format()(dim(score[score$assignment=='Maternal',])[1]),')'),
                                                         paste0('Paternal (N=',comma_format()(dim(score[score$assignment=='Paternal',])[1]),')'))) +
  scale_y_continuous(breaks=c(0, 2500, 5000, 7500, 10000), labels=c(0, 2500, 5000, 7500, 10000)) +
  guides(fill = guide_legend(title.position = "left")); g6




a1<-plot_grid(g1, g2, ncol=2)
a3<-plot_grid(g52, NULL, g6, ncol=1, rel_heights = c(1,0.1, 1.2)); a3
plot_grid(a1, NULL, a3, ncol=1, rel_heights = c(1,0.1,2.3))
ggsave('Plots/sib_score.call.jpg', width=8,  height=10)



















