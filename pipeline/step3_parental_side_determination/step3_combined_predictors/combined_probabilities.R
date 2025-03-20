'%ni%'=Negate('%in%')
library(dplyr)
library(parallel)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(ggdist)
library(ggh4x)
library(scales)
library(tidyverse)
library(ggdist)
library(ggh4x)
library(scales)
library(ggforce)
axis_text_size=8
axis_title_size=10
title_size=12
facet_text_size=8
col_selection<-c("#6699CC","#CC6677","#88CCEE", "#CC6677", "#DDCC77", "#44AA99", "#AA4499", 
        "#999933", "#882255", "#661100", "#888888","#117755", "#332288")


system('mkdir -p data/')
system('mkdir -p Plots/')



### *** Read data

# chr X
dtx<-read.table('../step0_chrX_ibd/data/chrx_accuracy_derived_prediction.txt', hea=T); dtx<-dtx[!is.na(dtx$target),]

# sibs score
dts<-read.table('../step1_sibling_score/data/Sib_score.processed.txt', hea=T);

# mtDNA
dtM<-read.table('../step2_mtDNA/data/mtdna_accuracy_derived_prediction.uniq_per_degree.txt', hea=T) ; dtM<-dtM[!is.na(dtM$target),]



# Build summary data
dt<-data.frame(target=unique(c(dtx$target, dtM$target, dts$target)))

dt$chrx<-dtx$chrx[match(dt$target, dtx$target)]
dt$chrx_group<-dtx$maternal_group[match(dt$target, dtx$target)]
dt$chrx_accuracy<-dtx$accuracy[match(dt$target, dtx$target)]
dt$chrx_accuracy_type<-dtx$accuracy_type[match(dt$target, dtx$target)]
dt$chrx_prob<-dtx$maternal_probability[match(dt$target, dtx$target)]


dtm<-dtM[dtM$degree=='2nd',]
dt$mtdna2<-dtm$mtdna[match(dt$target, dtm$target)]
dt$mtdna2_group<-dtm$maternal_group[match(dt$target, dtm$target)]
dt$mtdna2_accuracy<-dtm$accuracy[match(dt$target, dtm$target)]
dt$mtdna2_accuracy_type<-dtm$accuracy_type[match(dt$target, dtm$target)]
dt$mtdna2_prob<-dtm$maternal_probability[match(dt$target, dtm$target)]



dtm<-dtM[dtM$degree=='3rd',]
dt$mtdna3<-dtm$mtdna[match(dt$target, dtm$target)]
dt$mtdna3_group<-dtm$maternal_group[match(dt$target, dtm$target)]
dt$mtdna3_accuracy<-dtm$accuracy[match(dt$target, dtm$target)]
dt$mtdna3_accuracy_type<-dtm$accuracy_type[match(dt$target, dtm$target)]
dt$mtdna3_prob<-dtm$maternal_probability[match(dt$target, dtm$target)]


dtm<-dtM[dtM$degree=='4th',]
dt$mtdna4<-dtm$mtdna[match(dt$target, dtm$target)]
dt$mtdna4_group<-dtm$maternal_group[match(dt$target, dtm$target)]
dt$mtdna4_accuracy<-dtm$accuracy[match(dt$target, dtm$target)]
dt$mtdna4_accuracy_type<-dtm$accuracy_type[match(dt$target, dtm$target)]
dt$mtdna4_prob<-dtm$maternal_probability[match(dt$target, dtm$target)]



dt$sib_score<-dts$score[match(dt$target, dts$target)]
dt$sib_group<-dts$assignment[match(dt$target, dts$target)]
dt$sib_prob<-dts$prob[match(dt$target, dts$target)]
dt$sib_accuracy<-dts$prob[match(dt$target, dts$target)]


dt<-dt[!is.na(dt$target),]
dt[dt=='G']<-'G1' # when only one group has been used, THORIN reports "G" and not "G1"


# Keep the records with the highest accuracy per target
# when several records have the same accuracy, keep all. Then, check if the group matches.
dt <- as.data.frame( dt %>%
                       rowwise() %>%
                       mutate(
                         max_value = max(c(chrx_accuracy, mtdna2_accuracy, mtdna3_accuracy, mtdna4_accuracy, sib_accuracy), na.rm=T),
                         max_name = paste(c("chrx", "mtdna_2nd", "mtdna_3rd", "mtdna_4th", "sib_score")[c(chrx_accuracy, mtdna2_accuracy, mtdna3_accuracy, mtdna4_accuracy, sib_accuracy) == max_value][!is.na(c("chrx", "mtdna_2nd", "mtdna_3rd", "mtdna_4th", "sib_score")[c(chrx_accuracy, mtdna2_accuracy, mtdna3_accuracy, mtdna4_accuracy, sib_accuracy) == max_value])], collapse = ","),
                         maternal_prob = paste(unique(c(chrx_prob, mtdna2_prob, mtdna3_prob, mtdna4_prob, sib_prob)[c(chrx_accuracy, mtdna2_accuracy, mtdna3_accuracy, mtdna4_accuracy, sib_accuracy) == max_value][!is.na(c(chrx_prob, mtdna2_prob, mtdna3_prob, mtdna4_prob, sib_prob)[c(chrx_accuracy, mtdna2_accuracy, mtdna3_accuracy, mtdna4_accuracy, sib_accuracy) == max_value])]), collapse = ","),
                         maternal_prob_n = length(c(chrx_prob, mtdna2_prob, mtdna3_prob, mtdna4_prob, sib_prob)[c(chrx_accuracy, mtdna2_accuracy, mtdna3_accuracy, mtdna4_accuracy, sib_accuracy) == max_value][!is.na(c(chrx_prob, mtdna2_prob, mtdna3_prob, mtdna4_prob, sib_prob)[c(chrx_accuracy, mtdna2_accuracy, mtdna3_accuracy, mtdna4_accuracy, sib_accuracy) == max_value])]),
                         maternal_prob_group = paste(c(chrx_group, mtdna2_group, mtdna3_group, mtdna4_group, sib_group)[c(chrx_accuracy, mtdna2_accuracy, mtdna3_accuracy, mtdna4_accuracy, sib_accuracy) == max_value][!is.na(c(chrx_group, mtdna2_group, mtdna3_group, mtdna4_group, sib_group)[c(chrx_accuracy, mtdna2_accuracy, mtdna3_accuracy, mtdna4_accuracy, sib_accuracy) == max_value])], collapse = ","),
                         maternal_prob_group_n = length(unique((c(chrx_group, mtdna2_group, mtdna3_group, mtdna4_group, sib_group)[c(chrx_accuracy, mtdna2_accuracy, mtdna3_accuracy, mtdna4_accuracy, sib_accuracy) == max_value][!is.na(c(chrx_group, mtdna2_group, mtdna3_group, mtdna4_group, sib_group)[c(chrx_accuracy, mtdna2_accuracy, mtdna3_accuracy, mtdna4_accuracy, sib_accuracy) == max_value])]))))
                     %>%
                       ungroup() )



# when multiple groups, check concordance rate.
t<-data.frame(table(dt$maternal_prob_group))
t$N<-unlist(lapply(as.character(t$Var1), FUN=function(x){length(unlist(strsplit(x,',')))}))
t$N_unique<-unlist(lapply(as.character(t$Var1), FUN=function(x){length(unique(unlist(strsplit(x,','))))}))

tt<-t[t$N>1,]
a<-sum(tt$Freq[tt$N_unique==1])
b<-sum(tt$Freq[tt$N_unique==2])
concordance<-a/(a+b)


# Remove target with both group assign to the same parent 
dt0<-dt[dt$maternal_prob_group_n==2,]; dim(dt0); length(unique(dt0$target))
dt1<-dt[dt$maternal_prob_group_n==1,]; dim(dt1); length(unique(dt1$target))
prop_removed<-(dim(dt0)[1]/(dim(dt0)[1] + dim(dt1)[1]))*100; prop_removed
prop_kept<-(dim(dt1)[1]/(dim(dt0)[1] + dim(dt1)[1]))*100; prop_kept





# for dt1, keep only one occurence (eg. reduce "G1,G1,G1" to only "G1" in the data.frame )
reduce_repetitions <- function(group) {
  elements <- str_split(group, ",")[[1]]
  unique_element <- unique(elements)
  if (length(unique_element) == 1) {
    return(unique_element)
  } else {
    return(group)
  }
}

dt1 <- as.data.frame( dt1 %>%
                        mutate(most_frequent = sapply(maternal_prob_group, reduce_repetitions)))





# For target with more than one group, keep the most frequent group if it's not only once occurrence of G1 and one occurrence of G1
more_frequent_value <- function(group) {
  elements <- str_split(group, ",")[[1]]
  freq_table <- table(elements)
  if(length(elements) > 2) {
    if(freq_table["G1"] > freq_table["G2"]) {
      return("G1")
    } else if(freq_table["G2"] > freq_table["G1"]) {
      return("G2")
    } else {
      return(NA)
    }
  } else {
    return(NA)
  }
}

dt0 <- dt0 %>%
  mutate(most_frequent = sapply(maternal_prob_group, more_frequent_value))

dt0mf<-dt0[!is.na(dt0$most_frequent),] # those with only "G1,G2"
dt0mf$maternal_prob_group<-dt0mf$most_frequent



###
### Start comment: Due to the small number of individuals when using the 1KGP example data, this part does not work since all target are already correctly assigned.
###

# For target with only G1 and G2, prioritizing models: 1-chrx; 2-sibs with high score; 3- 2nd degree
dt0<-dt0[is.na(dt0$most_frequent),] # those with only "G1,G2"


# For target with both group assign to the same parent, re-run by prioritizing models: 1-chrx; 2-sibs with high score; 3- 2nd degree

# chrx
dt0x<-as.data.frame(dt0 %>% filter(grepl("chrx", max_name)))
dt0x<-dt0x[dt0x$chrx>5,]
dt0x$maternal_prob_group<-dt0x$chrx_group
dt0x$most_frequent<-dt0x$chrx_group
dt0x$maternal_prob<-dt0x$chrx_prob
if (dim(dt0x)[1]!=0){dt0x$max_name<-'chrx'}
dt0<-dt0[!(dt0$target %in% dt0x$target),]


#sibs-score
dt0s<-as.data.frame(dt0 %>% filter(grepl("sib_score", max_name)))
dt0s<-dt0s[abs(dt0s$sib_score)>10,]
dt0s$maternal_prob_group<-dt0s$sib_group
dt0s$most_frequent<-dt0s$sib_group
dt0s$maternal_prob<-dt0s$sib_prob
if (dim(dt0s)[1]!=0){dt0s$max_name<-'sib_score'}
dt0<-dt0[!(dt0$target %in% dt0s$target),]


# 2nd degree
dt02<-as.data.frame(dt0 %>% filter(grepl("mtdna_2nd", max_name)))
dt02$maternal_prob_group<-dt02$mtdna2_group
dt02$most_frequent<-dt02$mtdna2_group
dt02$maternal_prob<-dt02$mtdna2_prob
if (dim(dt02)[1]!=0){dt02$max_name<-'mtdna_2nd'}
dt0<-dt0[!(dt0$target %in% dt02$target),]


# Remain un-assigned:
table(dt0$max_name)
sum(table(dt0$max_name))

###
### End comment: Due to the small number of individuals when using the 1KGP example data, this part does not work since all target are already correctly assigned.
###




# Final calls.
dt<-rbind(dt1, dt0mf, dt0x, dt0s, dt02)
write.table(dt, 'data/joint_probabilities.txt', quote=F, col.names=T, row.names=F, sep='\t')


















### ***
### *** Plots
### ***
library("viridis")
col_selection<-c("#DDCC77","#999933","#44AA59","#416755", "#008080", "#6699CC")

col_selection<-viridis(7)





dt<-read.table('data/joint_probabilities.txt', hea=T)

trios<-read.table('../../step1_surrogate_parents/data/Trios.ped')
duos<-read.table('../../step1_surrogate_parents/data/Duos.ped')
NN<-length(unique(c(dt$target, trios$V1, duos$V1))); NN


## Probability selection
t<-as.data.frame(table(dt$max_name))
t<-t[order(t$Freq),]
t$Var1<-factor(t$Var1, levels = t$Var1)
more_than_two_elements<-grepl(".*,.*,.*", t$Var1)
t2<-t[more_than_two_elements, ]
t<-t[!more_than_two_elements, ]
t$Var1<-factor(t$Var1, levels = t$Var1)
t<-rbind(t, data.frame(Var1="parental genomes", Freq=(length(unique(trios$V1))+length(unique(duos$V1)))))
t$Freq<-as.numeric(t$Freq)
t<-t[order(t$Freq),]
t<-rbind(data.frame(Var1=">=3", Freq=sum(t2$Freq)),t)
t$Var1<-factor(t$Var1, levels=t$Var1)

g2<-ggplot(t, aes(x=Var1, y=Freq, fill=Var1))+ geom_bar(stat='identity', col='black', alpha=.9) +
  theme_classic() +
  labs(title=paste0('b. Probability selection for N=',comma_format()(NN),' individuals'),
       x='Parent-of-Origin determined from', y='Number of individuals', fill='We predict') +
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
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.2)) +
  scale_x_discrete(guide = guide_axis(angle = 45),, labels=c('>=3 predictors','mtdna 2nd + 3rd degree', 'mtdna 2nd + 4th degree', 'mtdna 3rd + 4th degree',
                                                             'mtdna 2nd degree + sib score (inter-)', 'mtdna 4th degree + sib score (inter-)', 'mtdna 3rd degree + sib score (inter-)',
                                                             'chr X + sib score (inter-)', 'chr X + mtdna 2nd degree', 'parental genomes', 'chr X + mtdna 4th degree', 'chr X + mtdna 3rd degree',
                                                             'mtdna 2nd degree', 'sib score (intra-)', 'sib score (inter-)', 'mtdna 3rd degree', 'chr X', 'mtdna 4th degree')) +
  scale_fill_manual(values=c(rep('grey', 9), col_selection[7], 'grey','grey', col_selection[1:6]  )); g2

t_order<-t

## Probability frequency binned and colored by selection
d_data<-dt
d_data$max_name[!(d_data$max_name %in% c('mtdna_4th', 'chrx', 'mtdna_3rd', 'sib_score', 'mtdna_2nd'))]<-'>1 predictor'
d_data$bin<-NA
step=0.01
for (bin in seq(0.5, 1, step)){
  d_data$bin[d_data$maternal_prob>bin & d_data$maternal_prob<=(bin+step)]<-paste0(bin,'-',bin+step)
}

d<-d_data[,c(26,32)]; d<-d[!is.na(d$bin),]
d<-rbind(d,data.frame(max_name=rep('parental genomes', (length(unique(trios$V1))+length(unique(duos$V1)))), bin=rep('0.99-1', (length(unique(trios$V1))+length(unique(duos$V1))))))

d$max_name<-factor(d$max_name, levels = c(c('mtdna_4th', 'chrx', 'mtdna_3rd', 'sib_score', 'mtdna_2nd','parental genomes', '>1 predictor')))

t_order<-t_order[t_order$Var1 %in% d$max_name,]
d$max_name<-factor(d$max_name, levels=c(rev(as.character(t_order$Var1)), '>1 predictor'))

ylab_breaks=seq(0,125000,25000)
ylab_labels=unlist(lapply(ylab_breaks, FUN=function(x){comma_format()(x)}))

g4<-ggplot(d, aes(x=factor(bin), fill=max_name)) + geom_bar() +
  theme_classic() +
  labs(title=paste0('a. Probability distribution for N=',comma_format()(NN),' individuals'), x='Probability bin', y='Number of individuals', fill='Assigned from: ') +
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
        legend.direction="horizontal") +
  scale_fill_manual(values=c(rev(col_selection[1:6]), col_selection[7], 'grey')) +
  scale_y_continuous(breaks=ylab_breaks, labels=ylab_labels) +
  scale_x_discrete(guide = guide_axis(angle = 45)); g4


plot_grid(g4, NULL, g2,  ncol=1, rel_heights = c(1,0.05,1))


ggsave('Plots/joint_probability.jpg', width=8.2, height=10)



######################################################################################################################################################
##### Bring back Benchmark set , Trios, Duos ####################################################################################
######################################################################################################################################################
call<-read.table('data/joint_probabilities.txt', hea=T)
bench<-read.table('../../step1_surrogate_parents/data/benchmark/Relatives.benchmark.side', hea=T)
trios<-read.table('../../step1_surrogate_parents/data/Trios.ped')
duos<-read.table('../../step1_surrogate_parents/data/Duos.ped')


call<-call[,c(1, 29, 27, 26)]
colnames(call)<-c('target','maternal_group','maternal_probability','predictors')

# ped phasing --> h0=father ; scaffold phasing --> h0=G1

call<-call[!(call$target %in% c(trios$V1, duos$V1)),]
call<-rbind(call, data.frame(target=trios$V1, maternal_group='G2', maternal_probability=1, predictors='trios'))
call<-rbind(call, data.frame(target=duos$V1, maternal_group='G2', maternal_probability=1, predictors='duos'))

call$probability_G1_maternal<-NA
call$probability_G1_maternal[call$maternal_group=='G1']<-call$maternal_probability[call$maternal_group=='G1']
call$probability_G1_maternal[call$maternal_group=='G2']<- 1-call$maternal_probability[call$maternal_group=='G2']
call<-call[complete.cases(call),]

ggplot(call, aes(x=maternal_probability, y=probability_G1_maternal, col=maternal_group)) + geom_point() # visual check that reversing the probabilities worked


write.table(call, 'PofO_probability.txt', quote=F, col.names=T, row.names=F, sep='\t')


