setwd('~/Dropbox/Ressources/Git_repository/THORIN/pipeline/Parent_of_Origin/step5_interchromosomal_phasing/')
library(dplyr)
library(ggplot2)
library(cowplot)



# Read HMM output file.
hmm<-as.data.frame(data.table::fread('data/THORIN/benchmark/KGP.chr20.benchmark.thorin.prob', hea=T))
n<-unlist(lapply(colnames(hmm), FUN=function(x){unlist(strsplit(x,'_'))[1]})) ; colnames(hmm)<-n #rename cols

# Read group file - to know the number of groups to plot
groups<-as.data.frame(data.table::fread('../step1_surrogate_parents/data/benchmark/Relatives.benchmark.group', hea=F))

# Read IBD tracks file
ibd<-as.data.frame(data.table::fread('data/THORIN/benchmark/KGP.chr20.benchmark.thorin.prob.ibd', hea=T))




plot_target<-function(t){
  
  # get the number of groups for this target.
  ng<-2;
  if (groups$V3[groups$V1==t]==''){ng<-1}
  
  # subet hmm data
  sub<-cbind(hmm[, 1:4], hmm[colnames(hmm)==t])
  
  # rename cols according to groups
  if (ng==1){
    colnames(sub)[5:8]<-c('H0G1','H0U','H1G1','H1U')
    h0<-data.frame(pos=rep(sub$POS,2), prob=c(sub$H0G1, sub$H0U), type=c(rep('G1', dim(sub)[1]), rep('U', dim(sub)[1])) )
    h1<-data.frame(pos=rep(sub$POS,2), prob=c(sub$H1G1, sub$H1U), type=c(rep('G1', dim(sub)[1]), rep('U', dim(sub)[1])) )
    cols<-c('red','grey')
  }
  
  if (ng==2){
    colnames(sub)[5:10]<-c('H0G1', 'H0G2', 'H0U','H1G1', 'H1G2', 'H1U')
    h0<-data.frame(pos=rep(sub$POS,3), prob=c(sub$H0G1, sub$H0G2, sub$H0U), type=c(rep('G1', dim(sub)[1]), rep('G2', dim(sub)[1]), rep('U', dim(sub)[1])) )
    h1<-data.frame(pos=rep(sub$POS,3), prob=c(sub$H1G1, sub$H1G2, sub$H1U), type=c(rep('G1', dim(sub)[1]), rep('G2', dim(sub)[1]), rep('U', dim(sub)[1])) )
    cols<-c('red', 'blue','grey')
  }
  
  # subset ibd tracks
  ibd_t<-ibd[ibd$target==t,]
  
  # plot
  g0<-ggplot() +
    geom_point(data=h0, aes(x=pos, y=prob, col=type)) + 
    geom_line(data=h0, aes(x=pos, y=prob, col=type)) +
    scale_color_manual(values=cols) +
    theme_classic() +
    geom_vline(xintercept = ibd_t$start) +
    geom_text(data=ibd_t, aes(x=start+(end-start)/2, y=1, label=Prob))
  
  g1<-ggplot() +
    geom_point(data=h1, aes(x=pos, y=prob, col=type)) + 
    geom_line(data=h1, aes(x=pos, y=prob, col=type)) +
    scale_color_manual(values=cols) +
    theme_classic() +
    geom_vline(xintercept = ibd_t$start) +
    geom_text(data=ibd_t, aes(x=start+(end-start)/2, y=1, label=Prob))
  
  a<-plot_grid(g0, g1, ncol=1)
  return(a)
}

targets<-unique(n[5:length(n)])


plot_target('HG00420')
plot_target('NA19685')






