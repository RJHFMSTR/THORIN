library(dplyr)
library(ggplot2)
library(cowplot)



# Read HMM output file.
hmm<-as.data.frame(data.table::fread('data/THORIN/KGP.chr20.thorin.prob', hea=T))
n<-unlist(lapply(colnames(hmm), FUN=function(x){unlist(strsplit(x,'_'))[1]})) ; colnames(hmm)<-n #rename cols

# Read group file - to know the number of groups to plot
groups<-as.data.frame(data.table::fread('../step1_surrogate_parents/data/Relatives.group', hea=F))

# Read IBD tracks file
ibd<-as.data.frame(data.table::fread('data/THORIN/KGP.chr20.thorin.prob.ibd', hea=T))




plot_target<-function(t){
  
  # get the number of groups for this target.
  ng<-2;
  if (groups$V3[groups$V1==t]==''){ng<-1}
  
  # subet hmm data
  sub<-cbind(hmm[, 1:4], hmm[colnames(hmm)==t])
  
  # rename cols according to groups
  if (ng==1){
    colnames(sub)[5:8]<-c('H0G1','H0U','H1G1','H1U')
    h0<-data.frame(pos=rep(sub$POS,2), prob=c(sub$H0G1, sub$H0U), type=c(rep('G1', dim(sub)[1]), rep('Unr', dim(sub)[1])) )
    h1<-data.frame(pos=rep(sub$POS,2), prob=c(sub$H1G1, sub$H1U), type=c(rep('G1', dim(sub)[1]), rep('Unr', dim(sub)[1])) )
    cols<-c('red','grey')
    cols2<-rep('white',2)
  }
  
  if (ng==2){
    colnames(sub)[5:10]<-c('H0G1', 'H0G2', 'H0U','H1G1', 'H1G2', 'H1U')
    h0<-data.frame(pos=rep(sub$POS,3), prob=c(sub$H0G1, sub$H0G2, sub$H0U), type=c(rep('G1', dim(sub)[1]), rep('G2', dim(sub)[1]), rep('Unr', dim(sub)[1])) )
    h1<-data.frame(pos=rep(sub$POS,3), prob=c(sub$H1G1, sub$H1G2, sub$H1U), type=c(rep('G1', dim(sub)[1]), rep('G2', dim(sub)[1]), rep('Unr', dim(sub)[1])) )
    cols<-c('red', 'blue','grey')
    cols2<-rep('white',3)
  }
  
  # subset ibd tracks
  ibd_t<-ibd[ibd$target==t,]
  
  # plot
  g0<-ggplot() +
    geom_point(data=h0, aes(x=pos, y=prob, col=type)) + 
    geom_line(data=h0, aes(x=pos, y=prob, col=type)) +
    scale_color_manual(values=cols) +
    theme_classic() +
    theme(axis.text.x = element_text(colour='white'),
          axis.title.x = element_text(colour='white'),
          axis.ticks.x = element_line(colour='white')) +
    labs(y='Haplotype 0 IBD probability', col='shared with')
  
  #+  geom_text(data=ibd_t, aes(x=start+(end-start)/2, y=1, label=Prob))
  
  g1<-ggplot() +
    geom_point(data=h1, aes(x=pos, y=prob, col=type)) + 
    geom_line(data=h1, aes(x=pos, y=prob, col=type)) +
    scale_color_manual(values=cols) +
    theme_classic() +
    labs(y='Haplotype 1 IBD probability', x='chromosome positions (bp)', col='shared with')
  
  g2<-ggplot() +
    geom_point(data=h1, aes(x=pos, y=prob, col=type)) + 
    scale_color_manual(values=cols2) +
    theme_classic() + labs(col='shared with') +
    geom_vline(xintercept = ibd_t$start) +
    geom_text(data=ibd_t, aes(x=start+(end-start)/2, y=0.65, label=Prob), size=4) +
    geom_text(data=ibd_t[ibd_t$length_CM>=3 & (ibd_t$Prob=='A' | ibd_t$Prob=='B'),], aes(x=start+(end-start)/2, y=0.35, label=paste0('l=',round(length_CM,1),'cM')), size=3) +
    theme(axis.text.y = element_text(colour='white'),
          axis.title.y = element_text(colour='white'),
          axis.ticks.y = element_line(colour='white'),
          axis.text.x = element_text(colour='white'),
          axis.title.x = element_text(colour='white'),
          axis.ticks.x = element_line(colour='white'),
          axis.line = element_line(colour='white'),
          legend.text = element_text(colour='white'),
          legend.title = element_text(colour='white'))
  
  #+ geom_text(data=ibd_t, aes(x=start+(end-start)/2, y=1, label=Prob))
  
  a<-plot_grid(g0, g2, g1, ncol=1, rel_heights = c(1,0.5,1))
  return(a)
}

targets<-unique(n[5:length(n)])




system('mkdir -p Plots/')
t='HG00544'
g<-plot_target('HG00544')
ggsave(paste0('Plots/thorin_IBD_plot.',t,'.png'), g, width=10, height=6)

#for (t in targets[1:100]){
#  g<-plot_target(t)
#  ggsave(paste0('Plots/thorin_IBD_plot.',t,'.png'), g, width=10, height=6)
#}



