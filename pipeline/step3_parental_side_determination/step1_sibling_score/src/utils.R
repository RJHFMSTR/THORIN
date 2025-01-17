#setwd('~/Dropbox/Ressources/Git_repository/THORIN/pipeline/step3_parental_side_determination/step1_sibling_score/src/')
library(ggplot2)
library(cowplot)
library(parallel)
library(dplyr)
library(data.table)






### *** define functions.

ibd2co<-function(data, t, ibd_length, CHR){
  tmp<-data[data$CHR==paste0('chr',CHR),]
  tmp<-tmp[tmp$target==t & tmp$length_CM>=ibd_length & (tmp$Prob=='A' | tmp$Prob=='B'),]
  if (dim(tmp)[1]!=0){
    tmp$target_haplotype<-'H0'
    tmp$target_haplotype[tmp$Prob=='B']<-'H1'
    res<-data.frame(target=t, CHR=CHR, bp=c(tmp$start, tmp$end), haplotype=rep(tmp$target_haplotype,2))
    return(res)
  }
}


read_rec_map<-function(CHR, sex){
  tmp<-read.table(paste0('data/Refined_EUR_genetic_map_b37/',sex,'_chr',CHR,'.txt'), hea=T)
  tmp$chr<-CHR; tmp$Morgan<-tmp$cM/100
  return(tmp)
}


bp_window<-function(bps, window_size, map){
  window_down<-bps-window_size
  window_up<-bps+window_size
  window_down[window_down<min(map$pos)]<-min(map$pos); window_down[window_down>max(map$pos)]<-max(map$pos)
  window_up[window_up<min(map$pos)]<-min(map$pos); window_up[window_up>max(map$pos)]<-max(map$pos)
  return(list(window_down, window_up))
}



bp2morgan<-function(bps, map){
  return(approx(x = map$pos, y = map$Morgan, xout = bps, rule = 2)$y)
}


convert2log<-function(array){
  v<-unlist(lapply(array, log10))
  return(sum(v))
}




chr2expectedCO<-function(CHR, map_males, map_females){
  return(mean(c(max(map_males$Morgan), max(map_females$Morgan))))
}


find_outliers<-function(data, value){
  q1=quantile(data,0.25)
  q3=quantile(data,0.75)
  IQR=q3-q1
  UB=q3+(value*IQR)
  LB=q1-(value*IQR)
  return(data[data>UB | data<LB])
}




