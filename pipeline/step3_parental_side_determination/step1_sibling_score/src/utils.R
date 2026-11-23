#setwd('~/Dropbox/Ressources/Git_repository/THORIN/pipeline/step3_parental_side_determination/step1_sibling_score/src/')
library(ggplot2)
library(cowplot)
library(parallel)
library(dplyr)
library(data.table)






### *** define functions.

## ---------------------------------------------------------------------------
## ibd2co: call crossovers on each target haplotype from a one-group THORIN run
## (G1 = sibling), using the IBD labels written by THORIN since commit 707b325
## (January 2025; v1.2.0 release binary and v1.3):
##   A = only h0 shares with the group   (H0G * H1U)
##   B = only h1 shares with the group   (H1G * H0U)
##   C = neither haplotype shares
##   D = both haplotypes share (IBD2)
## Sharing of h0 changes only at crossovers on h0's parental side, and sharing
## of h1 only at crossovers on h1's side. So we track each haplotype on its own
## (h0 shared = A or D; h1 shared = B or D) and call a crossover wherever that
## haplotype's own sharing changes. An A->D boundary is therefore a crossover
## on h1, not h0.
##
## Short runs: on each haplotype track, consecutive segments with the same
## state are merged into runs. Runs shorter than `ibd_length` cM are treated
## as noise and removed, and the remaining neighbouring runs with the same
## state are merged again. A crossover is called between two consecutive
## remaining runs with different states, at the midpoint between the end of
## the left run and the start of the right run. Chromosome ends are never
## counted as crossovers.
##
## Returns data.frame(target, CHR, bp, haplotype) with haplotype in H0/H1,
## or NULL if no crossover is found (same contract as the previous version).
## ---------------------------------------------------------------------------

segments_to_runs<-function(start, end, length_cm, state){
  r<-rle(state)
  idx_end<-cumsum(r$lengths)
  idx_start<-c(1, head(idx_end, -1) + 1)
  data.frame(state=r$values,
             start=start[idx_start],
             end=end[idx_end],
             length_CM=mapply(function(a, b) sum(length_cm[a:b]), idx_start, idx_end))
}

merge_runs<-function(runs){
  if (nrow(runs) < 2) return(runs)
  r<-rle(runs$state)
  idx_end<-cumsum(r$lengths)
  idx_start<-c(1, head(idx_end, -1) + 1)
  data.frame(state=r$values,
             start=runs$start[idx_start],
             end=runs$end[idx_end],
             length_CM=mapply(function(a, b) sum(runs$length_CM[a:b]), idx_start, idx_end))
}

runs_to_crossovers<-function(runs, ibd_length){
  runs<-runs[runs$length_CM >= ibd_length, , drop=FALSE]   # drop short runs (noise)
  if (nrow(runs) < 2) return(numeric(0))
  runs<-merge_runs(runs)                                    # re-merge same-state neighbours
  if (nrow(runs) < 2) return(numeric(0))
  left<-runs[-nrow(runs), ]; right<-runs[-1, ]
  round((left$end + right$start) / 2)
}

ibd2co<-function(data, t, ibd_length, CHR){
  tmp<-data[data$target==t & data$CHR %in% c(paste0('chr',CHR), as.character(CHR)), ]
  tmp<-tmp[!is.na(tmp$Prob) & tmp$Prob %in% c('A','B','C','D'), ]
  if (nrow(tmp) < 2) return(NULL)
  tmp<-tmp[order(tmp$start), ]
  res<-list()
  for (h in c('H0','H1')){
    shared<-if (h=='H0') tmp$Prob %in% c('A','D') else tmp$Prob %in% c('B','D')
    runs<-segments_to_runs(tmp$start, tmp$end, tmp$length_CM, shared)
    co<-runs_to_crossovers(runs, ibd_length)
    if (length(co) > 0) res[[h]]<-data.frame(target=t, CHR=CHR, bp=co, haplotype=h)
  }
  if (length(res) == 0) return(NULL)
  res<-do.call(rbind, res); res<-res[order(res$bp), ]; rownames(res)<-NULL
  return(res)
}


## Legacy version, kept only to reproduce results from THORIN builds that
## predate commit 707b325 (January 2025). In those builds, a one-group run read
## misaligned columns, which made A ~ "h0 shared" and B ~ "h1 shared"; this
## function relies on that behaviour. Do NOT use it with the v1.2.0 release
## binary or with v1.3.
ibd2co_v12<-function(data, t, ibd_length, CHR){
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




