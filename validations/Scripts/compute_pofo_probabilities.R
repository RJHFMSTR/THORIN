library(parallel)
library(dplyr)
args=commandArgs(trailingOnly=T)

# ***** Script parameters
N_CORE=2
ID=00
prob_file=paste0('Out/step2_thorin_autosomes/split/chr22/ukb.chr22.thorin.subset_',ID,'.prob.gz')
group_file=paste0('data/split/Relatives.maternal_side.subset_',ID,'.group')
prob_out='prob_out.bed'


prob_file=args[1]
group_file=args[2]
prob_out=args[3]
CHR=args[4]
N_CORE=args[5]


cat('-----------------------\n')
cat(paste0('THORIN file: ',prob_file, '\n'))
cat('-----------------------\n')


group<-as.data.frame(data.table::fread(group_file, hea=F, fill=T, sep='\t')); group<-group[,1:3]
targets=group$V1




# ***** Function to get consecutive sequence of the same item in a vector
find_consecutive_sequences <- function(vector, item) {
  sequences <- data.frame(start = integer(), end = integer(), stringsAsFactors = FALSE)
  start <- NULL

  for (i in seq_along(vector)) {
    if (vector[i] == item) {
      if (is.null(start)) {
        start <- i
      }
    } else if (!is.null(start)) {
      sequences <- rbind(sequences, c(start, i - 1))
      start <- NULL
    }
  }

  if (!is.null(start)) {
    sequences <- rbind(sequences, c(start, length(vector)))
  }

  colnames(sequences) <- c("start", "end")
  return(sequences)
}
# *************************************







# ***** Function to identify IBD tracks and build a scaffold for the phasing (3rd degree relatives)
get_prob_tracks<-function(t){
#	print(t)
	COLS_BASE=c('#CHROM','POS','CM')
	t_group<-group[group$V1==t,]
	t_group[is.na(t_group)]<-'' # solve cases when G2 is NA


	# check the number of group. Adapt if n=1 and group=G2
	
        	# G1 and not G2
	if (t_group$V2!='' & t_group$V3==''){
		n_group=1
		COLS<-c(COLS_BASE, paste0(t,'_G1_0'),paste0(t,'_HOLE_0'),paste0(t,'_G1_1'),paste0(t,'_HOLE_1'))
		name_cols<-c('CHR','POS','CM','H1G1','H1U','H2G1','H2U')
	
		d<-as.data.frame(data.table::fread(prob_file, hea=T, select=COLS))
        	colnames(d)<-name_cols
		d$H1G2<-0; d$H2G2<-0
	}

	# G2 and not G1
       if (t_group$V2=='' & t_group$V3!=''){
                n_group=1
                COLS<-c(COLS_BASE, paste0(t,'_G2_0'),paste0(t,'_HOLE_0'),paste0(t,'_G2_1'),paste0(t,'_HOLE_1'))
                name_cols<-c('CHR','POS','CM','H1G2','H1U','H2G2','H2U')
        
                d<-as.data.frame(data.table::fread(prob_file, hea=T, select=COLS))
                colnames(d)<-name_cols
                d$H1G1<-0; d$H2G1<-0
        }


       # G1 and G2
       if (t_group$V2!='' & t_group$V3!=''){
                n_group=1
		COLS<-c(COLS_BASE, paste0(t,'_G1_0'),paste0(t,'_G2_0'),paste0(t,'_HOLE_0'),paste0(t,'_G1_1'),paste0(t,'_G2_1'),paste0(t,'_HOLE_1'))
		name_cols<-c('CHR','POS','CM','H1G1','H1G2','H1U','H2G1','H2G2','H2U')
                d<-as.data.frame(data.table::fread(prob_file, hea=T, select=COLS))
                colnames(d)<-name_cols
        }
	



	## *** Prob tracks detection:

	# Get probabilities of H1-G1 and H2-G2

        d$PrA<-d$H1G1*d$H2G2 + d$H1G1*d$H2U + d$H1U*d$H2G2
        d$PrB<-d$H1G2*d$H2G1 + d$H1G2*d$H2U + d$H2G1*d$H1U
        d$PrC<-d$H1U*d$H2U
	d$PrD<-d$H1G1*d$H2G1 + d$H1G2*d$H2G2
	




	# Normalize the probabilities
	d$PrAnorm<-d$PrA/(d$PrA+d$PrB+d$PrC+d$PrD)
	d$PrBnorm<-d$PrB/(d$PrA+d$PrB+d$PrC+d$PrD)
	d$PrCnorm<-d$PrC/(d$PrA+d$PrB+d$PrC+d$PrD)
	d$PrDnorm<-d$PrD/(d$PrA+d$PrB+d$PrC+d$PrD)
	
	# Find the max probability at each position
	m<-as.data.frame(as.matrix(d[colnames(d) %in% c('PrAnorm','PrBnorm','PrCnorm','PrDnorm')]))
	pr<-c('A','B','C','D')
	d$max_pr<-unlist(apply(m, 1, FUN=function(x){pr[which.max(x)]}))
	
	# Probabilities:
	# A --> H1=G1 (and H2=G2, if G2 used)
        # B --> H2=G1 (and H1=G2, if G2 used)
        # C --> H1=Unr and H2=Unr
        # D --> H1=G1 and H2=G1 (or H1=G2 and H2=G2, if G2 used; correspond to a disomy event, or a deletion)


	# Find consecutive sequence of the same probability. Get the length in CM of that sequence. We will then use the length [CM] to keep only large IBD tracks.
	SEQ<-data.frame()
	for (pr in c('A','B','C','D')){
	seq<-find_consecutive_sequences(d$max_pr, pr); if (dim(seq)[1]>0){seq$CM_start=d$CM[seq$start]; seq$CM_end=d$CM[seq$end]; seq$CM_length<-seq$CM_end-seq$CM_start; seq$POS_start=d$POS[seq$start]; seq$POS_end=d$POS[seq$end]; seq$Pr<-pr; SEQ<-rbind(SEQ, seq)}
	}
	SEQ$target<-t
	SEQ<-SEQ[order(SEQ$start),]



	#*********************************************************
        # ADD THIS TO CATEGORIZE UPD INTO MAT AND PAT
	SEQ$UPD<-'U'
	d$PrD_G1<-d$H1G1*d$H2G1
        d$PrD_G2<-d$H1G2*d$H2G2
	if ('D' %in% SEQ$Pr){
		tmp<-SEQ[SEQ$Pr=='D',]
		for (i in 1:dim(tmp)[1]){
			print(tmp[i,])
			g1_upd_prob=mean(d$PrD_G1[d$POS>=tmp$POS_start[i] & d$POS<=tmp$POS_end[i]])
                        g2_upd_prob=mean(d$PrD_G2[d$POS>=tmp$POS_start[i] & d$POS<=tmp$POS_end[i]])
			if (g1_upd_prob>g2_upd_prob){SEQ$UPD[SEQ$POS_start==tmp$POS_start[i]]<-'G1'}
                        if (g1_upd_prob<g2_upd_prob){SEQ$UPD[SEQ$POS_start==tmp$POS_start[i]]<-'G2'}

		}
	
	}
        #*********************************************************

	return(SEQ)

}



# run in parrellel for each target-relative
x<-mclapply(targets, get_prob_tracks, mc.cores=N_CORE)

# transform the results into a bed file
xx<-cbind(bind_rows(x))


#d<-data.frame(CHR=CHR, start=xx$POS_start, end=xx$POS_end, Prob=xx$Pr, length_CM=xx$CM_length, target=xx$target)
#write.table(d, prob_out, quote=F, col.names=T, row.names=F, sep='\t')


d<-data.frame(CHR=CHR, start=xx$POS_start, end=xx$POS_end, Prob=xx$Pr, length_CM=xx$CM_length, target=xx$target, UPD=xx$UPD)
write.table(d, prob_out, quote=F, col.names=T, row.names=F, sep='\t')






