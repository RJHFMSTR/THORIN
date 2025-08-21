
# New clustering, handle better MZ twins to increase sample size
grouping_ped<-function(t, ped=data.frame(), mz){
	sex_target=sex$sex[sex$ID1==t]
	age_target=age$age[age$ID1==t]

	# Get IDs of PO relationships
	PO<-unique(c(as.character(rel$ID2[rel$ID1==t & rel$Kinship>=po_down & rel$Kinship<po_up & rel$IBS0<0.0012]),as.character(rel$ID1[rel$ID2==t & rel$Kinship>=po_down & rel$Kinship<po_up & rel$IBS0<0.0012])))
	
	if (length(PO)>0){
		dpo<-data.frame(PO=PO);
		dpo$sex<-sex$sex[match(dpo$PO, sex$ID1)]
		dpo$age<-age$age[match(dpo$PO, age$ID1)]
		dpo<-dpo[(dpo$age-age_target)>15,] # keeps only parents, remove kids (i.e younger).
	
		if (dim(dpo)[1]>0){
		
			# check for MZ among parents
			dpo$mz<-0; dpo$mz[dpo$PO %in% mz]<-1
			if (1 %in% dpo$mz){
				if (d$Kinship[(d$ID1==dpo$PO[dpo$mz==1][1] & d$ID2==dpo$PO[dpo$mz==1][2]) | (d$ID1==dpo$PO[dpo$mz==1][2] & d$ID2==dpo$PO[dpo$mz==1][1])] > 0.4){
					dpomz1<-dpo[dpo$mz==1,]; dpomz1<-dpomz1[1,] # keep only one MZ parent, it does not matter which one, same genetics ...
					dpo<-rbind(dpo[dpo$mz==0,], dpomz1)
				}
			}

			if (dim(dpo)[1]==1 & dpo$sex[1]=='Male'){ped<-data.frame(target=t, father=NA, mother=dpo$PO[1])} #duos father
			if (dim(dpo)[1]==1 & dpo$sex[1]=='Female'){ped<-data.frame(target=t, father=dpo$PO[1], mother=NA)} # duos mother
			if (dim(dpo)[1]==2 & (dpo$sex[1]!=dpo$sex[2])){ped<-data.frame(target=t, father=dpo$PO[dpo$sex=='Male'], mother=dpo$PO[dpo$sex=='Female'])} # trios, opposite sex
			if (dim(dpo)[1]>2){print(t)} # what is this?
			if (dim(dpo)[1]==2 & (dpo$sex[1]==dpo$sex[2])){print(t)} # what is this?
		}
	}
	return(ped)
}

## Old clustering
#grouping_ped<-function(t, ped=data.frame()){
#        sex_target=sex$sex[sex$ID1==t]
#        PO<-unique(c(as.character(rel$ID2[rel$ID1==t & rel$Kinship>=po_down & rel$Kinship<po_up & rel$IBS0<0.0012]),as.character(rel$ID1[rel$ID2==t & rel$Kinship>=po_down & rel$Kinship<po_up & rel$IBS0<0.0012])))
#        if (length(PO)==1){
#                age1<-age$age[age$ID1==PO[1]]
#                if (age$age[age$ID1==t]+15<age1){ # makes sure that the diff in age between parent-offspring in more than 15 years. This is stringent but avoids selecting sibs pairs as parent-offspring.
#                        parent1=PO[1]
#                        sex1=sex$sex[sex$ID1==PO[1]]
#                        if (sex1=='Male'){father=PO[1]; mother='NA'}
#                        if (sex1=='Female'){father='NA'; mother=PO[1]}
#                        ped<-data.frame(target=t, father=father, mother=mother)
#        }}
#        if (length(PO)==2){
#                age1<-age$age[age$ID1==PO[1]]
#                age2<-age$age[age$ID1==PO[2]]
#                if ((age$age[age$ID1==t]+15<age1) & (age$age[age$ID1==t]+15<age2)){
#                        sex1=sex$sex[sex$ID1==PO[1]]
#                        sex2=sex$sex[sex$ID1==PO[2]]
#                        if (sex1!=sex2){ # make sure that the two parent have not the same genetic sex
#                                Tmp<-data.frame(ID1=c(PO[1], PO[2]), sex=c(sex1, sex2))
#                                father=as.character(Tmp$ID1[Tmp$sex=='Male']); mother=as.character(Tmp$ID1[Tmp$sex=='Female'])
#                                ped<-data.frame(target=t, father=father, mother=mother)
#        }}}
#        return(ped)
#}




grouping_sibs<-function(s, d_sibs=data.frame()){

        sibs<-c(rel$ID2[rel$ID1==s], rel$ID1[rel$ID2==s])
        for (s2 in sibs){

                ages<-c(age$V2[age$V1==s], age$V2[age$V1==s2])
                if (max(ages)-min(ages)<15){# avoid picking up parent-offspring relationships...
                        d_sibs<-rbind(d_sibs, data.frame(target=s, sib=s2))
                }
        }
        return(d_sibs)
}









grouping_rel<-function(t, out=data.frame(), full_relatives=d){

        print(t)
                # 1.2 Grouping of second- and third-degree relatives

        g1<-""; g2<-""
        relatives<-unique(c(as.character(rel$ID2[rel$ID1==t & rel$Kinship<0.1767]),as.character(rel$ID1[rel$ID2==t & rel$Kinship<0.1767])))


        # ** small update: for second degree relatives, use only older ones (likely uncle/aunt or grandparents)
        for (r in relatives){
                k<-c(rel$InfType[rel$ID1==t & rel$ID2==r], rel$InfType[rel$ID2==t & rel$ID1==r])
                if (k=='2nd'){ # 2nd degree relative --> check diff in age.
                        age_t<-age$age[age$ID1==t]
                        age_r<-age$age[age$ID1==r]
                        if (age_t>(age_r+20)){ # Target is still older than the relative+20, relative likely an nephew. avoid using the nephew. Remove only extreme cases. Might still have some error here, but it's impossible to do that perfectly.
                                relatives<-relatives[relatives!=r]
                        }
                }
        }


        # ** end update





        if (length(relatives)!=0){

                d_pairs<-data.frame()
                for (relat in relatives){
                        relatives2<-unique(c(as.character(full_relatives$ID2[full_relatives$ID1==relat]),as.character(full_relatives$ID1[full_relatives$ID2==relat])))
                        relatives2<-relatives2[relatives2 %in% relatives]; relatives2<-relatives2[relatives2!=t]
                        d_pairs<-rbind(d_pairs, data.frame(t1=rep(relat, length(relatives2)), t2=relatives2))
                        if (length(relatives2)==0){d_pairs<-rbind(d_pairs, data.frame(t1=relat, t2=relat)) }}

                m<-as.matrix(d_pairs); g<-simplify(graph_(m, from_edgelist(), directed = FALSE));
                dg <- decompose.graph(g)

                p='T'
                if (length(dg)==0){p='F'}
                if (length(dg)==1) { g1=vertex_attr(dg[[1]])$name; g2=""}
                if (length(dg)==2){ g1=vertex_attr(dg[[1]])$name; g2=vertex_attr(dg[[2]])$name}
                if (length(dg)>2){# more than two groups. keep only one relative. Preferentially 2nd degree
                        r2<-rbind(rel[rel$ID1==t,], rel[rel$ID2==t,])
                        r2<-r2[order(-r2$Kinship),]
                        g1<-c(r2$ID1[1], r2$ID2[1]); g1<-g1[g1!=t]
                }


                if (g1[1]=="" & g2[1]==""){ p='F'}
                if (g1[1]!="" & g2[1]==""){ G1<-paste0("G1=",paste0(g1, collapse=';')); G2=""; out<-data.frame(target=t,G1=G1, G2=G2)}
                if (g1[1]!="" & g2[1]!="") { G1<-paste0("G1=",paste0(g1, collapse=';'));  G2<-paste0("G2=",paste0(g2, collapse=';')); out<-data.frame(target=t,G1=G1, G2=G2)}


                # this part can be used to produce a plot of the close relatives "network"
#               d_pairs2<-d_pairs
#               for (relat in relatives){d_pairs2<-rbind(d_pairs2, data.frame(t1=t, t2=relat))}
#
#               m2<-as.matrix(d_pairs2); g2<-simplify(graph_(m2, from_edgelist(), directed = FALSE))
#               V(g2)$color[names(V(g2))==t]<-'green'
#               png(paste0(t,'_grouping.png'))
#               plot(g2)
#               dev.off()
        }

        return(out)

}















get_side<-function(t){
        sub<-rel[rel$V1==t,]
        father<-ped$V2[ped$V1==t]
        mother<-ped$V3[ped$V1==t]
        
        g1=unlist(strsplit(gsub('G1=','',sub$V2),';'))
        g2=unlist(strsplit(gsub('G2=','',sub$V3),';'))
        
        tmp<-data.frame()
        for (r in g1){
                relatives<-unique(c(as.character(d$ID2[d$ID1==r]),as.character(d$ID1[d$ID2==r])))
                if (mother %in% relatives){tmp<-rbind(tmp, data.frame(target=t, relative=r, group='G1',side='maternal'))}
                if (father %in% relatives){tmp<-rbind(tmp, data.frame(target=t, relative=r, group='G1',side='paternal'))}
        }
        
        if (length(g2)!=0){
                for (r in g2){
                        relatives<-unique(c(as.character(d$ID2[d$ID1==r]),as.character(d$ID1[d$ID2==r])))
                        if (mother %in% relatives){tmp<-rbind(tmp, data.frame(target=t, relative=r, group='G2',side='maternal'))}
                        if (father %in% relatives){tmp<-rbind(tmp, data.frame(target=t, relative=r, group='G2',side='paternal'))}
                }
        }
        return(tmp)
}


