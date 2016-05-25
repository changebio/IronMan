
Peaks.Distribution<- function(peaks,refhg19=refhg19){
  peaks$state<- "Intergenic"
  tp.promoter<- overlap(peaks,refhg19$PROMOTER)
  peaks$state[tp.promoter[[1]]>0]<- "Promoter"
  tp.exon<- overlap(peaks[peaks$state=="Intergenic",],refhg19$EXON)
  peaks$state[peaks$state=="Intergenic"][tp.exon[[1]]>0]<- "Exon"
  tp.intron<- overlap(peaks[peaks$state=="Intergenic",],refhg19$INTRON)
  peaks$state[peaks$state=="Intergenic"][tp.intron[[1]]>0]<- "Intron"
  return(peaks)
}


Toppeaks<- function(x,top=20000){
  pk<- x[order(-x[,7]),]
  return(pk[1:min(top,nrow(x)),])
}


#peak targeted by genes
overlap.target.exp<- function(rbed,cbed,col=4,exp,lap=TRUE){
  rlap<- rep(NA,nrow(rbed))
  rlap.list<- rep(NA,nrow(rbed))
  temp<- rep(NA,nrow(rbed))
  for(chr in intersect(levels(rbed[,1]),levels(cbed[,1]))){
    rind<- which(rbed[,1]==chr)
    cind<- which(cbed[,1]==chr)
    cname<- cbed[cind,col]
    cexp<- exp[cind]
    if(length(rind)>0 && length(cind)>0){
      rmts<- rep.col(as.integer(rbed[rind,2]),length(cind))
      rmte<- rep.col(as.integer(rbed[rind,3]),length(cind))
      cmts<- rep.row(as.integer(cbed[cind,2]),length(rind))
      cmte<- rep.row(as.integer(cbed[cind,3]),length(rind))
      if(lap){
        s<- rmts<=cmte
        e<- rmte>=cmts
        se<- s==e
      }else{
        s<- rmts<=cmts
        e<- rmte>=cmte
        se<- s==e
      }
      
      rlap[rind]<- apply(se, 1, function(x){as.character(cname[x>0][which.max(cexp[x>0])])})
      rlap.list[rind]<- apply(se, 1, function(x){as.character(cname[x>0])})
      temp[rind]<- apply(se, 1, function(x){cexp[x>0]})
    }
  }
  return(list(rlap,rlap.list,temp))
}

overlap.target.dist<- function(rbed,cbed,col=4,lap=TRUE){
  rlap<- rep(NA,nrow(rbed))
  rlap.list<- rep(NA,nrow(rbed))
  summit<- rowMeans(rbed[,2:3])
  tss<- rowMeans(cbed[,2:3])
  temp<- rep(NA,nrow(rbed))
  for(chr in intersect(levels(rbed[,1]),levels(cbed[,1]))){
    rind<- which(rbed[,1]==chr)
    cind<- which(cbed[,1]==chr)
    cname<- cbed[cind,col]
    ctss<- tss[cind]
    if(length(rind)>0 && length(cind)>0){
      rmts<- rep.col(as.integer(rbed[rind,2]),length(cind))
      rmte<- rep.col(as.integer(rbed[rind,3]),length(cind))
      cmts<- rep.row(as.integer(cbed[cind,2]),length(rind))
      cmte<- rep.row(as.integer(cbed[cind,3]),length(rind))
      if(lap){
        s<- rmts<=cmte
        e<- rmte>=cmts
        se<- s==e
      }else{
        s<- rmts<=cmts
        e<- rmte>=cmte
        se<- s==e
      }
      
      rlap[rind]<- lapply(1:nrow(se),function(i){as.character(cname[se[i,]>0])[which.min(abs(ctss[se[i,]>0]-summit[rind[i]]))]})
      rlap.list[rind]<- apply(se, 1, function(x){as.character(cname[x>0])})
      temp[rind]<- lapply(1:nrow(se),function(i){abs(ctss[se[i,]>0]-summit[rind[i]])})
    }
  }
  return(list(rlap,rlap.list,temp))
}
#genes targeted by peaks