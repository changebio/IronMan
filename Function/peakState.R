### =========================================================================
### peakState()
### -------------------------------------------------------------------------


peakStatebylap <- function(query,malist,subjlist){
  if(!is.null(malist)){
    qm<-sapply(malist,function(x){
      qm<- rep(NA,length(query))
      qs.hit<- as.data.frame(GenomicRanges::findOverlaps(query,x))
      qm[qs.hit[,1]]<- x$M_value[qs.hit[,2]]
      return(qm)
    })
  }
  if(!is.null(subjlist)){
    qs<-sapply(subjlist,function(x)GenomicRanges::countOverlaps(query,x)>0)
  }
  
  GenomicRanges::mcols(query)<- cbind(GenomicRange::mcols(query),DataFrame(qm,qs))
  return(query)
}

#tpp<- peakStatebylap(k562.h1.k4me3.manorm,c(H3K4me1=k562.h1.k4me1.manorm,H3K27ac=k562.h1.k27ac.manorm),c(K562SE=K562.SE,H1SE=H1.SE))
#tpp


