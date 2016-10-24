### =========================================================================
### peakState()
### -------------------------------------------------------------------------
#' this is a silly function
#' @param input1 this is an input to our function
#' @param input2 this is another input
#' @return some value
#' @export

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


setMAvalbylap<- function(query,masub,type="M_value"){
  qm <- rep(NA,length(query))
  qs.hit <- as.data.frame(GenomicRanges::findOverlaps(query,masub))
  if(is.factor(mcols(masub)[[type]])){
    qm[qs.hit[,1]] <- as.character(mcols(masub)[[type]])[qs.hit[,2]]
  }else{
    qm[qs.hit[,1]] <- mcols(masub)[[type]][qs.hit[,2]]
  }
  return(qm)
}

setMAvalbylap(k562.h1.k4me3.manorm,k562.h1.k4me1.manorm)

setMtofactors<- function(mval){
  m<- rep("M",length(mval))
  m[mval>=1]<- "B"
  m[-1>=mval]<- "S"
  m[is.na(mval)]<- "N"
  return(m)
}