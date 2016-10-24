### =========================================================================
### peakState()
### -------------------------------------------------------------------------
##' GO Enrichment Analysis of a gene set.
##' Given a vector of genes, this function will return the enrichment GO
##' categories after FDR control.
##'
##'
##' @param gene a vector of entrez gene id.
##' @param OrgDb OrgDb
##' @param keytype keytype of input gene
##' @param ont One of "MF", "BP", and "CC" subontologies.
##' @param pvalueCutoff Cutoff value of pvalue.
##' @param pAdjustMethod one of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
##' @param universe background genes
##' @param qvalueCutoff qvalue cutoff
##' @param minGSSize minimal size of genes annotated by Ontology term for testing.
##' @param maxGSSize maximal size of genes annotated for testing
##' @param readable whether mapping gene ID to gene Name
##' @return A \code{enrichResult} instance.
##' @importClassesFrom DOSE enrichResult
##' @importMethodsFrom DOSE show
##' @importMethodsFrom DOSE summary
##' @importMethodsFrom DOSE plot
##' @importFrom DOSE setReadable
##' @seealso \code{\link{enrichResult-class}}, \code{\link{compareCluster}}
##' @keywords manip
##' @export
##' @author Guangchuang Yu \url{http://ygc.name}
##' @examples
##' \dontrun{
##' 	data(gcSample)
##' 	yy <- enrichGO(gcSample[[1]], 'org.Hs.eg.db', ont="BP", pvalueCutoff=0.01)
##' 	head(summary(yy))
##' 	plot(yy)
##' }

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


