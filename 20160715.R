#Dscore in different genomic features

#-----
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(GenomicFeatures)

ucsc.prom<- promoters(txdb,upstream = 2000,downstream = 2000,columns=c("tx_id", "tx_name","gene_id"))
ucsc.prom<- keepSeqlevels(ucsc.prom,paste0("chr",c(1:22,"X")))
ucsc.prom<- unique(ucsc.prom)
CTSS.prom <- readRDS("data/CTSS_promoters.rds")
CTSS.prom <- unique(CTSS.prom)
only.prom<- ucsc.prom[countOverlaps(ucsc.prom,CTSS.prom,ignore.strand=TRUE)==1]
ucsc.gene<- genes(txdb)




only.prom<- only.prom[countOverlaps(only.prom,ucsc.gene,ignore.strand=TRUE)<2]
require(genomation)
Dscore<- function(prom){
  tss.sml<-ScoreMatrixList(gro.seq,prom,weight.col = "V5")
  prom <- prom[as.numeric(Reduce(intersect,lapply(tss.sml,rownames)))]
  tss.sml<- intersectScoreMatrixList(tss.sml)
  tss.rs<- lapply(tss.sml,rowSums)
  prom$wTAP <- (tss.rs[[2]]-tss.rs[[1]])/(tss.rs[[2]]+tss.rs[[1]])
  prom$nTAP <- (tss.rs[[4]]-tss.rs[[3]])/(tss.rs[[4]]+tss.rs[[3]])
  prom$cage <- (tss.rs[[6]]-tss.rs[[5]])/(tss.rs[[6]]+tss.rs[[5]])
  prom$Cage <- NA
  prom$Cage[tss.rs[[6]]+tss.rs[[5]]>=6] <- ((tss.rs[[6]]-tss.rs[[5]])/(tss.rs[[6]]+tss.rs[[5]]))[tss.rs[[6]]+tss.rs[[5]]>=6]
  prom$TAP <- (tss.rs[[4]]-tss.rs[[2]]-tss.rs[[3]]+tss.rs[[1]])/(tss.rs[[4]]+tss.rs[[2]]-tss.rs[[3]]-tss.rs[[1]])
  mcols(prom) <- cbind(mcols(prom),DataFrame(tss.rs))
  return(prom)
}


#Dscore distribution in region 500bp-----------
ucsc.gene.pt<- promoters(ucsc.gene,upstream = 250,downstream = 250)
gro.ds<- Dscore(ucsc.gene.pt)

temp<- as.data.frame(only.prom)
View(temp)
a<- (k4me3.gro1[[6]]-k4me3.gro1[[5]])/(k4me3.gro1[[6]]+k4me3.gro1[[5]])
hist(a,xlab = "D score",main = "Histogram of D score in K562-unique promoter H3K4me3")
k4me3.smt.gro<- lapply(k4me3.smt,function(x)ScoreMatrixBin(gro.seq[[1]],x,bin.num = 50,weight.col = "V5"))
k4me3.smt.gro1<- sapply(k4me3.smt.gro,colMeans)