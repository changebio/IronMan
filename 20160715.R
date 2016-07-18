require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(GenomicFeatures)
ucsc.prom<- promoters(txdb,upstream = 2000,downstream = 2000,columns=c("tx_id", "tx_name","gene_id"))
ucsc.prom<- keepSeqlevels(ucsc.prom,paste0("chr",c(1:22,"X")))
ucsc.prom<- unique(ucsc.prom)
CTSS.prom <- readRDS("data/CTSS_promoters.rds")
CTSS.prom <- unique(CTSS.prom)
only.prom<- ucsc.prom[countOverlaps(ucsc.prom,CTSS.prom)==1]
ucsc.gene<- genes(txdb)
only.prom<- only.prom[countOverlaps(only.prom,ucsc.gene)<2]
require(genomation)
k4me3.gro<-ScoreMatrixList(gro.seq,only.prom,weight.col = "V5")
k4me3.gro<- intersectScoreMatrixList(k4me3.gro)
k4me3.gro1<- lapply(k4me3.gro,rowSums)


a<- (k4me3.gro1[[6]]-k4me3.gro1[[5]])/(k4me3.gro1[[6]]+k4me3.gro1[[5]])
hist(a,xlab = "D score",main = "Histogram of D score in K562-unique promoter H3K4me3")
k4me3.smt.gro<- lapply(k4me3.smt,function(x)ScoreMatrixBin(gro.seq[[1]],x,bin.num = 50,weight.col = "V5"))
k4me3.smt.gro1<- sapply(k4me3.smt.gro,colMeans)