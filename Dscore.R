#Dscore in different genomic features

#load packages-----
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(GenomicFeatures)
require(genomation)

#Dscore function-----
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
ucsc.gene<- genes(txdb)
ucsc.gene.pt<- promoters(ucsc.gene,upstream = 250,downstream = 250)
gro.ds<- Dscore(ucsc.gene.pt)

gene.ds <- as.data.frame(gro.ds)

gene.ds.gp<-melt(gene.ds[,c(5,7:9)])
ggplot(gene.ds.gp)+geom_histogram(data=gene.ds.gp[gene.ds.gp$strand=="+",],aes(x=value,y=..count..,fill=strand))+
  geom_histogram(data=gene.ds.gp[gene.ds.gp$strand=="-",],aes(x=value,y=-..count..,fill=strand))+
  facet_grid(. ~ variable)+
  labs(x = "",y = " ",title="the distribution of Dscore") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size=16)
  )

ggplot(gene.ds.gp[gene.ds.gp$strand=="+",])+geom_histogram(aes(x=value,y=..count..,fill=strand))+
  geom_histogram(data=gene.ds.gp[gene.ds.gp$strand=="-",],aes(x=value,y=-..count..,fill=strand))

apply(as.data.frame(mcols(gro.ds)[,2]),2,function(x)print(hist(x,main = "Histogram of D score of genes in K562 GROcap wTAP")))
apply(as.data.frame(mcols(gro.ds)[,3]),2,function(x)print(hist(x,main = "Histogram of D score of genes in K562 GROcap nTAP")))
apply(as.data.frame(mcols(gro.ds)[,4]),2,function(x)print(hist(x,main = "Histogram of D score of genes in K562 CAGE")))


grl.k4me3.ma.c <- lapply(k4me3.smt,function(x)resize(x,width = 240,fix="center"))
grl.k4me3.ma.ds<- lapply(grl.k4me3.ma.c,Dscore)
lapply(grl.k4me3.ma.ds,function(x)print(hist(as.numeric(mcols(x)[,2]),main = "Histogram of D score of genes in K562 GROcap wTAP")))
temp<-sapply(grl.k4me3.ma.ds,function(x)as.numeric(mcols(x)[,3]))
temp<- melt(temp)
ggplot(temp,aes(value))+geom_histogram()+facet_grid(.~L1)+
  labs(x = "",y = " ",title="The Dscore of H3K4me3 in K562 CAGE",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )

ggplot(temp,aes(value))+geom_density()+facet_grid(.~L1)

ucsc.prom<- promoters(txdb,upstream = 2000,downstream = 2000,columns=c("tx_id", "tx_name","gene_id"))
ucsc.prom<- keepSeqlevels(ucsc.prom,paste0("chr",c(1:22,"X")))
ucsc.prom<- unique(ucsc.prom)
CTSS.prom <- readRDS("data/CTSS_promoters.rds")
CTSS.prom <- unique(CTSS.prom)
only.prom<- ucsc.prom[countOverlaps(ucsc.prom,CTSS.prom,ignore.strand=TRUE)==1]
only.prom<- only.prom[countOverlaps(only.prom,ucsc.gene,ignore.strand=TRUE)<2]




only.prom<- only.prom[countOverlaps(only.prom,ucsc.gene,ignore.strand=TRUE)<2]
