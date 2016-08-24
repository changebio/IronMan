#Dscore in different genomic features

#load packages-----
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(GenomicFeatures)
require(genomation)

#gene signal calculation
##gro and cage signal distribution in promoters
ucsc.gene<- genes(txdb)
ucsc.gene.pt<- promoters(ucsc.gene,upstream = 250,downstream = 250)
gro.sg.m<-lapply(gro.seq,function(x)ScoreMatrixBin(x,ucsc.gene.pt[strand(ucsc.gene.pt)=="+"],bin.num = 50,weight.col = "V5"))
gro.sg.p<-lapply(gro.seq,function(x)ScoreMatrixBin(x,ucsc.gene.pt[strand(ucsc.gene.pt)=="-"],bin.num = 50,weight.col = "V5"))

gro.sg<- melt(rbind(sapply(gro.sg.p,colMeans),sapply(gro.sg.m,colMeans)))
gro.sg$gene<- c(rep("+",50),rep("-",50))
gro.sg$read<- c(rep("minus",100),rep("plus",100))
gro.sg$type <- c(rep("wTAP",200),rep("nTAP",200),rep("CAGE",200))
gro.sg$Var1<- c(-24.5:24.5)*10

ggplot(gro.sg)+geom_line(data=gro.sg[gro.sg$gene=="+",],aes(x=Var1,y=value,colour=gene,linetype=as.factor(read)))+
  geom_line(data=gro.sg[gro.sg$gene=="-",],aes(x=Var1,y=-value,colour=gene,linetype=as.factor(read)))+
  facet_grid(. ~ type)+
  labs(x = "",y = " ",title="The average signal of Promoters with 500bp width in K562",linetype="read") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size = 16)
  )

###cage from ENCODE
cage.seq<- readRDS("/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/data/cage_gr.rds")
cage.sg.m<-lapply(cage.seq,function(x)ScoreMatrixBin(x,ucsc.gene.pt[strand(ucsc.gene.pt)=="+"],bin.num = 50,weight.col = "X"))
cage.sg.p<-lapply(cage.seq,function(x)ScoreMatrixBin(x,ucsc.gene.pt[strand(ucsc.gene.pt)=="-"],bin.num = 50,weight.col = "X"))

cage.sg<- melt(rbind(sapply(cage.sg.p,colMeans),sapply(cage.sg.m,colMeans)))
cage.sg$gene<- c(rep("+",50),rep("-",50))
cage.sg$read<- c(rep("minus",100),rep("plus",100))
cage.sg$type <- substring(cage.sg$Var2,first = 1,last = nchar(as.character(cage.sg$Var2))-c(rep(6,100),rep(5,100)))
#cage.sg$region<- c(rep("cell",400),rep("chromatin",200),rep("cytosol",800),rep("nucleolus",200),rep("nucleoplasm",200),rep("nucleus",600),rep("polysome",200))
cage.sg$Var1<- c(-24.5:24.5)*10


gro.cage.sg<- rbind(gro.sg,cage.sg)
ggplot(gro.cage.sg)+geom_line(data=gro.cage.sg[gro.cage.sg$gene=="+",],aes(x=Var1,y=value,colour=gene,linetype=as.factor(read)))+
  geom_line(data=gro.cage.sg[gro.cage.sg$gene=="-",],aes(x=Var1,y=-value,colour=gene,linetype=as.factor(read)))+
  facet_wrap( ~ type)+
  labs(x = "",y = " ",title="The average signal of Promoters with 500bp width in K562",linetype="read") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size = 16)
  )
##other signal distribution (H3K4me3,H3K4me1,....)



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
  labs(x = "",y = " ",title="Histogram of D score of Promoters with 500bp width in K562") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size=16)
  )



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
