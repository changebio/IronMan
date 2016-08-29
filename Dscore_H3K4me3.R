##load packages------------------
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(GenomicRanges)
require(ChIPseeker)
#preprocess H3K4me3 after running MAnorm by K562 and H1======
k562.h1.k4me3.manorm<- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me3/K562.H1.rep1.top2k.20160520/K562.H1.rep1.top2k.20160520_all_peak_MAvalues.xls",header = TRUE)

CTSS.prom <- readRDS("data/CTSS_promoters.rds")
H1.H3K4me3.state<- function(k562.h1.k4me3.manorm){
  k562.h1.k4me3.manorm<- dropSeqlevels(k562.h1.k4me3.manorm,"chrY")
  k562.h1.k4me3.manorm$K562<- k562.h1.k4me3.manorm$M_value>=-1
  k562.h1.k4me3.manorm$H1hesc<- k562.h1.k4me3.manorm$M_value<=1
  k562.h1.k4me3.manorm$Promoter<-countOverlaps(k562.h1.k4me3.manorm,promoters(txdb,upstream = 2000,downstream = 2000))>0
  k562.h1.k4me3.manorm$nonPromoter<-countOverlaps(k562.h1.k4me3.manorm,CTSS.prom)==0
  return(k562.h1.k4me3.manorm)
}
k562.h1.k4me3.manorm<- H1.H3K4me3.state(k562.h1.k4me3.manorm)


k562.h1.k4me3.manorm$Prom<-"N"
k562.h1.k4me3.manorm$Prom[countOverlaps(k562.h1.k4me3.manorm,CTSS.prom)>0]<-"V"
k562.h1.k4me3.manorm$Prom[countOverlaps(k562.h1.k4me3.manorm,promoters(txdb,upstream = 2000,downstream = 2000))>0]<-"P"

k562.h1.k4me3.manorm<- k562.h1.k4me3.manorm[k562.h1.k4me3.manorm$Prom!="V"]
grid.newpage()
T<-venn.diagram(list(Promoter=which(k562.h1.k4me3.manorm$Promoter),Non_Promoter=which(k562.h1.k4me3.manorm$nonPromoter),K562=which(k562.h1.k4me3.manorm$K562),
                     H1hesc=which(k562.h1.k4me3.manorm$H1hesc)),fill=c('darkorange', 'dodgerblue', 'hotpink', 'limegreen'), alpha=c(0.5,0.5,0.5,0.5), cex=2, filename=NULL)
grid.draw(T)

k562.h1.k4me3.manorm$Peaktype<- setMtofactors(k562.h1.k4me3.manorm$M_value)
k562.h1.k4me3.manorm$State<- with(k562.h1.k4me3.manorm,interaction(Prom,Peaktype))
grl.k4me3.ma<- split(k562.h1.k4me3.manorm,as.factor(k562.h1.k4me3.manorm$State))


#read DNase narrow peaks------
require(genomation)
readNarrowPeak <- function(file){
  readGeneric(file,
              strand=6,
              meta.cols=list(name=4,
                             score=5,
                             signalValue=7,
                             pvalue=8,
                             qvalue=9,
                             peak=10),
              header=FALSE,
              zero.based=TRUE)
}

dnase.K562 <- readNarrowPeak("/mnt/local-disk1/rsgeno2/MAmotif/RACK7/routput/wgEncodeAwgDnaseUwdukeK562UniPk.narrowPeak/wgEncodeAwgDnaseUwdukeK562UniPk.narrowPeak")

dnase.K562.NB <- subsetByOverlaps(dnase.K562,grl.k4me3.ma$N.B)
dnase.K562.NB<- resize(dnase.K562.NB,width = 500,fix = "center")
names(dnase.K562.NB)<- 1:length(dnase.K562.NB)
## Dscore of non-promoter H3K4me3

names(gro.seq)<- substring(names(gro.seq),first = 1,last = nchar(names(gro.seq))-9)
gro.cage.seq<- c(gro.seq,cage.seq)

region.Dscore<- function(gr){
  gc.ds<- lapply(seq(from=1,to=length(gro.cage.seq),by=2),function(i,y)Dscore(y[[i+1]],y[[i]],gr,wt = "V5"),y=gro.cage.seq)
  names(gc.ds)<- substring(names(gro.cage.seq),first = 1,last = nchar(names(gro.cage.seq))-c(6,5))[seq(from=1,to=length(gro.cage.seq),by=2)]
  gc.ds.idx<- Reduce(intersect,lapply(gc.ds,names))
  gene.ds<- sapply(gc.ds,function(x)x[gc.ds.idx]$Dscore)
  gene.ds<- as.data.frame(gene.ds)
  gene.ds$strand<- as.character(strand(gr[gc.ds.idx]))
  gene.ds.gp<-melt(gene.ds)
  return(gene.ds.gp)
}
dnase.K562.NB.gp<- region.Dscore(dnase.K562.NB)


ggplot(dnase.K562.NB.gp)+geom_histogram(aes(x=value,y=..count..))+
  facet_wrap( ~ variable)+
  labs(x = "",y = " ",title="Histogram of D score of non-promoter peaks in K562") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size=16)
  )


###

region.signal<- function(gr){
  cage.sg<-lapply(gro.cage.seq,function(x)ScoreMatrixBin(x,gr,bin.num = 50,weight.col = "V5"))
  cage.sg<- melt(sapply(cage.sg,colMeans))
  cage.sg$gene<- rep("*",50)
  cage.sg$read<- c(rep("minus",50),rep("plus",50))
  cage.sg$type <- substring(cage.sg$Var2,first = 1,last = nchar(as.character(cage.sg$Var2))-c(rep(6,50),rep(5,50)))
  cage.sg$Var1<- c(-24.5:24.5)*10
  return(cage.sg)
}

dnase.K562.NB.sg<- region.signal(dnase.K562.NB)




ggplot(dnase.K562.NB.sg)+geom_line(aes(x=Var1,y=value,linetype=as.factor(read)))+
  facet_wrap( ~ type)+
  labs(x = "",y = " ",title="The average signal of non-promoter peaks in K562",linetype="read") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size = 16)
  )


###
all.gr.m<- readRDS("/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/data/hg19.ctss_all_minus_pool.rds")
all.gr.p<- readRDS("/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/data/hg19.ctss_all_plus_pool.rds")
dnase.K562.NB.all.gs<-lapply(list(all.gr.m,all.gr.p),function(x)ScoreMatrixBin(x,dnase.K562.NB,bin.num = 50,weight.col = "V5"))

all.sg<- melt(sapply(dnase.K562.NB.all.gs,colMeans))
all.sg$gene<- rep("*",50)
all.sg$read<- c(rep("minus",50),rep("plus",50))
all.sg$Var1<- c(-24.5:24.5)*10

ggplot(all.sg)+geom_line(data=all.sg,aes(x=Var1,y=value,linetype=as.factor(read)))+
  labs(x = "",y = " ",title="The average signal of non-promoter peaks",linetype="read") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size = 16)
  )

dnase.K562.NB.all.ds<- Dscore(all.gr.p,all.gr.m,dnase.K562.NB,wt="V5")
ggplot(as.data.frame(dnase.K562.NB.all.ds))+geom_freqpoly(aes(x=Dscore,y=..count..))+
  labs(x = "",y = " ",title="Histogram of D score of non-promoter peaks by merged CAGE") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size=16)
  )

dnase.K562.NB.all.ds<- Dscore(all.gr.p,all.gr.m,ucsc.gene.pt,wt="V5")
ggplot(as.data.frame(dnase.K562.NB.all.ds))+geom_histogram(aes(x=Dscore,y=..count..))+
  labs(x = "",y = " ",title="Histogram of D score of promoters by merged CAGE") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size=16)
  )