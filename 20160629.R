### =========================================================================
### H3K4me1 and H3K27ac overlap with K562-unique non-promoter peaks
### -------------------------------------------------------------------------

##load packages------------------
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(GenomicRanges)
#preprocess H3K4me3 after running MAnorm by K562 and H1======
k562.h1.k4me3.manorm<- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me3/K562.H1.rep1.top2k.20160520/K562.H1.rep1.top2k.20160520_all_peak_MAvalues.xls",header = TRUE)
k562.h1.k4me3.manorm<- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me3/K562.H1.rep2.top2k.20160524/K562.H1.rep2.top2k.20160524_all_peak_MAvalues.xls",header = TRUE)

k562.gm12878.manorm <- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me3/K562.Gm12878.H3K4me3.rep1.top20k.20160808/K562.Gm12878.H3K4me3.rep1.top20k.20160808_all_peak_MAvalues.xls",header=TRUE)
k562.gm12878.manorm <- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me3/K562.Gm12878.H3K4me3.rep2.top20k.20160808/K562.Gm12878.H3K4me3.rep2.top20k.20160808_all_peak_MAvalues.xls",header=TRUE)
temp <- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/Roadmap/H3K4me3T20k/Mesendoderm.H1.top20k.20160811/Mesendoderm.H1.top20k.20160811_all_peak_MAvalues.xls",header=TRUE)
k562.h1.k4me3.manorm<- temp

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

require(VennDiagram)
require(gridExtra)
h1.k4me3.file<-list.files("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me3",pattern = "H1.H3K4me3.*?.xls",full.names = TRUE,recursive = TRUE)
h1.k4me3.pks<- lapply(h1.k4me3.file,function(x)readPeakFile(x,header=TRUE))
names(h1.k4me3.pks)<- list.files("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me3",pattern = "H1.H3K4me3.",full.names = FALSE)

gs<- lapply(seq_along(h1.k4me3.pks), function(x, n, i){k562.h1.k4me3.manorm<- H1.H3K4me3.state(x[[i]]);
                                                  temp<- list(Promoter=which(k562.h1.k4me3.manorm$Promoter),Non_Promoter=which(k562.h1.k4me3.manorm$nonPromoter),H1hesc=which(k562.h1.k4me3.manorm$H1hesc));
                                                  temp[[unlist(strsplit(n[[i]],"[.]"))[1]]]<- which(k562.h1.k4me3.manorm$K562);
                                                  grobTree(venn.diagram(temp,fill=c('darkorange', 'dodgerblue', 'hotpink', 'limegreen'), alpha=c(0.5,0.5,0.5,0.5), filename=NULL)
                                                  )} , x=h1.k4me3.pks, n=names(h1.k4me3.pks))
grid.arrange(grobs=gs,ncol=5)
file.remove(list.files(pattern = "VennDiagram*"))


# grid.newpage()
# T<-venn.diagram(list(Promoter=which(k562.h1.k4me3.manorm$Promoter),Non_Promoter=which(k562.h1.k4me3.manorm$nonPromoter),K562=which(k562.h1.k4me3.manorm$K562),
#                      H1hesc=which(k562.h1.k4me3.manorm$H1hesc)),fill=c('darkorange', 'dodgerblue', 'hotpink', 'limegreen'), alpha=c(0.5,0.5,0.5,0.5), cex=2, filename=NULL)
# grid.draw(T)


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

#H3K4me3 peak annotations
require(ChIPseeker)
k562.h1.k4me3.manorm.anno <- annotatePeak(k562.h1.k4me3.manorm, tssRegion=c(-2000, 2000), TxDb =txdb, annoDb="org.Hs.eg.db")
saveRDS(k562.h1.k4me3.manorm.anno,file = "data/annotated_manorm_H3K4me3.rds")

k562.h1.k4me3.manorm.anno<- readRDS("data/annotated_manorm_H3K4me3.rds")
k562.h1.k4me3.manorm<- k562.h1.k4me3.manorm.anno@anno
k562.h1.k4me3.manorm$Feature<- sub(" .*","",k562.h1.k4me3.manorm$annotation)
#------------
table(k562.h1.k4me3.manorm$Feature,k562.h1.k4me3.manorm$Promoter)
tp <- sort(k562.h1.k4me3.manorm)
tp$score=as.numeric(tp$Promoter)
export.ucsc(tp,con="data/temp_H3K4me3.bed")

#------------
boxplot(M_value~State,data=mcols(k562.h1.k4me3.manorm))
abline(h=1)
abline(h=-1)

#read H3K4me1 macs peaks======
k4me1.mc<- list(K562=readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/macs.20151224/H3K4me1.top20k/K562.H3k4me1.Rep2.macs.20160512_Top20K_peaks.xls"),
             H1hesc=readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/macs.20151224/H3K4me1.top20k/H1.H3k4me1.Rep2.macs.20160512_Top20K_peaks.xls"))
# k4me1.mc<- list(K562=readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/macs.20151224/H3K27ac.top20k/K562.H3k27ac.Rep2.macs.20160512_Top20K_peaks.xls"),
#                 H1hesc=readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/macs.20151224/H3K27ac.top20k/H1.H3k27ac.rep2.macs.20151224_Top20K_peaks.xls"))

#k4me1.mc<- lapply(list.files("/mnt/local-disk1/rsgeno2/MAmotif/macs.20151224/H3K27ac.top20k/",pattern = ".xls"),function(x)readPeakFile(paste0("/mnt/local-disk1/rsgeno2/MAmotif/macs.20151224/H3K27ac.top20k/",x)))

k4me1.d<- as.data.frame(sapply(k4me1.mc,function(x)countOverlaps(grl.k4me3.ma,x)))
k4me1.d$len<- sapply(grl.k4me3.ma,length)
k4me1.d$K562p <- k4me1.d$K562/k4me1.d$len
k4me1.d$H1p<- k4me1.d$H1hesc/k4me1.d$len
k4me1.d$type<- row.names(k4me1.d)
#barplot----
require(ggplot2)
require(reshape2)
ggplot(melt(k4me1.d[,4:6]),aes(x=type,y=value,fill=variable))  + geom_bar(stat = "identity",position ="dodge" ) +
  labs(x = "",y = "",title="The persentage of H3K4me1 in different kinds of H3K4me3",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        #axis.title.y = element_text(color="#993333", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = "NA", colour = "black", size = 2)
  )

#read H3K4me1 MAnorm peaks
k562.h1.K4me1.manorm<- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me1/K562.H1.H3K4me1.rep2.top20k.20160525/K562.H1.H3K4me1.rep2.top20k.20160525_all_peak_MAvalues.xls",header = TRUE)
k562.h1.K4me1.manorm<- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K27ac/K562.H1.H3K27ac.rep2.top20k.20160525/K562.H1.H3K27ac.rep2.top20k.20160525_all_peak_MAvalues.xls",header = TRUE)

k562.h1.K4me1.manorm$H3K4me3<- setMAvalbylap(k562.h1.K4me1.manorm,k562.h1.k4me3.manorm[k562.h1.k4me3.manorm$Prom=="N"],type="State")
k562.h1.K4me1.manorm$H3K4me3<- setMAvalbylap(k562.h1.K4me1.manorm,k562.h1.k4me3.manorm,type="State")
k562.h1.K4me1.manorm$H3K4me3.mval<- setMAvalbylap(k562.h1.K4me1.manorm,k562.h1.k4me3.manorm)
with(na.omit(mcols(k562.h1.K4me1.manorm)),boxplot(M_value~H3K4me3))
#boxplot----
ggplot(as.data.frame(na.omit(mcols(k562.h1.K4me1.manorm))),aes(H3K4me3,M_value))  + geom_boxplot(colour = "black",alpha=0.5,size = 2,width = 0.5) + geom_jitter(width = 0.25,colour = "#3366FF")+
  labs(x = "",y = "",title="The M value of H3K27ac in different kinds of H3K4me3") +
  theme(plot.title = element_text(color="black", size=20, face="bold"),
        axis.title.x = element_text( face="bold",size=14),
        #axis.title.y = element_text(color="#993333", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = "NA", colour = "black", size = 2)
  )

ggplot(as.data.frame(na.omit(mcols(k562.h1.K4me1.manorm))),aes(x=H3K4me3.mval,y=M_value))  + geom_point(aes(colour = factor(H3K4me3))) +
  labs(x = "",y = "",title="The relationship of M value between H3K27ac and H3K4me3",colour="") +
  theme(plot.title = element_text(color="black", size=20, face="bold"),
        axis.title.x = element_text( face="bold",size=14),
        #axis.title.y = element_text(color="#993333", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = "NA", colour = "black", size = 2)
  )