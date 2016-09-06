
#classification and annotation dnase
##load packages------------------
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(GenomicRanges)
require(genomation)
require(ChIPseeker)
#preprocess H3K4me3 after running MAnorm by K562 and H1======
k562.h1.manorm<- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me3/K562.H1.rep1.top2k.20160520/K562.H1.rep1.top2k.20160520_all_peak_MAvalues.xls",header = TRUE)
k562.h1.manorm$state<- "None"
k562.h1.manorm$state[k562.h1.manorm$M_value>=1]<- "K562"
k562.h1.manorm$state[k562.h1.manorm$M_value<1 & k562.h1.manorm$M_value> -1]<- "Common"
k562.h1.manorm$state[k562.h1.manorm$M_value <= -1]<- "H1hesc"
##confident
k562.both.ma<- list(H3K4me3 = readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me3/K562.H3K4me3.both.top20k.20160824/K562.H3K4me3.both.top20k.20160824_all_peak_MAvalues.xls"),
                    H3K4me1 = readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me1/K562.H3K4me1.both.top20k.20160602/K562.H3K4me1.both.top20k.20160602_all_peak_MAvalues.xls"),
                    H3K27ac = readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K27ac/K562.H3K27ac.both.top20k.20160607/K562.H3K27ac.both.top20k.20160607_all_peak_MAvalues.xls"))
trimMA<- function(gr){
  gr<- gr[gr$M_value<1 & gr$M_value>-1]
  gr$score<- gr$summit
  return(gr)
}
k562.both.ma<- lapply(k562.both.ma,trimMA)

dnase.K562 <- readNarrowPeak("/mnt/local-disk1/rsgeno2/MAmotif/RACK7/routput/wgEncodeAwgDnaseUwdukeK562UniPk.narrowPeak/wgEncodeAwgDnaseUwdukeK562UniPk.narrowPeak")
##annotation
dnase.K562.anno <- annotatePeak(dnase.K562, tssRegion=c(-2000, 2000), TxDb =txdb, annoDb="org.Hs.eg.db")
plotAnnoPie(dnase.K562.anno)
dnase.K562.anno<- as.GRanges(dnase.K562.anno)
dnase.K562.anno$H3K4me3<- "None"
dnase.K562.anno$H3K4me3[countOverlaps(dnase.K562.anno,k562.both.ma$H3K4me3)>0]<- "H3K4me3"
dnase.K562.anno$H3K27ac<- "None"
dnase.K562.anno$H3K27ac[countOverlaps(dnase.K562.anno,k562.both.ma$H3K27ac)>0]<- "H3K27ac"
dnase.K562.anno$H3K4me1<- "None"
dnase.K562.anno$H3K4me1[countOverlaps(dnase.K562.anno,k562.both.ma$H3K4me1)>0]<- "H3K4me1"

dnase.K562.anno$Promoter<- FALSE
dnase.K562.anno$Promoter[dnase.K562.anno$annotation=="Promoter (<=1kb)"]<- TRUE
dnase.K562.anno$Promoter[dnase.K562.anno$annotation=="Promoter (1-2kb)"]<- TRUE
dnase.K562.anno$annotation[grep("Intron",dnase.K562.anno$annotation)]<- "Intron"
dnase.K562.anno$annotation[grep("Exon",dnase.K562.anno$annotation)]<- "Exon"
dnase.K562.anno$annotation[grep("Downstream",dnase.K562.anno$annotation)]<- "Downstream"
dnase.K562.anno$State<- "None"
ind <- as.matrix(findOverlaps(dnase.K562.anno,k562.h1.manorm))
dnase.K562.anno$State[ind[,1]]<- k562.h1.manorm$state[ind[,2]]
dnase.K562.anno$Group<-with(dnase.K562.anno,paste(annotation,H3K4me3,H3K27ac,H3K4me1,State,sep = "."))
saveRDS(dnase.K562.anno,file = "data/Dnase_K562_anno.rds")
require(VennDiagram)
grid.newpage()
T<-venn.diagram(with(dnase.K562.anno,list(H3K4me3=which(H3K4me3),Promoter=which(Promoter),H3K27ac=which(H3K27ac),H3K4me1=which(H3K4me1))),fill=c('darkorange', 'dodgerblue', 'hotpink', 'limegreen'), alpha=c(0.5,0.5,0.5,0.5), cex=2, filename=NULL)
grid.draw(T)

dnase.K562.anno.s<- resize(dnase.K562.anno,width = 500,fix = "center")
dnase.K562.anno.sl<- split(dnase.K562.anno.s,dnase.K562.anno.s$Group)
dnase.K562.anno.sl<- dnase.K562.anno.sl[sapply(dnase.K562.anno.sl,function(x)length(x)>500)]



#####
sml = ScoreMatrixBin(cage.seq[[2]],windows = dnase.K562.anno.sl[[7]],weight.col = "V5",bin.num = 50)
sml.scaled <- scaleScoreMatrix(sml)
multiHeatMatrix(sml,
                clustfun=function(x) kmeans(x, centers=3)$cluster,
                cex.axis=0.8,xcoords=c(-25,25),
                winsorize=c(0,99),
                legend.name=names(sml),xlab="region around TSS")

plotMeta(sml)

require(rtracklayer)
tf.bw.files <- list.files("/mnt/local-disk1/rsgeno2/MAmotif/ENCODE/Tfbs/SydhTfbsK562/",pattern = "bigWig")[1:20]

lapply(1:length(tf.bw.files),function(i,tf){
  Mazab<- import.bw(paste0("/mnt/local-disk1/rsgeno2/MAmotif/ENCODE/Tfbs/SydhTfbsK562/",tf[i]),asRle = TRUE)
  sml = ScoreMatrixBin(Mazab,windows = dnase.K562.anno.sl[[7]],bin.num = 50)
  sml.scaled <- scaleScoreMatrix(sml)
  
  heatMatrix(sml,
             clustfun=function(x) kmeans(x, centers=3)$cluster,
             cex.axis=0.8,xlab=tf[i])
  plotMeta(sml,xlab = tf[i])
  heatMatrix(sml.scaled,
             clustfun=function(x) kmeans(x, centers=3)$cluster,
             cex.axis=0.8,xlab=tf[i])
  plotMeta(sml.scaled,xlab = tf[i])
  
},tf=tf.bw.files)

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



multiHeatMatrix(sml,
                clustfun=function(x) kmeans(x, centers=3)$cluster,
                cex.axis=0.8,xcoords=c(-25,25),
                winsorize=c(0,99),
                legend.name=names(sml),xlab="region around TSS")

plotMeta(sml)