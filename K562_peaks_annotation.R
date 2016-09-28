require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(ChIPseeker)
require(clusterProfiler)
require(genomation)

##read manorm results
ma.h3k4me1<-readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me1/K562.H3K4me1.both.top30k.20160607/K562.H3K4me1.both.top30k.20160607_all_peak_MAvalues.xls")
ma.h3k4me1<- ma.h3k4me1[ma.h3k4me1$M_value>= -1 & ma.h3k4me1$M_value<= 1]
ma.h3k4me3<- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me3/K562.H3K4me3.both.top20k.20160824/K562.H3K4me3.both.top20k.20160824_all_peak_MAvalues.xls")
ma.h3k4me3<- ma.h3k4me3[ma.h3k4me3$M_value>= -1 & ma.h3k4me3$M_value<= 1]
ma.h3k27ac<- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K27ac/K562.H3K27ac.both.top20k.20160607/K562.H3K27ac.both.top20k.20160607_all_peak_MAvalues.xls")
ma.h3k27ac<- ma.h3k27ac[ma.h3k27ac$M_value>= -1 & ma.h3k27ac$M_value<= 1]
h3k27me3<- readPeakFile("/mnt/local-disk1/rsgeno2/huangyin/PRC2/K562_H3K27me3/wgEncodeSydhHistoneK562bH3k27me3bUcdPk.narrowPeak.gz")[1:30000]

dnase.K562 <- readNarrowPeak("/mnt/local-disk1/rsgeno2/MAmotif/RACK7/routput/wgEncodeAwgDnaseUwdukeK562UniPk.narrowPeak/wgEncodeAwgDnaseUwdukeK562UniPk.narrowPeak")


##annotation H3K4me3
ma.h3k4me3.anno <- annotatePeak(ma.h3k4me3, tssRegion=c(-2000, 2000), TxDb =txdb, annoDb="org.Hs.eg.db")
ma.h3k4me3.anno<- as.GRanges(ma.h3k4me3.anno)[,8:19]
ma.h3k4me3.anno$H3K4me1<- FALSE
ma.h3k4me3.anno$H3K4me1[countOverlaps(ma.h3k4me3.anno,ma.h3k4me1)>0]<- TRUE
ma.h3k4me3.anno$H3K27ac<- FALSE
ma.h3k4me3.anno$H3K27ac[countOverlaps(ma.h3k4me3.anno,ma.h3k27ac)>0]<- TRUE
CTSS.prom <- readRDS("data/CTSS_promoters.rds")
ma.h3k4me3.anno$Prom<-"Non"
ma.h3k4me3.anno$Prom[countOverlaps(ma.h3k4me3.anno,CTSS.prom)>0]<-"V"
ma.h3k4me3.anno$Prom[countOverlaps(ma.h3k4me3.anno,promoters(txdb,upstream = 2000,downstream = 2000))>0]<-"Pro"
ma.h3k4me3.anno<- ma.h3k4me3.anno[ma.h3k4me3.anno$Prom!="V"]
ma.h3k4me3.anno$Promoter<- FALSE
ma.h3k4me3.anno$Promoter[ma.h3k4me3.anno$annotation=="Promoter (<=1kb)"]<- TRUE
ma.h3k4me3.anno$Promoter[ma.h3k4me3.anno$annotation=="Promoter (1-2kb)"]<- TRUE
ma.h3k4me3.anno$State<- "None"
ma.h3k4me3.anno$State[ma.h3k4me3.anno$H3K4me1]<- "H3K4me1"
ma.h3k4me3.anno$State[ma.h3k4me3.anno$H3K27ac]<- "H3K27ac"
ma.h3k4me3.anno$State[ma.h3k4me3.anno$H3K4me1 & ma.h3k4me3.anno$H3K27ac]<- "Both"

ma.h3k4me3.anno<- ma.h3k4me3.anno[countOverlaps(ma.h3k4me3.anno,h3k27me3)==0]
ma.h3k4me3.anno<- ma.h3k4me3.anno[!(ma.h3k4me3.anno$Prom=="Pro" & !ma.h3k4me3.anno$Promoter)]

dnase.me3.hits<- findOverlaps(ma.h3k4me3.anno,dnase.K562)
dnase.me3<- dnase.K562[dnase.me3.hits@subjectHits]
dnase.me3$qurey<- dnase.me3.hits@queryHits
dnase.me3<- sort(dnase.me3,by = ~ qurey + signalValue,decreasing = TRUE)
dnase.me3<- dnase.me3[!duplicated(dnase.me3$qurey)]
mcols(dnase.me3)<- mcols(ma.h3k4me3.anno[dnase.me3$qurey])
dnase.me3$State[dnase.me3$Promoter]<- "Promoter"

###annotate non-promoter enhancer
enh.k562<- reduce(c(ma.h3k4me1[,1],ma.h3k27ac[,1]))
enh.k562<- enh.k562[countOverlaps(enh.k562,ma.h3k4me3)==0]
enh.k562<- enh.k562[countOverlaps(enh.k562,h3k27me3)==0]
enh.k562.anno <- annotatePeak(enh.k562, tssRegion=c(-2000, 2000), TxDb =txdb, annoDb="org.Hs.eg.db")
enh.k562.anno<- as.GRanges(enh.k562.anno)
enh.k562.anno$Promoter<- FALSE
enh.k562.anno$Promoter[enh.k562.anno$annotation=="Promoter (<=1kb)"]<- TRUE
enh.k562.anno$Promoter[enh.k562.anno$annotation=="Promoter (1-2kb)"]<- TRUE
enh.k562.anno$Prom<-"Non"
enh.k562.anno$Prom[countOverlaps(enh.k562.anno,CTSS.prom)>0]<-"V"
enh.k562.anno$Prom[countOverlaps(enh.k562.anno,promoters(txdb,upstream = 2000,downstream = 2000))>0]<-"Pro"
enh.k562.anno<- enh.k562.anno[enh.k562.anno$Prom!="V"]
enh.k562.anno<- enh.k562.anno[enh.k562.anno$Prom=="Non" & !enh.k562.anno$Promoter]

dnase.enh.hits<- findOverlaps(enh.k562.anno,dnase.K562)
dnase.enh<- dnase.K562[dnase.enh.hits@subjectHits]
dnase.enh$qurey<- dnase.enh.hits@queryHits
dnase.enh<- sort(dnase.enh,by = ~ qurey + signalValue,decreasing = TRUE)
dnase.enh<- dnase.enh[!duplicated(dnase.enh$qurey)]
mcols(dnase.enh)<- mcols(enh.k562.anno[dnase.enh$qurey])
dnase.enh$State<- "Enhancer"

k562.pk.dnase<- c(dnase.me3[,c(1:12,17)],dnase.enh[,c(1:12,15)])
k562.pk.dnase<- sort(k562.pk.dnase)
k562.pk.dnase<- sort(k562.pk.dnase,by = ~ State)
k562.pk.dnase.5h<- resize(k562.pk.dnase,width = 500,fix="center")
k562.pk.dnase.5h<- k562.pk.dnase.5h[countOverlaps(k562.pk.dnase.5h,k562.pk.dnase.5h)==1]
k562.5h.s<- split(k562.pk.dnase.5h,k562.pk.dnase.5h$State)





