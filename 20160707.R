### =========================================================================
### peak character(active or non-active) 
### calculate the density of a bed file in different groups of peaks 
### -------------------------------------------------------------------------
k4me3.smt<- lapply(grl.k4me3.ma,function(x)GRanges(seqnames = seqnames(x),ranges = IRanges(start = start(x)+x$summit-1000,width = 2000)))
k4me3.smt<- lapply(k4me3.smt,function(x)dropSeqlevels(x,c("chrY","chrM")))
k4me3.smt<- GRangesList(k4me3.smt)
k4me3_gr<- stack(k4me3.smt)
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


#gro1<- import.bed("/mnt/local-disk1/rsgeno2/MAmotif/For_huang_K562_CAGE_GROseq/K562_biol_rep123.hg19.ctss_minus_pool.bed")

#get score matrix in peak regions, which must be the same length
k4me3.smt.dnase<- lapply(k4me3.smt,function(x)ScoreMatrix(dnase.K562,x))
k4me3.smt.dnase1<- sapply(k4me3.smt.dnase,colMeans)
ggplot(melt(k4me3.smt.dnase1[,c(1,3,5)]))+geom_line(aes(x=Var1,y=value,colour=Var2))+
  labs(x = "",y = " ",title="The average DNase signal in different enhancer",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )

gro.seq<- lapply(list.files("/mnt/local-disk1/rsgeno2/MAmotif/For_huang_K562_CAGE_GROseq/",pattern = ".bed$")[1:6],function(x)readPeakFile(paste0("/mnt/local-disk1/rsgeno2/MAmotif/For_huang_K562_CAGE_GROseq/",x)))
names(gro.seq)<-list.files("/mnt/local-disk1/rsgeno2/MAmotif/For_huang_K562_CAGE_GROseq/",pattern = ".bed$")[1:6]
gro.seq<- lapply(gro.seq,function(x)dropSeqlevels(x,c("chrY","chrM")))

# k4me3.gro<- lapply(gro.seq,function(gro){
#   lapply(k4me3.smt,function(x)ScoreMatrix(gro,x))
# })
k4me3.gro<-lapply(k4me3.smt,function(x)ScoreMatrixList(gro.seq,x,weight.col = "V5"))
k4me3.gro<- lapply(k4me3.gro,intersectScoreMatrixList)
k4me3.gro1<- lapply(k4me3.gro,function(gro){
  lapply(gro,rowSums)
})

a<-lapply(k4me3.gro1,function(x){(x[[2]]-x[[1]])/(x[[2]]+x[[1]])})
hist(a$P.B,xlab = "D score",main = "Histogram of D score in K562-unique promoter H3K4me3")
k4me3.smt.gro<- lapply(k4me3.smt,function(x)ScoreMatrixBin(gro.seq[[1]],x,bin.num = 50,weight.col = "V5"))
k4me3.smt.gro1<- sapply(k4me3.smt.gro,colMeans)
tp<- k4me3.smt.gro1
ggplot(melt((tp+k4me3.smt.gro1)/2))+geom_line(aes(x=Var1,y=value,colour=Var2))+
  labs(x = "",y = " ",title="The average GRO-Seq signal in different enhancer",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )

k4me3.smt.dnase.bin<- ScoreMatrixBin(dnase.K562,k4me3.smt,bin.num = 50)

#negtive control
require(rtracklayer)
k27me3.K562<- import.bed("/mnt/local-disk1/rsgeno2/huangyin/PRC2/K562_H3K27me3/GSM733658_hg19_wgEncodeBroadHistoneK562H3k27me3StdAlnRep1.bed")
k4me3.smt.k27me3<- lapply(k4me3.smt,function(x)ScoreMatrix(k27me3.K562,x))
k4me3.smt.k27me31<- sapply(k4me3.smt.k27me3,colMeans)
ggplot(melt(k4me3.smt.k27me31[,c(1,3,5)]))+geom_line(aes(x=Var1,y=value,colour=Var2))+
  labs(x = "",y = " ",title="The average H3K27me3 signal in different enhancer",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )

suz12.K562<- import.bed("/mnt/local-disk1/rsgeno2/huangyin/PRC2/ChIP_raw_data/K562_FB_EED_GCCAAT_L007_R1.bed")
k4me3.smt.suz12<- lapply(k4me3.smt,function(x)ScoreMatrix(a1,x))
k4me3.smt.suz121<- sapply(k4me3.smt.suz12,colMeans)
ggplot(melt(k4me3.smt.suz121))+geom_line(aes(x=Var1,y=value,colour=Var2))+
  labs(x = "",y = " ",title="The average SUZ12 signal in different enhancer",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )


# ---------------------------------------------------------------------------- #
#' Get common rows from all matrices in a ScoreMatrixList object
#' 
#' Returns a intersection of rows for each matrix in a ScoreMatrixList object. 
#' This is done using the rownames of each element in the list.
#' 
k4me3.smt.dnase.lst<- ScoreMatrixList(dnase.K562,k4me3.smt)




# ---------------------------------------------------------------------------- #
#' Heatmap for meta-region profiles
#' 
#' Function calculates meta-profile(s) from a ScoreMatrix or a ScoreMatrixList, then
#' produces a heatmap or a set of stacked heatmaps for meta-region profiles
#' 


#######DNase signal in enhancer
require(genomation)

ma.h3k4me1.enh<- resize(ma.h3k4me1.anno.s[[1]],width = 2000,fix = "center")
ma.h3k4me1.enh.lst<- split(ma.h3k4me1.enh,ma.h3k4me1.enh$State)
ma.h3k4me1.pro<- resize(ma.h3k4me1.anno.s[[2]],width = 2000,fix = "center")
ma.h3k4me1.pro.lst<- split(ma.h3k4me1.pro,ma.h3k4me1.pro$State)
###TF and DNase files
require(rtracklayer)
bed.files<-list(DNase.rep1="/mnt/local-disk1/rsgeno2/MAmotif/ENCODE/2.DNase_Duke_hg19/wgEncodeOpenChromDnaseK562AlnRep1.bed",
                DNase.rep2="/mnt/local-disk1/rsgeno2/MAmotif/ENCODE/2.DNase_Duke_hg19/wgEncodeOpenChromDnaseK562AlnRep2.bed",
                H3K4me3.rep1="/mnt/local-disk1/rsgeno2/MAmotif/3.Histone_Broad_hg19/H3K4me3/K562/wgEncodeBroadHistoneK562H3k4me3StdAlnRep1.bed",
                H3K4me3.rep2="/mnt/local-disk1/rsgeno2/MAmotif/3.Histone_Broad_hg19/H3K4me3/K562/wgEncodeBroadHistoneK562H3k4me3StdAlnRep2.bed",
                H3K4me1.rep1="/mnt/local-disk1/rsgeno2/MAmotif/3.Histone_Broad_hg19/H3K4me1/wgEncodeBroadHistoneK562H3k4me1StdAlnRep1.bed",
                H3K4me1.rep2="/mnt/local-disk1/rsgeno2/MAmotif/3.Histone_Broad_hg19/H3K4me1/wgEncodeBroadHistoneK562H3k4me1StdAlnRep2.bed")

DNase.bed<- lapply(bed.files,import.bed)


#bam.files<- c("/mnt/local-disk1/rsgeno2/MAmotif/RACK7/wgEncodeOpenChromChipK562Pol2AlnRep1.bam","/mnt/local-disk1/rsgeno2/MAmotif/RACK7/wgEncodeOpenChromChipK562Pol2AlnRep2.bam")
ma.h3k4me1.enh.sml.lst<- lapply(ma.h3k4me1.enh.lst,ScoreMatrixList,targets=DNase.bed, bin.num=100, cores=2)
ma.h3k4me1.pro.sml.lst<- lapply(ma.h3k4me1.pro.lst,ScoreMatrixList,targets=DNase.bed, bin.num=100, cores=2)
dnase.rep1<- sapply(ma.h3k4me1.enh.sml.lst,function(x)colMeans(x[[1]]))
#dnase.rep1<- sapply(ma.h3k4me1.pro.sml.lst,function(x)colMeans(x[[1]]))
ggplot(melt(dnase.rep1))+geom_line(aes(x=Var1,y=value,colour=Var2))+
  labs(x = "",y = " ",title="The average DNase signal in different enhancer",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )

#dnase.rep2<- sapply(ma.h3k4me1.enh.sml.lst,function(x)colMeans(x[[2]]))
dnase.rep2<- sapply(ma.h3k4me1.pro.sml.lst,function(x)colMeans(x[[2]]))
ggplot(melt(dnase.rep2))+geom_line(aes(x=Var1,y=value,colour=Var2))+
  labs(x = "",y = " ",title="The average DNase signal in different enhancer",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )


ma.h3k4me1.enh.scale<- scaleScoreMatrixList(ma.h3k4me1.enh.sml)
# k-means algorithm with 2 clusters
cl1 <- function(x) kmeans(x, centers=4)$cluster
multiHeatMatrix(ma.h3k4me1.enh.scale, xcoords=c(-1000, 1000), group = factor(ma.h3k4me1.enh$State[-c(1,2)]),clustfun = cl1)


genomationDataPath = system.file('extdata',package='genomationData')
bam.files = list.files(genomationDataPath, full.names=TRUE, pattern='bam$')
bam.files = bam.files[!grepl('Cage', bam.files)]

ctcf.peaks = readBroadPeak(file.path(genomationDataPath, 
                                     'wgEncodeBroadHistoneH1hescCtcfStdPk.broadPeak.gz'))
ctcf.peaks = ctcf.peaks[seqnames(ctcf.peaks) == 'chr21']
ctcf.peaks = ctcf.peaks[order(-ctcf.peaks$signalValue)]
ctcf.peaks = resize(ctcf.peaks, width=1000, fix='center')

sml = ScoreMatrixList(bam.files, ctcf.peaks, bin.num=50, type='bam', cores=2)

# descriptions of file that contain info. about transcription factors
sampleInfo = read.table(system.file('extdata/SamplesInfo.txt',
                                    package='genomationData'),header=TRUE, sep='\t')
names(sml) = sampleInfo$sampleName[match(names(sml),sampleInfo$fileName)]



sml = ScoreMatrixList(bam.files, ctcf.peaks, bin.num=50, type='bam', cores=2)

# descriptions of file that contain info. about transcription factors
sampleInfo = read.table(system.file('extdata/SamplesInfo.txt',
                                    package='genomationData'),header=TRUE, sep='\t')
names(sml) = sampleInfo$sampleName[match(names(sml),sampleInfo$fileName)]
sml.scaled = scaleScoreMatrixList(sml)

# k-means algorithm with 2 clusters
cl1 <- function(x) kmeans(x, centers=2)$cluster
multiHeatMatrix(ma.h3k4me1.enh.sml, xcoords=c(-1000, 1000), clustfun = cl1)


covplot(ma.h3k4me1, weightCol="normalized_read_density_in_wgEncodeBroadHistoneK562H3k4me1StdAlnRep1")
covplot(peak, weightCol="V5", chrs=c("chr17", "chr18"), xlim=c(4.5e7, 5e7))
ma.h3k4me3.anno<- annotatePeak(ma.h3k4me3, tssRegion=c(-2000, 2000), TxDb =txdb, annoDb="org.Hs.eg.db")
promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)
tagMatrix <- getTagMatrix(ma.h3k27ac, windows=promoter)

tagHeatmap(tagMatrix, xlim=c(-2000, 2000), color="red")
plotAvgProf(tagMatrix, xlim=c(-2000, 2000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAvgProf(tagMatrix, xlim=c(-2000, 2000), conf = 0.95, resample = 500)
plotAnnoPie(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
#########
ma.h3k4me1.enh.bp <- lapply(ma.h3k4me1.enh.lst,function(x)enrichGO(x$geneId, ont="BP", readable=TRUE))
head(summary(bp), n=3)
dotplot(ma.h3k4me1.enh.bp$H3K4me3, showCategory=20)
ma.h3k4me1.pro.bp <- lapply(ma.h3k4me1.pro.lst,function(x)enrichGO(x$geneId, ont="BP", readable=TRUE))
dotplot(ma.h3k4me1.pro.bp$H3K4me3, showCategory=20)
