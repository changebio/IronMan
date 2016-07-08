### =========================================================================
### calculate the density of a bed file in different groups of peaks 
### -------------------------------------------------------------------------

#read DNase narrow peaks------
dnase.K562 <- readNarrowPeak("/mnt/local-disk1/rsgeno2/MAmotif/RACK7/routput/wgEncodeAwgDnaseUwdukeK562UniPk.narrowPeak/wgEncodeAwgDnaseUwdukeK562UniPk.narrowPeak")

gro1<- import.bed("/mnt/local-disk1/rsgeno2/MAmotif/For_huang_K562_CAGE_GROseq/K562_biol_rep123.hg19.ctss_minus_pool.bed")
k4me3_NB<- lapply(grl.k4me3.ma,function(x)GRanges(seqnames = seqnames(x),ranges = IRanges(start = start(x)+x$summit-1000,width = 2000)))

#get score matrix in peak regions, which must be the same length
k4me3_NB.dnase<- lapply(k4me3_NB,function(x)ScoreMatrix(dnase.K562,x))
k4me3_NB.dnase1<- sapply(k4me3_NB.dnase,colMeans)
ggplot(melt(k4me3_NB.dnase1[,c(1,3,5)]))+geom_line(aes(x=Var1,y=value,colour=Var2))+
  labs(x = "",y = " ",title="The average DNase signal in different enhancer",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )

k4me3_NB.gro<- lapply(k4me3_NB,function(x)ScoreMatrixBin(gro1,x,bin.num = 50))
k4me3_NB.gro1<- sapply(k4me3_NB.gro,colMeans)
ggplot(melt(k4me3_NB.gro1[,c(1,3,5)]))+geom_line(aes(x=Var1,y=value,colour=Var2))+
  labs(x = "",y = " ",title="The average GRO-Seq signal in different enhancer",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )

k4me3_NB.dnase.bin<- ScoreMatrixBin(dnase.K562,k4me3_NB,bin.num = 50)


# ---------------------------------------------------------------------------- #
#' Get common rows from all matrices in a ScoreMatrixList object
#' 
#' Returns a intersection of rows for each matrix in a ScoreMatrixList object. 
#' This is done using the rownames of each element in the list.
#' 
k4me3_NB.dnase.lst<- ScoreMatrixList(dnase.K562,k4me3_NB)


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
