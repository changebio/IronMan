require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(ChIPseeker)
require(clusterProfiler)

##read manorm results
ma.h3k4me1<-readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me1/K562.H3K4me1.both.top30k.20160607/K562.H3K4me1.both.top30k.20160607_all_peak_MAvalues.xls")
ma.h3k4me1<- ma.h3k4me1[ma.h3k4me1$M_value>= -1 & ma.h3k4me1$M_value<= 1]
ma.h3k4me3<- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me3/K562.H3K4me3.both.top20k.20160607/K562.H3K4me3.both.top20k.20160607_all_peak_MAvalues.xls")
ma.h3k4me3<- ma.h3k4me3[ma.h3k4me3$M_value>= -1 & ma.h3k4me3$M_value<= 1]
ma.h3k27ac<- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K27ac/K562.H3K27ac.both.top20k.20160607/K562.H3K27ac.both.top20k.20160607_all_peak_MAvalues.xls")
ma.h3k27ac<- ma.h3k27ac[ma.h3k27ac$M_value>= -1 & ma.h3k27ac$M_value<= 1]


##annotation
ma.h3k4me1.anno <- annotatePeak(ma.h3k4me1, tssRegion=c(-2000, 2000), TxDb =txdb, annoDb="org.Hs.eg.db")
ma.h3k4me1.anno<- as.GRanges(ma.h3k4me1.anno)
ma.h3k4me1.anno$H3K4me3<- FALSE
ma.h3k4me1.anno$H3K4me3[countOverlaps(ma.h3k4me1.anno,ma.h3k4me3)>0]<- TRUE
ma.h3k4me1.anno$H3K27ac<- FALSE
ma.h3k4me1.anno$H3K27ac[countOverlaps(ma.h3k4me1.anno,ma.h3k27ac)>0]<- TRUE
#promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)
ma.h3k4me1.anno$Promoter<- FALSE
ma.h3k4me1.anno$Promoter[ma.h3k4me1.anno$annotation=="Promoter (<=1kb)"]<- TRUE
ma.h3k4me1.anno$Promoter[ma.h3k4me1.anno$annotation=="Promoter (1-2kb)"]<- TRUE
ma.h3k4me1.anno$State<- "None"
ma.h3k4me1.anno$State[ma.h3k4me1.anno$H3K4me3]<- "H3K4me3"
ma.h3k4me1.anno$State[ma.h3k4me1.anno$H3K27ac]<- "H3K27ac"
ma.h3k4me1.anno$State[ma.h3k4me1.anno$H3K4me3 & ma.h3k4me1.anno$H3K27ac]<- "Both"
require(VennDiagram)
grid.newpage()
T<-venn.diagram(with(ma.h3k4me1.anno,list(H3K4me3=which(H3K4me3),Promoter=which(Promoter),H3K27ac=which(H3K27ac),H3K4me1=1:length(ma.h3k4me1.anno))),fill=c('darkorange', 'dodgerblue', 'hotpink', 'limegreen'), alpha=c(0.5,0.5,0.5,0.5), cex=2, filename=NULL)
grid.draw(T)
ma.h3k4me1.anno.s<- split(ma.h3k4me1.anno,ma.h3k4me1.anno$Promoter)

ma.h3k4me1.enh.anno<- annotatePeak(ma.h3k4me1.anno.s[[1]][,1], tssRegion=c(-2000, 2000), TxDb =txdb, annoDb="org.Hs.eg.db")
plotAnnoPie(ma.h3k4me1.enh.anno)
ma.h3k4me1.pro.anno<- annotatePeak(ma.h3k4me1.anno.s[[2]][,1], tssRegion=c(-2000, 2000), TxDb =txdb, annoDb="org.Hs.eg.db")
plotAnnoPie(ma.h3k4me1.pro.anno)

##relationship with gene expression
all.rpkm<- read.table("/mnt/local-disk1/rsgeno2/MAmotif/2.Processing/5.RNAseq_hg19_DEseq/2x75/rpkms_result/all_cells_rpkms_hg19.txt",row.names = NULL)
all.rpkm<- all.rpkm[!duplicated(all.rpkm$row.names),]
row.names(all.rpkm)<- all.rpkm$row.names
all.rpkm<- all.rpkm[,-1]

h3k4me1.rpkm<- all.rpkm[ma.h3k4me1.anno.s[[1]]$SYMBOL,13:14]
colnames(h3k4me1.rpkm)<- c("Rep1","Rep2")
h3k4me1.rpkm$state<- ma.h3k4me1.anno.s[[1]]$State

require(ggplot2)
require(reshape2)
ggplot(melt(h3k4me1.rpkm))+geom_boxplot(aes(x=variable,y=value,fill=state))+
  labs(x = "",y = "RPKM",title="The gene expression in different enhancer",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )+ ylim(1,60)


k562.array<- read.table("/mnt/local-disk1/rsgeno2/MAmotif/RACK7/H1ES_K562_HelaS3_hg19_new/K562_gene_expression.txt",header = TRUE,row.names = 1)
h3k4me1.array<- k562.array[ma.h3k4me1.anno.s[[1]]$SYMBOL,]
h3k4me1.array$state<- ma.h3k4me1.anno.s[[1]]$State
h3k4me1.array<- na.omit(h3k4me1.array)


ggplot(melt(h3k4me1.array))+geom_boxplot(aes(x=variable,y=value,fill=state))+
  labs(x = "",y = " ",title="The gene expression in different enhancer",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )+ ylim(1,500)+coord_flip()

## in promoter

h3k4me1.rpkm<- all.rpkm[ma.h3k4me1.anno.s[[2]]$SYMBOL,13:14]
colnames(h3k4me1.rpkm)<- c("Rep1","Rep2")
h3k4me1.rpkm$state<- ma.h3k4me1.anno.s[[2]]$State

ggplot(melt(h3k4me1.rpkm))+geom_boxplot(aes(x=variable,y=value,fill=state))+
  labs(x = "",y = "RPKM",title="The gene expression in different promoter",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )+ ylim(1,60)


h3k4me1.array<- k562.array[ma.h3k4me1.anno.s[[2]]$SYMBOL,]
h3k4me1.array$state<- ma.h3k4me1.anno.s[[2]]$State
h3k4me1.array<- na.omit(h3k4me1.array)


ggplot(melt(h3k4me1.array))+geom_boxplot(aes(x=variable,y=value,fill=state))+
  labs(x = "",y = " ",title="The gene expression in different promoter",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )+ ylim(1,500)+coord_flip()

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