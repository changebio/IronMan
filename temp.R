require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(ChIPseeker)
require(clusterProfiler)

##read manorm results
ma.h3k4me1<-readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me1/K562.H3K4me1.both.top30k.20160607/K562.H3K4me1.both.top30k.20160607_all_peak_MAvalues.xls")
ma.h3k4me1<- ma.h3k4me1[ma.h3k4me1$M_value>= -1 & ma.h3k4me1$M_value<= 1]
ma.h3k4me3<- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me3/K562.H3K4me3.both.top20k.20160824/K562.H3K4me3.both.top20k.20160824_all_peak_MAvalues.xls")
ma.h3k4me3<- ma.h3k4me3[ma.h3k4me3$M_value>= -1 & ma.h3k4me3$M_value<= 1]
ma.h3k27ac<- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K27ac/K562.H3K27ac.both.top20k.20160607/K562.H3K27ac.both.top20k.20160607_all_peak_MAvalues.xls")
ma.h3k27ac<- ma.h3k27ac[ma.h3k27ac$M_value>= -1 & ma.h3k27ac$M_value<= 1]


##confident
ma.k4me3.fs<-list.files(path = "/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me3",pattern =".both.*?.xls",recursive = TRUE,full.names = TRUE)
ma.k4me3 <- lapply(ma.k4me3.fs,readPeakFile)
names(ma.k4me3)<-sapply(list.files(path = "/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me3",pattern =".both.*?.xls",recursive = TRUE),function(x)strsplit(x,split = "[.]")[[1]][1])
trimMA<- function(gr){
  gr<- gr[gr$M_value<1 & gr$M_value>-1]
  gr$score<- gr$summit
  return(gr)
}
ma.k4me3<- lapply(ma.k4me3,trimMA)
lapply(1:length(ma.k4me3),function(i)export.bedGraph(ma.k4me3[[i]],con = paste0("/mnt/local-disk1/rsgeno2/MAmotif/RACK7/routput/matrim/",names(ma.k4me3)[i],".trim.abs1.bed")))
                
                
                
                
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


############
Tfbs<- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/Track/wgEncodeRegTfbsClusteredV3.bed")
Tfbs<- sort(Tfbs)
Tfbs<- Tfbs[order(Tfbs$V4)]
Tfbs.anno<- annotatePeak(Tfbs, tssRegion=c(-2000, 2000), TxDb =txdb, annoDb="org.Hs.eg.db")


Tfbs.cell<- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/Track/wgEncodeRegTfbsClusteredWithCellsV3.bed")
Tfbs.cell <- sort(Tfbs.cell)
Tfbs.cell <- Tfbs.cell[order(Tfbs.cell$V4)]

Tfbs.K562 <- Tfbs.anno@anno[grep("K562",Tfbs.cell$V6)]


Dnase.cluster<- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/Track/wgEncodeRegDnaseClusteredV3.bed")


## Dscore of non-promoter H3K4me3

names(gro.seq)<- substring(names(gro.seq),first = 1,last = nchar(names(gro.seq))-9)
gro.cage.seq<- c(gro.seq,cage.seq)


dnase.K562.NB.gp<- region.cage.Dscore(dnase.K562.NB,cage.seq[1:4])


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

########

cage.seq.ds<- lapply(dnase.K562.anno.sl,function(x)region.cage.Dscore(x,cage.seq = cage.seq))
saveRDS(cage.seq.ds,file = "CAGE_seq_ds.rds")
gro.seq.ds<- lapply(dnase.K562.anno.sl,function(x)region.cage.Dscore(x,cage.seq = gro.seq))
saveRDS(gro.seq.ds,file = "GRO_seq_ds.rds")
all.seq.ds<- lapply(dnase.K562.anno.sl,function(x)region.cage.Dscore(x,cage.seq = list(all.gr.m,all.gr.p)))
saveRDS(all.seq.ds,file = "ALL_seq_ds.rds")

a<- Reduce(rbind,all.seq.ds)
a$group<- rep(names(all.seq.ds),as.numeric(sapply(all.seq.ds,nrow)))

ggplot(a)+geom_histogram(aes(x=value,y=..count..))+
  facet_wrap( ~ group)+
  labs(x = "",y = " ",title="Histogram of D score of Dnase groups in K562") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size=8)
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


###Dnase overlap
macs.dnase.f<- list.files("/mnt/local-disk1/rsgeno2/MAmotif/DNase/",pattern = "xls")[2:3]
macs.dnase<- lapply(macs.dnase.f,function(x)readPeakFile(paste0("/mnt/local-disk1/rsgeno2/MAmotif/DNase/",x)))
names(macs.dnase)<- macs.dnase.f
sapply(macs.dnase,function(x)table(countOverlaps(dnase.K562,x)))
enc.dnase.f<- list.files("/mnt/local-disk1/rsgeno2/MAmotif/ENCODE/2.DNase_Duke_hg19/",pattern = "bed")[1:4]
enc.dnase<- lapply(enc.dnase.f,function(x)import.bed(paste0("/mnt/local-disk1/rsgeno2/MAmotif/ENCODE/2.DNase_Duke_hg19/",x)))
names(enc.dnase)<- enc.dnase.f
sapply(enc.dnase, function(x)table(countOverlaps(dnase.K562,x)))

###Roadmap broadpeaks
bp.pk<- readRDS("/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/Roadmap_H3K4me3_BroadPeak_annotated.rds")
mcols(bp.pk$`ES-I3_Cell_Line`)[,c(3,)]
