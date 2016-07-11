require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(ChIPseeker)
require(clusterProfiler)
require(rtracklayer)
#MAnorm result between K562 and H1
k562.h1.manorm<- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me3/K562.H1.rep1.top2k.20160520/K562.H1.rep1.top2k.20160520_all_peak_MAvalues.xls",header = TRUE)
#k562.h1.manorm<- read.table("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me3/K562.H1.rep2.top2k.20160524/K562.H1.rep2.top2k.20160524_all_peak_MAvalues.xls",header = TRUE)
k562.h1.k27ac.manorm<- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K27ac/K562.H1.H3K27ac.rep2.top20k.20160525/K562.H1.H3K27ac.rep2.top20k.20160525_all_peak_MAvalues.xls",header = TRUE)

k562.h1.K4me1.manorm<- readPeakFile("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me1/K562.H1.H3K4me1.rep2.top20k.20160525/K562.H1.H3K4me1.rep2.top20k.20160525_all_peak_MAvalues.xls",header = TRUE)


k562.h1.manorm<- peakStatebylap(k562.h1.manorm,c(H3K27acK562=k562.h1.k27ac.manorm[k562.h1.k27ac.manorm$M_value>1],H3K27accomm=k562.h1.k27ac.manorm[abs(k562.h1.k27ac.manorm$M_value)<=1],H3K27acH1=k562.h1.k27ac.manorm[-1 > k562.h1.k27ac.manorm$M_value],
                                                 H3K4me1K562=k562.h1.K4me1.manorm[k562.h1.K4me1.manorm$M_value>1],H3K4me1comm=k562.h1.K4me1.manorm[abs(k562.h1.K4me1.manorm$M_value)<=1],H3K4me1H1=k562.h1.K4me1.manorm[-1 > k562.h1.K4me1.manorm$M_value],
                                                 Prom=promoters(txdb,upstream = 2000,downstream = 2000),K562SE=K562.SE,H1hescSE=H1.SE))


tp <- c(H3K4me1=k562.h1.k4me1.manorm,H3K27ac=k562.h1.k27ac.manorm,Prom=promoters(txdb,upstream = 2000,downstream = 2000),
           K562SE=K562.SE,H1hescSE=H1.SE)

k562.h1.manorm$K562<- k562.h1.manorm$M_value>=-1
k562.h1.manorm$H1hesc<- k562.h1.manorm$M_value<=1
k562.h1.manorm$H3K4me1[countOverlaps(k562.h1.manorm,k562.h1.k4me1.manorm)>0]<- TRUE
k562.h1.manorm$H3K27ac[countOverlaps(k562.h1.manorm,k562.h1.k27ac.manorm)>0]<- TRUE
k562.h1.manorm$Promoter[countOverlaps(k562.h1.manorm,promoters(txdb,upstream = 2000,downstream = 2000))>0]<- TRUE
k562.h1.manorm$monPromoter[countOverlaps(k562.h1.manorm,CTSS.prom)==0]<- TRUE

k562.h1.manorm.anno <- annotatePeak(k562.h1.manorm, tssRegion=c(-2000, 2000), TxDb =txdb, annoDb="org.Hs.eg.db")
k562.h1.manorm<- k562.h1.manorm.anno@anno

K562.SE<- import.bed("/mnt/local-disk1/rsgeno2/MAmotif/RACK7/routput/K562.bed")
k562.h1.manorm$K562SE[countOverlaps(k562.h1.manorm,K562.SE)>0]<- TRUE
H1.SE<- import.bed("/mnt/local-disk1/rsgeno2/MAmotif/RACK7/routput/H1.bed")
k562.h1.manorm$H1SE[countOverlaps(k562.h1.manorm,H1.SE)>0]<- TRUE


require(VennDiagram)
grid.newpage()
T<-venn.diagram(list(K562=which(k562.h1.manorm$K562),H1hesc=which(k562.h1.manorm$H1hesc),Promoter=which(k562.h1.manorm$Prom),
                     H3K27ac=which(k562.h1.manorm$H3K27ac),H3K4me1=which(k562.h1.manorm$H3K4me1)),fill=c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow'), alpha=c(0.5,0.5,0.5,0.5,0.5), cex=2, filename=NULL)
grid.draw(T)

grid.newpage()
T<-venn.diagram(list(K562=which(k562.h1.manorm$K562),H1hesc=which(k562.h1.manorm$H1hesc),Promoter=which(k562.h1.manorm$Prom)
                     ),fill=c('darkorange', 'dodgerblue', 'hotpink'), alpha=c(0.5,0.5,0.5), cex=2, filename=NULL)
grid.draw(T)

grid.newpage()
T<-venn.diagram(list(K562=which(k562.h1.manorm$K562 & !k562.h1.manorm$H1hesc),Promoter=which(k562.h1.manorm$Prom),
                     H3K27ac=which(k562.h1.manorm$H3K27ac),H3K4me1=which(k562.h1.manorm$H3K4me1)),fill=c('darkorange', 'dodgerblue', 'hotpink', 'limegreen'), alpha=c(0.5,0.5,0.5,0.5), cex=2, filename=NULL)
grid.draw(T)

grid.newpage()
T<-venn.diagram(list(H1hesc=which(k562.h1.manorm$H1hesc & !k562.h1.manorm$K562),Promoter=which(k562.h1.manorm$Prom),
                     H3K27ac=which(k562.h1.manorm$H3K27ac),H3K4me1=which(k562.h1.manorm$H3K4me1)),fill=c('darkorange', 'dodgerblue', 'hotpink', 'limegreen'), alpha=c(0.5,0.5,0.5,0.5), cex=2, filename=NULL)
grid.draw(T)




grid.newpage()
T<-venn.diagram(list(K562=which(k562.h1.manorm$K562),K562SE=which(k562.h1.manorm$K562SE),Promoter=which(k562.h1.manorm$Prom),
                     H3K27ac=which(k562.h1.manorm$H3K27ac),H3K4me1=which(k562.h1.manorm$H3K4me1)),fill=c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow'), alpha=c(0.5,0.5,0.5,0.5,0.5), cex=2, filename=NULL)
grid.draw(T)




k562.h1.manorm<- Peaks.Promoter(k562.h1.manorm,CTSS.promoter.merge)

grid.newpage()
T<-venn.diagram(list(K562=row.names(k562.h1.manorm)[k562.h1.manorm$K562],H1hesc=row.names(k562.h1.manorm)[k562.h1.manorm$H1hesc],Promoter=row.names(k562.h1.manorm)[k562.h1.manorm$state=="Promoter"]),fill=c('darkorange', 'dodgerblue', 'hotpink'), alpha=c(0.5,0.5,0.5), cex=2, filename=NULL)
grid.draw(T)
k562.h1.manorm$Gene.dist<- overlap.target.dist(k562.h1.manorm,enhancer.bed)
k562.h1.manorm$exp<- res[as.character(k562.h1.manorm$Gene.dist),2]
###8 state
k562.h1.manorm$group<- NA
k562.h1.manorm$group[k562.h1.manorm$state=="Intergenic" & k562.h1.manorm$K562 & !k562.h1.manorm$H1hesc]<- "K562\n-promoter\n-enhancer"
k562.h1.manorm$group[k562.h1.manorm$state=="Intergenic" & k562.h1.manorm$K562 & !k562.h1.manorm$H1hesc & k562.h1.manorm$H3K27ac>0 & k562.h1.manorm$H3K4me1]<- "K562\n-promoter\n+enhancer"
k562.h1.manorm$group[k562.h1.manorm$state=="Intergenic" & k562.h1.manorm$K562 & k562.h1.manorm$H1hesc]<- "Common\n-promoter\n-enhancer"
k562.h1.manorm$group[k562.h1.manorm$state=="Intergenic" & k562.h1.manorm$K562 & k562.h1.manorm$H1hesc & k562.h1.manorm$H3K27ac>0 & k562.h1.manorm$H3K4me1]<- "Common\n-promoter\n+enhancer"
k562.h1.manorm$group[k562.h1.manorm$state=="Intergenic" & !k562.h1.manorm$K562 & k562.h1.manorm$H1hesc]<- "H1hesc\n-promoter\n-enhancer"
k562.h1.manorm$group[k562.h1.manorm$state=="Intergenic" & !k562.h1.manorm$K562 & k562.h1.manorm$H1hesc & k562.h1.manorm$H3K27ac>0 & k562.h1.manorm$H3K4me1]<- "H1hesc\n-promoter\n+enhancer"

k562.h1.manorm$group[k562.h1.manorm$state=="Promoter" & k562.h1.manorm$K562 & !k562.h1.manorm$H1hesc]<- "K562\n+promoter\n-enhancer"
k562.h1.manorm$group[k562.h1.manorm$state=="Promoter" & k562.h1.manorm$K562 & !k562.h1.manorm$H1hesc & k562.h1.manorm$H3K27ac>0 & k562.h1.manorm$H3K4me1]<- "K562\n+promoter\n+enhancer"
k562.h1.manorm$group[k562.h1.manorm$state=="Promoter" & k562.h1.manorm$K562 & k562.h1.manorm$H1hesc]<- "Common\n+promoter\n-enhancer"
k562.h1.manorm$group[k562.h1.manorm$state=="Promoter" & k562.h1.manorm$K562 & k562.h1.manorm$H1hesc & k562.h1.manorm$H3K27ac>0 & k562.h1.manorm$H3K4me1]<- "Common\n+promoter\n+enhancer"
k562.h1.manorm$group[k562.h1.manorm$state=="Promoter" & !k562.h1.manorm$K562 & k562.h1.manorm$H1hesc]<- "H1hesc\n+promoter\n-enhancer"
k562.h1.manorm$group[k562.h1.manorm$state=="Promoter" & !k562.h1.manorm$K562 & k562.h1.manorm$H1hesc & k562.h1.manorm$H3K27ac>0 & k562.h1.manorm$H3K4me1]<- "H1hesc\n+promoter\n+enhancer"

boxplot(exp~group,data=k562.h1.manorm[k562.h1.manorm$state=="Promoter",])
boxplot(exp~group,data=k562.h1.manorm[k562.h1.manorm$state=="Intergenic",])
write.table(unique(k562.h1.manorm[k562.h1.manorm$state=="Intergenic" & k562.h1.manorm$K562 & !k562.h1.manorm$H1hesc,16]),file = "/mnt/local-disk1/rsgeno2/MAmotif/RACK7/routput/K562_unique_intergenic_peak_targenes_dist.txt",quote = FALSE,row.names = F,col.names = F)


k562.h1.manorm1<- k562.h1.manorm







peakfile = system.file("extdata", "sample_peaks.txt", package="ChIPseeker")
peakfile

peakAnno <- annotatePeak(peakfile, tssRegion=c(-3000, 1000), TxDb =txdb, annoDb="org.Hs.eg.db")

files <- getSampleFiles()
print(files)
peak <- readPeakFile(files[[4]])
peak

covplot(peak, weightCol="V5")
covplot(peak, weightCol="V5", chrs=c("chr17", "chr18"), xlim=c(4.5e7, 5e7))

promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)
tagMatrix <- getTagMatrix(peak, windows=promoter)
##
## to speed up the compilation of this vignettes, we use a precalculated tagMatrix
#data("tagMatrixList")
#tagMatrix <- tagMatrixList[[4]]
tagHeatmap(tagMatrix, xlim=c(-2000, 2000), color="red")
plotAvgProf(tagMatrix, xlim=c(-2000, 2000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
plotAvgProf(tagMatrix, xlim=c(-2000, 2000), conf = 0.95, resample = 500)
plotAnnoPie(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)

bp <- enrichGO(as.data.frame(peakAnno)$geneId, OrgDb='org.Hs.eg.db', ont="BP", readable=TRUE)
head(summary(bp), n=3)


#######mutiple chip peak data
data("tagMatrixList")
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))

## resample = 500 by default, here use 100 to speed up the compilation of this vignette.
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=100, facet="row")

compKEGG <- compareCluster(geneCluster   = as.data.frame(peakAnno)$geneId, 
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "BH")
plot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb, 
                       tssRegion=c(-3000, 3000), verbose=FALSE)