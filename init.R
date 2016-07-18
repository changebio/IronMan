

# setting color -----------------------------------------------------------

tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
par(pch=19)
plot(c(1,2,3),c(1,2,3),col=11)
# # load refgene hg19 -------------------------------------------------------
# 
# load("/opt/rstudio/rstudio1/rstudio1/proj/refdata/refhg19.RData")
# gtf<- read.table("/mnt/local-disk1/rsgeno2/huangyin/PRC2/Refgenes/genes.gtf",sep = "\t")

#merge multiple version promoter
ref.name<- c("bin","name","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount","exonStarts","exonEnds","score","name2","cdsStartStat","cdsEndStat","exonFrames")
ref.files<- list.files("/mnt/local-disk1/rsgeno2/MAmotif/UCSC_brower_2016-2-17", pattern = "txt")
hg19_refGene<- lapply(paste0("/mnt/local-disk1/rsgeno2/MAmotif/UCSC_brower_2016-2-17/",ref.files), function(x)read.table(x,sep="\t"))
names(hg19_refGene)<- substr(ref.files,1,nchar(ref.files)-4)
hg19_refGene$knownGene[,2:13]<- hg19_refGene$knownGene
# refGene<- paste0(hg19_refGene$refGene$V3,hg19_refGene$refGene$V4,hg19_refGene$refGene$V5)
# knowGene<- paste0(hg19_refGene$knownGene$V2,hg19_refGene$knownGene$V3,hg19_refGene$knownGene$V4)
# gencode<- paste0(hg19_refGene$hg19_GENCODE_V14$V3,hg19_refGene$hg19_GENCODE_V14$V4,hg19_refGene$hg19_GENCODE_V14$V5)
# ensembl<- paste0(hg19_refGene$hg19_Ensembl_Genes$V3,hg19_refGene$hg19_Ensembl_Genes$V4,hg19_refGene$hg19_Ensembl_Genes$V5)
# grid.newpage()
# T<-venn.diagram(list(RefGene=refGene,KnowGene=knowGene,GENCODE=gencode,Ensembl=ensembl),fill=c('darkorange', 'dodgerblue', 'hotpink', 'limegreen'), alpha=c(0.5,0.5,0.5,0.5), cex=2, filename=NULL)
# grid.draw(T)

#combine the TSS from different refgene 
CTSS<- rbind(hg19_refGene$knownGene[,3:6],hg19_refGene$hg19_Ensembl_Genes[,3:6],hg19_refGene$hg19_GENCODE_V14[,3:6],hg19_refGene$refGene[,3:6])
CTSS<- CTSS[!duplicated(CTSS),]
CTSS<- CTSS[nchar(as.character(CTSS$V3))<6,]
CTSS<- CTSS[,c(1,3,4,2)]
colnames(CTSS)<- c("chr","start","end","strand")
CTSS<- makeGRangesFromDataFrame(CTSS)
start(CTSS)<- start(CTSS) + 1
CTSS.prom<- promoters(CTSS,upstream = 2000,downstream = 2000)
saveRDS(CTSS.prom,file = "data/CTSS_promoters.rds")
# pstart<- CTSS$V5-2000
# pstart[CTSS$V4=="-"]<- CTSS$V6[CTSS$V4=="-"]-2000
# pend<- CTSS$V5+2000
# pend[CTSS$V4=="-"]<- CTSS$V6[CTSS$V4=="-"]+2000
# CTSS.promoter<- data.frame(chr=CTSS$V3,start=pstart,end=pend)
# CTSS.promoter.merge<- merge.peaks(CTSS.promoter)
# load hg19 cpg island  ---------------------------------------------------

cpghg19<- read.table("/mnt/local-disk1/rsgeno2/huangyin/PRC2/hg19.cpgIslandExt.txt",header=TRUE,row.names = NULL,sep = "\t")
cpghg19<- cpghg19[,-1]

# load two H3K4me3 ChIP Seq data from ENCODE ------------------------------
h3k4me3.files <- list.files("/mnt/local-disk1/rsgeno2/MAmotif/macs.20151224/H3K4me3", pattern = "xls",include.dirs = TRUE)
h3k4me3<- lapply(paste0("/mnt/local-disk1/rsgeno2/MAmotif/macs.20151224/H3K4me3/",h3k4me3.files),function(x)read.table(x,header = TRUE))
names(h3k4me3)<- substr(h3k4me3.files,1,nchar(h3k4me3.files)-24)
h3k4me3<- lapply(h3k4me3,function(x)Peaks.Distribution(x,refhg19))
save(h3k4me3, file = "/opt/rstudio/rstudio1/rstudio1/proj/IronMan/data/h3k4me3.Rdata")

h3k4me3.roadmap.files <- list.files("/mnt/local-disk1/rsgeno2/huangyin/RoadMap/MAmotif/1.Histone_Roadmap_hg19_peaks/H3K4me3", pattern = "xls",include.dirs = TRUE)
h3k4me3.roadmap<- lapply(paste0("/mnt/local-disk1/rsgeno2/huangyin/RoadMap/MAmotif/1.Histone_Roadmap_hg19_peaks/H3K4me3/",h3k4me3.roadmap.files),function(x)read.table(x,header = TRUE))
names(h3k4me3.roadmap)<- substr(h3k4me3.roadmap.files,1,nchar(h3k4me3.roadmap.files)-24)
h3k4me3.roadmap<- lapply(h3k4me3.roadmap,function(x)Peaks.Distribution(x,refhg19))
h3k4me3.ctss<- lapply(h3k4me3,function(x)Peaks.Promoter(x,CTSS.promoter.merge))
# load two H3K27ac ChIP Seq data from ENCODE ------------------------------
h3k27ac.files <- list.files("/mnt/local-disk1/rsgeno2/MAmotif/macs.20151224/H3K27ac", pattern = "xls",include.dirs = TRUE)
h3k27ac<- lapply(paste0("/mnt/local-disk1/rsgeno2/MAmotif/macs.20151224/H3K27ac/",h3k27ac.files),function(x)read.table(x,header = TRUE))
names(h3k27ac)<- substr(h3k27ac.files,1,nchar(h3k27ac.files)-24)
h3k27ac<- lapply(h3k27ac,function(x)Peaks.Distribution(x,refhg19))
save(h3k27ac, file = "/opt/rstudio/rstudio1/rstudio1/proj/IronMan/data/h3k27ac.Rdata")
# load two H3K4me1 ChIP Seq data from ENCODE ------------------------------
h3k4me1.files <- list.files("/mnt/local-disk1/rsgeno2/MAmotif/macs.20151224/H3K4me1", pattern = "xls",include.dirs = TRUE)
h3k4me1<- lapply(paste0("/mnt/local-disk1/rsgeno2/MAmotif/macs.20151224/H3K4me1/",h3k4me1.files),function(x)read.table(x,header = TRUE))
names(h3k4me1)<- substr(h3k4me1.files,1,nchar(h3k4me1.files)-24)
h3k4me1<- lapply(h3k4me1,function(x)Peaks.Distribution(x,refhg19))
save(h3k4me1, file = "/opt/rstudio/rstudio1/rstudio1/proj/IronMan/data/h3k4me1.Rdata")
# read RNA Seq data of K562 and H1
H1.K562.exp<- read.table("/mnt/local-disk1/rsgeno2/MAmotif/2.Processing/5.RNAseq_hg19_DEseq/2x75/read_count/K562_read_count_H1hesc_read_count.txt")
Nhek.exp<-read.table("/mnt/local-disk1/rsgeno2/MAmotif/2.Processing/5.RNAseq_hg19_DEseq/2x75/read_count/Nhek_read_count.txt")
countData<- H1.K562.exp
colData<- data.frame(condition=c("K562","K562","H1","H1","H1","H1"),type="paired-end")
row.names(colData)<- colnames(H1.K562.exp)
library("DESeq2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.2")
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ condition)
dds<- DESeq(dds)
res <- results(dds)
res <- as.data.frame(res)
plotMA(res, main="DESeq2", ylim=c(-2,2))
res <- results(dds, contrast=c("condition","K562","H1"))
promoter.bed$exp<- res[promoter.bed$name,2]

library("DESeq")	

# give the sample label
Nrep1 = 2
Nrep2 = 4
condition = c(rep("K562", times=Nrep1), rep("H1", times=Nrep2))
countTable<- H1.K562.exp
cds = newCountDataSet(countTable,condition)

# perfor normalization #
cds = estimateSizeFactors(cds)
sizeFactors(cds)
head(counts(cds,normalized=TRUE))
normalized_cds=counts(cds,normalized=TRUE)

# estimate variance #
if(Nrep1==1 || Nrep2==1) cds=estimateDispersions(cds,method='blind',sharingMode="fit-only") else cds=estimateDispersions(cds)
#plotDispEsts(cds)

# call differential expression
res = nbinomTest(cds,"H1","K562")
res = cbind(res,normalized_cds)  #add the normalized read count.

#load MAnorm2 result
Mar6v4<- read.csv("/mnt/local-disk1/rsgeno2/MAmotif/RACK7/Mar6v4H3K4me3_merged_peak_bins_info.csv")
Mar6v4<- Peaks.Distribution(Mar6v4,refhg19)
Mar6v4$Genes<- overlap.target.exp(Mar6v4,promoter.bed,col = 4,exp = promoter.bed$exp)
Mar6v4$Genes2<- overlap.target.dist(Mar6v4,promoter.bed)
Mar6v4$eGenes<- overlap.target.exp(Mar6v4,enhancer.bed,exp=enhancer.bed$exp)
Mar6v4$eGenes2<- overlap.target.dist(Mar6v4,enhancer.bed)
Mar6v4.pro<- Mar6v4[Mar6v4$state=="Promoter",]
Mar6v4.pro$Exp<- res[as.character(Mar6v4.pro$Genes),2]
Mar6v4.pro$Exp2<- res[as.character(Mar6v4.pro$Genes2),2]
Mar6v4.pro$FC<- log2(rowMeans(Mar6v4.pro[,26:31])/rowMeans(Mar6v4.pro[,32:35]))
Mar6v4.pro$Genes<- as.character(Mar6v4.pro$Genes)
Mar6v4.pro$Genes2<- as.character(Mar6v4.pro$Genes2)
plot(Mar6v4.pro$FC,Mar6v4.pro$Exp,xlab="M value",ylab="log2 FC expression")
abline(lm(Mar6v4.pro$FC~Mar6v4.pro$Exp))
Mar6v4.pro.uni<- Mar6v4.pro[,38:44]
Mar6v4.pro.uni<- Mar6v4.pro.uni[order(Mar6v4.pro.uni$Genes,- Mar6v4.pro.uni$FC,- Mar6v4.pro.uni$Exp),]
Mar6v4.pro.dep<- Mar6v4.pro.uni[!duplicated(Mar6v4.pro.uni$Genes),]
plot(Mar6v4.pro.dep$FC,Mar6v4.pro.dep$Exp)
summary(lm(Mar6v4.pro.dep$Exp~Mar6v4.pro.dep$FC))
abline(lm(Mar6v4.pro.dep$Exp~Mar6v4.pro.dep$FC))
Mar6v4.enh<- Mar6v4[Mar6v4$state!="Promoter",]
Mar6v4.enh$Exp<- res[as.character(Mar6v4.enh$eGenes),2]
Mar6v4.enh$Exp2<- res[as.character(Mar6v4.enh$eGenes2),2]
Mar6v4.enh$FC<- log2(rowMeans(Mar6v4.enh[,26:31])/rowMeans(Mar6v4.enh[,32:35]))
Mar6v4.enh$eGenes<- as.character(Mar6v4.enh$eGenes)
Mar6v4.enh$eGenes2<- as.character(Mar6v4.enh$eGenes2)
plot(Mar6v4.enh$FC,Mar6v4.enh$Exp)
summary(lm(Mar6v4.enh$FC~Mar6v4.enh$Exp))
Mar6v4.enh.uni<- Mar6v4.enh[,38:44]
Mar6v4.enh.uni<- Mar6v4.enh.uni[order(Mar6v4.enh.uni$eGenes,-Mar6v4.enh.uni$FC,-Mar6v4.enh.uni$Exp),]
Mar6v4.enh.dep<- Mar6v4.enh.uni[!duplicated(Mar6v4.enh.uni$eGenes),]
plot(Mar6v4.enh.dep$FC,Mar6v4.enh.dep$Exp)
summary(lm(Mar6v4.enh.dep$Exp~Mar6v4.enh.dep$FC))
abline(lm(Mar6v4.pro.dep$Exp~Mar6v4.pro.dep$FC))


#MAnorm result between K562 and H1
k562.h1.manorm<- read.table("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me3/K562.H1.rep1.top2k.20160520/K562.H1.rep1.top2k.20160520_all_peak_MAvalues.xls",header = TRUE)
#k562.h1.manorm<- read.table("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me3/K562.H1.rep2.top2k.20160524/K562.H1.rep2.top2k.20160524_all_peak_MAvalues.xls",header = TRUE)
k562.h1.k27ac.manorm<- read.table("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K27ac/K562.H1.H3K27ac.rep2.top20k.20160525/K562.H1.H3K27ac.rep2.top20k.20160525_all_peak_MAvalues.xls",header = TRUE)
k562.h1.k4me1.manorm<- read.table("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me1/K562.H1.H3K4me1.rep2.top20k.20160525/K562.H1.H3K4me1.rep2.top20k.20160525_all_peak_MAvalues.xls",header = TRUE)
#k562.h1.manorm<- Peaks.Promoter(k562.h1.manorm,refhg19$PROMOTER)
k562.h1.manorm<- Peaks.Promoter(k562.h1.manorm,CTSS.promoter.merge)
k562.h1.manorm$K562<- k562.h1.manorm$M_value>=-1
k562.h1.manorm$H1hesc<- k562.h1.manorm$M_value<=1
k562.h1.manorm$H3K27ac<- overlap(k562.h1.manorm,k562.h1.k27ac.manorm)[[1]]
k562.h1.manorm$H3K4me1<- overlap(k562.h1.manorm,k562.h1.k4me1.manorm)[[1]]
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

grid.newpage()
T<-venn.diagram(list(K562=row.names(k562.h1.manorm)[k562.h1.manorm$K562],H1hesc=row.names(k562.h1.manorm)[k562.h1.manorm$H1hesc],Promoter=row.names(k562.h1.manorm)[k562.h1.manorm$state=="Promoter"],
                     H3K27ac=row.names(k562.h1.manorm)[k562.h1.manorm$H3K27ac>0],H3K4me1=row.names(k562.h1.manorm)[k562.h1.manorm$H3K4me1>0]),fill=c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow'), alpha=c(0.5,0.5,0.5,0.5,0.5), cex=2, filename=NULL)
grid.draw(T)

k562.h1.manorm1<- k562.h1.manorm

#######K562-unique non-promoter H3K4me3 peaks in Rep1 and Rep2
k562.rep1.uqe<- k562.h1.manorm1[k562.h1.manorm1$state=="Intergenic" & k562.h1.manorm1$K562 & !k562.h1.manorm1$H1hesc,]
k562.rep2.uqe<- k562.h1.manorm[k562.h1.manorm$state=="Intergenic" & k562.h1.manorm$K562 & !k562.h1.manorm$H1hesc,]
###The common target genes of K562-unique non-promoter peaks  two H3K4me3 Rep 
grid.newpage()
T<-venn.diagram(list(Rep1=k562.rep1.uqe$Gene.dist,Rep2=k562.rep2.uqe$Gene.dist),fill=c('darkorange', 'dodgerblue'), alpha=c(0.5,0.5), cex=2, filename=NULL)
grid.draw(T)
###The relationship between H3K27ac and H3K4me1 in K562-unique non-promoter peaks
grid.newpage()
T<-venn.diagram(list(H3K4me3=row.names(k562.rep1.uqe),H3K27ac=row.names(k562.rep1.uqe)[k562.rep1.uqe$H3K27ac>0],H3K4me1=row.names(k562.rep1.uqe)[k562.rep1.uqe$H3K4me1>0]),fill=c('darkorange', 'dodgerblue', 'hotpink'), alpha=c(0.5,0.5,0.5), cex=2, filename=NULL)
grid.draw(T)
fisher.test(table(k562.rep1.uqe$H3K27ac>0,k562.rep1.uqe$H3K4me1>0))
grid.newpage()
T<-venn.diagram(list(H3K4me3=row.names(k562.rep2.uqe),H3K27ac=row.names(k562.rep2.uqe)[k562.rep2.uqe$H3K27ac>0],H3K4me1=row.names(k562.rep2.uqe)[k562.rep2.uqe$H3K4me1>0]),fill=c('darkorange', 'dodgerblue', 'hotpink'), alpha=c(0.5,0.5,0.5), cex=2, filename=NULL)
grid.draw(T)
fisher.test(table(k562.rep2.uqe$H3K27ac>0,k562.rep2.uqe$H3K4me1>0))
####K562-unique enhancer H3K4me3 peaks in Rep1 and Rep2
k562.rep1.enh<- k562.rep1.uqe[k562.rep1.uqe$H3K27ac>0 & k562.rep1.uqe$H3K4me1>0,]
k562.rep2.enh<- k562.rep2.uqe[k562.rep2.uqe$H3K27ac>0 & k562.rep2.uqe$H3K4me1>0,]
grid.newpage()
T<-venn.diagram(list(Rep1=k562.rep1.enh$Gene.dist,Rep2=k562.rep2.enh$Gene.dist),fill=c('darkorange', 'dodgerblue'), alpha=c(0.5,0.5), cex=2, filename=NULL)
grid.draw(T)
####GO analysis for K562-unique target genes
write.table(intersect(unlist(k562.rep1.enh$Gene.dist),unlist(k562.rep2.enh$Gene.dist)),file = "/mnt/local-disk1/rsgeno2/MAmotif/RACK7/routput/K562_unique_enhancer_targenes_by_two_rep.txt",quote = FALSE,row.names = F,col.names = F)

####H1hesc-unique non-promoter H3K4me3 peaks in Rep1
h1.rep1.uqe<- k562.h1.manorm[k562.h1.manorm$state=="Intergenic" & !k562.h1.manorm$K562 & k562.h1.manorm$H1hesc,]
h1.rep1.enh<- h1.rep1.uqe[h1.rep1.uqe$H3K27ac>0 & h1.rep1.uqe$H3K4me1>0,]
plot(h1.rep1.enh$M_value,h1.rep1.enh$exp)
abline(h = 0)
boxplot(h1.rep1.enh$exp)
summary(lm(h1.rep1.enh$exp~h1.rep1.enh$M_value))
grid.newpage()
T<-venn.diagram(list(H3K4me3=row.names(h1.rep1.uqe),H3K27ac=row.names(h1.rep1.uqe)[h1.rep1.uqe$H3K27ac>0],H3K4me1=row.names(h1.rep1.uqe)[h1.rep1.uqe$H3K4me1>0]),fill=c('darkorange', 'dodgerblue', 'hotpink'), alpha=c(0.5,0.5,0.5), cex=2, filename=NULL)
grid.draw(T)
fisher.test(table(h1.rep1.uqe$H3K27ac>0,h1.rep1.uqe$H3K4me1>0))


####the relationship between K562-unique enhancer H3K4me3 and the expression of their target gene
k562.rep1.enh.name<- as.character(k562.rep1.enh$Gene.dist)
k562.rep1.enh.name[!(as.character(k562.rep1.enh$Gene.dist) %in% row.names(res))]<- NA
k562.rep1.enh$exp<- res[k562.rep1.enh.name,2]
#oncogenes
comic<- read.csv("/mnt/local-disk1/rsgeno2/MAmotif/RACK7/routput/Census_allFri_May27_2016.csv")
comic<- comic[comic$Molecular.Genetics=="Dom",]
comic<- comic[comic$Role.in.Cancer!="TSG",]
intersect(k562.h1.manorm[k562.h1.manorm$state=="Intergenic" & k562.h1.manorm$K562 & !k562.h1.manorm$H1hesc,16],comic$Gene.Symbol)


##How to prove overactivation?
Histone.state<- function(cell.lines,other.lines=NULL){
  temp<- merge.peaks(do.call(rbind.data.frame,cell.lines))
  temp1<- data.frame(sapply(c(cell.lines,other.lines),function(x)overlap.target(temp,x,col = 10)[[1]]))
  temp1$gene<- overlap.target.dist(temp,enhancer.bed)
  return(cbind(temp,temp1))
}
h3k4me1.lap<- Histone.state(list(Rep1=h3k4me1.top$K562.H3k4me1.Rep1,Rep2=h3k4me1.top$K562.H3k4me1.Rep2),list(H3K27ac=h3k27ac.top$K562.H3k27ac.Rep1,H3K4me3=h3k4me3.top$K562.H3k4me3.Rep1))
h3k4me1.enh<- h3k4me1.lap[!is.na(h3k4me1.lap$Rep1) & !is.na(h3k4me1.lap$Rep2),]
h3k4me1.enh<- h3k4me1.enh[h3k4me1.enh$Rep1!="Promoter" & h3k4me1.enh$Rep2!="Promoter",]
h3k4me1.enh$state<- "Both"
h3k4me1.enh$state[is.na(h3k4me1.enh$H3K27ac)]<- "None"
h3k4me1.enh$state[is.na(h3k4me1.enh$H3K27ac) & !is.na(h3k4me1.enh$H3K4me3)]<- "H3K4me3"
h3k4me1.enh$state[!is.na(h3k4me1.enh$H3K27ac) & is.na(h3k4me1.enh$H3K4me3)]<- "H3K27ac"

h3k4me1.exp<- log2(H1.K562.exp[as.character(h3k4me1.enh$gene),1:2])
colnames(h3k4me1.exp)<- c("Rep1","Rep2")
h3k4me1.exp$state<- h3k4me1.enh$state


ggplot(melt(h3k4me1.exp))+geom_boxplot(aes(x=variable,y=value,fill=state))+
  labs(x = "",y = "log2-transformed read count",title="The gene expression in different enhancer",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )

all.rpkm<- read.table("/mnt/local-disk1/rsgeno2/MAmotif/2.Processing/5.RNAseq_hg19_DEseq/2x75/rpkms_result/all_cells_rpkms_hg19.txt",row.names = NULL)
all.rpkm<- all.rpkm[!duplicated(all.rpkm$row.names),]
row.names(all.rpkm)<- all.rpkm$row.names
all.rpkm<- all.rpkm[,-1]

h3k4me1.rpkm<- all.rpkm[as.character(h3k4me1.enh$gene),13:14]
colnames(h3k4me1.rpkm)<- c("Rep1","Rep2")
h3k4me1.rpkm$state<- h3k4me1.enh$state


ggplot(melt(h3k4me1.rpkm))+geom_boxplot(aes(x=variable,y=value,fill=state))+
  labs(x = "",y = "RPKM",title="The gene expression in different enhancer",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )+ ylim(1,25)

h3k4me1.m<- read.table("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me1/K562.H3K4me1.both.top20k.20160602/K562.H3K4me1.both.top20k.20160602_all_peak_MAvalues.xls",header = TRUE)
h3k4me1.m<- h3k4me1.m[with(h3k4me1.m,order(chr,start,end)),]
h3k4me1.m<- h3k4me1.m[h3k4me1.m$M_value>= -1 & h3k4me1.m$M_value<= 1,-10]
h3k4me1.m<- Peaks.Promoter(h3k4me1.m,CTSS.promoter.merge)
h3k4me1.m<- h3k4me1.m[h3k4me1.m$state=="Intergenic",]
h3k4me1.m<- Histone.state(list(H3K4me1=h3k4me1.m),list(H3K27ac=h3k27ac.top$K562.H3k27ac.Rep1,H3K4me3=h3k4me3.top$K562.H3k4me3.Rep1))
h3k4me1.m$state<- "Both"
h3k4me1.m$state[is.na(h3k4me1.m$H3K27ac)]<- "None"
h3k4me1.m$state[is.na(h3k4me1.m$H3K27ac) & !is.na(h3k4me1.m$H3K4me3)]<- "H3K4me3"
h3k4me1.m$state[!is.na(h3k4me1.m$H3K27ac) & is.na(h3k4me1.m$H3K4me3)]<- "H3K27ac"

h3k4me1.exp<- log2(H1.K562.exp[as.character(h3k4me1.m$gene),1:2])
colnames(h3k4me1.exp)<- c("Rep1","Rep2")
h3k4me1.exp$state<- h3k4me1.m$state
######
h3k4me1.m<- read.table("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me1/K562.H3K4me1.both.top30k.20160607/K562.H3K4me1.both.top30k.20160607_all_peak_MAvalues.xls",header = TRUE)
h3k4me1.m<- h3k4me1.m[with(h3k4me1.m,order(chr,start,end)),]
h3k4me1.m<- h3k4me1.m[h3k4me1.m$M_value>= -1 & h3k4me1.m$M_value<= 1,-10]
h3k4me1.m<- Peaks.Promoter(h3k4me1.m,CTSS.promoter.merge)
hist(h3k4me1.m$M_value)
ReadMacs.State<- function(fp,ptr){
  temp.m<- read.table(fp,header = TRUE)
  temp.m<- temp.m[with(temp.m,order(chr,start,end)),]
  temp.m<- temp.m[temp.m$M_value>= -1 & temp.m$M_value<= 1,-10]
  temp.m<- Peaks.Promoter(temp.m,ptr)
  return(temp.m)
}
h3k4me1.m<- h3k4me1.m[h3k4me1.m$state=="Intergenic",]
hist(h3k4me1.m$M_value)
hist(h3k4me1.m$normalized_read_density_in_wgEncodeBroadHistoneK562H3k4me1StdAlnRep1)
h3k4me1.ms<- Histone.state(list(H3K4me1=h3k4me1.m),list(H3K27ac=ReadMacs.State("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K27ac/K562.H3K27ac.both.top20k.20160607/K562.H3K27ac.both.top20k.20160607_all_peak_MAvalues.xls",CTSS.promoter.merge),
                                                       H3K4me3=ReadMacs.State("/mnt/local-disk1/rsgeno2/MAmotif/manorm.20160520/H3K4me3/K562.H3K4me3.both.top20k.20160607/K562.H3K4me3.both.top20k.20160607_all_peak_MAvalues.xls",CTSS.promoter.merge)))
h3k4me1.m<- cbind(h3k4me1.m,h3k4me1.ms[,4:7])
h3k4me1.m$state<- "Both"
h3k4me1.m$state[is.na(h3k4me1.m$H3K27ac)]<- "None"
h3k4me1.m$state[is.na(h3k4me1.m$H3K27ac) & !is.na(h3k4me1.m$H3K4me3)]<- "H3K4me3"
h3k4me1.m$state[!is.na(h3k4me1.m$H3K27ac) & is.na(h3k4me1.m$H3K4me3)]<- "H3K27ac"
boxplot(h3k4me1.m$normalized_read_density_in_wgEncodeBroadHistoneK562H3k4me1StdAlnRep1~h3k4me1.m$state)

h3k4me1.gr<- makeGRangesFromDataFrame(h3k4me1.m,keep.extra.columns = TRUE)
h3k4me1.exp<- log2(H1.K562.exp[as.character(h3k4me1.m$gene),1:2])
colnames(h3k4me1.exp)<- c("Rep1","Rep2")
h3k4me1.exp$state<- h3k4me1.m$state

ggplot(melt(h3k4me1.exp))+geom_boxplot(aes(x=variable,y=value,fill=state))+
  labs(x = "",y = "log2-transformed read count",title="The gene expression in different enhancer",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )

all.rpkm<- read.table("/mnt/local-disk1/rsgeno2/MAmotif/2.Processing/5.RNAseq_hg19_DEseq/2x75/rpkms_result/all_cells_rpkms_hg19.txt",row.names = NULL)
all.rpkm<- all.rpkm[!duplicated(all.rpkm$row.names),]
row.names(all.rpkm)<- all.rpkm$row.names
all.rpkm<- all.rpkm[,-1]

h3k4me1.rpkm<- all.rpkm[as.character(h3k4me1.m$gene),13:14]
colnames(h3k4me1.rpkm)<- c("Rep1","Rep2")
h3k4me1.rpkm$state<- h3k4me1.m$state


ggplot(melt(h3k4me1.rpkm))+geom_boxplot(aes(x=variable,y=value,fill=state))+
  labs(x = "",y = "RPKM",title="The gene expression in different enhancer",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )+ ylim(1,25)


k562.array<- read.table("/mnt/local-disk1/rsgeno2/MAmotif/RACK7/H1ES_K562_HelaS3_hg19_new/K562_gene_expression.txt",header = TRUE,row.names = 1)
h3k4me1.array<- k562.array[as.character(h3k4me1.m$gene),]
h3k4me1.array$state<- h3k4me1.m$state





h3k4me1.rpkm<- k562.array[as.character(temp3$gene),]
colnames(h3k4me1.rpkm)<- c("Rep1","Rep2")
h3k4me1.rpkm$state<- temp3$state


ggplot(melt(h3k4me1.rpkm))+geom_boxplot(aes(x=variable,y=value,fill=state))+
  labs(x = "",y = "RPKM",title="The gene expression in different enhancer",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )+ ylim(1,25)