
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
#saveRDS(dnase.K562.anno,file = "data/Dnase_K562_anno.rds")

require(VennDiagram)
grid.newpage()
T<-venn.diagram(with(dnase.K562.anno,list(H3K4me3=which(H3K4me3=="None"),Promoter=which(Promoter),H3K27ac=which(H3K27ac=="None"),H3K4me1=which(H3K4me1=="None"))),fill=c('darkorange', 'dodgerblue', 'hotpink', 'limegreen'), alpha=c(0.5,0.5,0.5,0.5), cex=2, filename=NULL)
grid.draw(T)

dnase.K562.anno.s<- resize(dnase.K562.anno,width = 500,fix = "center")
dnase.K562.anno.sl<- split(dnase.K562.anno.s,dnase.K562.anno.s$Group)
dnase.K562.anno.sl<- dnase.K562.anno.sl[sapply(dnase.K562.anno.sl,function(x)length(x)>500)]



##Cage coverage--------------

cage.seq.cov<-sapply(1:length(cage.seq),function(i){
  cage.cov<- sapply(dnase.K562.anno.sl,function(x)ScoreMatrixBin(cage.seq[[i]],windows = x,bin.num = 1,weight.col = "V5")@.Data*500)
  return(sapply(cage.cov, function(x)sum(x>3)/length(x)))
})


cage.seq.cov<- as.data.frame(cage.seq.cov)
colnames(cage.seq.cov)<- names(cage.seq)
saveRDS(cage.seq.cov,file = "data/CAGE_Seq_cov.rds")
require(ggplot2)
require(reshape2)
cage.seq.cov$region<- row.names(cage.seq.cov)
gd<- melt(cage.seq.cov)
ggplot(gd)+geom_bar(aes(x=region,y=value,fill=variable),stat = "identity",position = "dodge")+coord_flip()+
  labs(x = "",y = " ",title="The CAGE distribution in the whole genome",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )


gro.seq.cov<-sapply(1:length(gro.seq),function(i){
  gro.cov<- sapply(dnase.K562.anno.sl,function(x)ScoreMatrixBin(gro.seq[[i]],windows = x,bin.num = 1,weight.col = "V5")@.Data*500)
  return(sapply(gro.cov, function(x)sum(x>3)/length(x)))
})
saveRDS(gro.seq.cov,file = "data/GRO_Seq_cov.rds")
gro.seq.cov<- as.data.frame(gro.seq.cov)
colnames(gro.seq.cov)<- substring(names(gro.seq),first = 1,last = nchar(names(gro.seq))-9)


all.seq<-list(all.gr.m,all.gr.p)
all.seq.cov<-sapply(1:length(all.seq),function(i){
  all.cov<- sapply(dnase.K562.anno.sl,function(x)ScoreMatrixBin(all.seq[[i]],windows = x,bin.num = 1,weight.col = "V5")@.Data*500)
  return(sapply(all.cov, function(x)sum(x>3)/length(x)))
})
saveRDS(all.seq.cov,file = "data/all_seq_cov.rds")
all.seq.cov<- as.data.frame(all.seq.cov)
colnames(all.seq.cov)<- c("CTSS_all_minus_pool","CTSS_all_plus_pool")
all.gro.cage.cov<- cbind(cage.seq.cov,gro.seq.cov)
all.gro.cage.cov<- cbind(all.gro.cage.cov,all.seq.cov)
saveRDS(all.gro.cage.cov,file = "data/all_gro_cage_cov.rds")

ggplot(melt(all.gro.cage.cov[,27:33]))+geom_bar(aes(x=region,y=value,fill=variable),stat = "identity",position = "dodge")+coord_flip()+
  labs(x = "",y = " ",title="The CAGE distribution in the whole genome (FANTOM)",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )

#======================================

##TF coverage--------------

require(rtracklayer)
tf.bw.files <- list.files("/mnt/local-disk1/rsgeno2/MAmotif/ENCODE/Tfbs/SydhTfbsK562/",pattern = "bigWig")
tf.bw.files.uni<- tf.bw.files[!duplicated(substring(tf.bw.files,first = 1,last=24))]
tf.bw.sbin<- lapply(21:40,function(i,tf){
  tf.bw<- import.bw(paste0("/mnt/local-disk1/rsgeno2/MAmotif/ENCODE/Tfbs/SydhTfbsK562/",tf[i]),asRle = TRUE)
  tf.cov<- lapply(dnase.K562.anno.sl,function(x)ScoreMatrixBin(tf.bw,windows = x,bin.num=50))
  return(tf.cov)},tf=tf.bw.files.uni)

names(tf.bw.sbin)<- tf.bw.files.uni[21:40]
#saveRDS(tf.bw.sbin,file = "/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/data/TF_bw_sbin_uni21-40.rds")


#saveRDS(tf.bw.sbin,file = "data/TF_bw_sbin_21-25.rds")


tf.bw.sbin.m<- 

Dnase.TFbin.ms<- function(tf.bw.sbin,scale=FALSE){
  tf.bw.sbin.m<- lapply(1:length(tf.bw.sbin),function(i,tf){
    if(scale){
      temp<- melt(as.data.frame(sapply(tf[[i]],function(x)colMeans(scaleScoreMatrix(x)))))
    }else{
      temp<- melt(as.data.frame(sapply(tf[[i]],colMeans)))
    }
    
    temp$Var1<- c(-24.5:24.5)*10
    temp$TF <- names(tf)[i]
    return(temp)
  },tf=tf.bw.sbin)
  tf.bw.sbin.m<- Reduce(rbind,tf.bw.sbin.m)
  return(tf.bw.sbin.m)
}

a<- Dnase.TFbin.ms(tf.bw.sbin,scale = TRUE)
b<- Dnase.TFbin.ms(tf.bw.sbin)
ggplot(b)+geom_line(aes(x=Var1,y=value,colour=TF))+
  facet_wrap(~variable)+
  labs(x = "",y = " ",title="The average signal in Dnase annotation") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size = 8)
  )



TFbw.geom_line<- function(tf.bw.sbin.m){
lapply(1:length(tf.bw.sbin.m),function(i,tf.bw.sbin.m){
temp<- as.data.frame(tf.bw.sbin.m[[i]])
temp<- melt(temp)
temp$Var1<- c(-24.5:24.5)*10
ggplot(temp)+geom_line(aes(x=Var1,y=value))+
  facet_wrap(~variable)+
  labs(x = "",y = " ",title=names(tf.bw.sbin.m)[i]) +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size = 8)
  )},tf.bw.sbin.m=tf.bw.sbin.m)
}

tf.bw.sbin.sm<-

tf.bw.cov<- lapply(1:length(tf.bw.files),function(i,tf){
  tf.bw<- import.bw(paste0("/mnt/local-disk1/rsgeno2/MAmotif/ENCODE/Tfbs/SydhTfbsK562/",tf[i]),asRle = TRUE)
  tf.cov<- sapply(dnase.K562.anno.sl,function(x)ScoreMatrixBin(tf.bw,windows = x,bin.num = 1,weight.col = "V5")@.Data*500)
  return(sapply(tf.cov, function(x)sum(x>3)/length(x)))},tf=tf.bw.files)
saveRDS(tf.bw.cov,file = "data/TF_bw_cov.rds")

tf.bw.cov<- as.data.frame(tf.bw.cov)
colnames(tf.bw.cov)<-substring(tf.bw.files,first = 9,last = nchar(tf.bw.files)-7)
tf.bw.cov$region<- row.names(tf.bw.cov)
gd<- melt(tf.bw.cov)
ggplot(gd[1:44,])+geom_bar(aes(x=region,y=value,fill=variable),stat = "identity",position = "dodge")+coord_flip()+
  labs(x = " ",y = " ",title="The TFs distribution in the whole genome",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )


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