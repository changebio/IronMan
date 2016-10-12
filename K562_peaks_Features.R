#here analize genomic feature of K562.pk.dnase which was annotated in K562_peaks_annotation.R

##the relationship between K562.pk.dnase and the expression of their target gene


#================================================
##packages
require(ggplot2)
require(reshape2)

##classical theme
hy.theme<-   theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
                   axis.title.x = element_text( face="bold",size=14,colour = "black"),
                   axis.title.y = element_text(color="black", size=14, face="bold"),
                   legend.title =element_text(face = "bold", size = 14, color = "white"),
                   legend.text = element_text(face = "bold", size = 12),
                   axis.text.x = element_text(face = "bold",size=14),
                   axis.text.y = element_text(face="bold", size=14),
                   strip.text.x = element_text(face = "bold",size = 16)
)
##
k562.pk.dnase<- readRDS("data/K562_annotated_dnase_peaks.rds")
k562.rpkm<- readRDS("data/CAGE_RNA_relationship.rds")
table(k562.pk.dnase$SYMBOL %in% k562.rpkm$Gene.Symbol)
#FALSE  TRUE 
#535 28763 
k562.pk.gexp<- k562.rpkm[k562.pk.dnase$SYMBOL,]
k562.pk.gexp$State<- mcols(k562.pk.dnase)[,13]
k562.pk.gexp$Index<- 1:length(k562.pk.dnase)
k562.pk.gexp<- na.omit(k562.pk.gexp)
ggplot(k562.pk.gexp)+
  geom_jitter(aes(x=State,y=wgEncodeCaltechRnaSeqK562R1x75dAlignsRep1V2.bed),width = 0.8,alpha=0.1)+
  geom_boxplot(aes(x=State,y=wgEncodeCaltechRnaSeqK562R1x75dAlignsRep1V2.bed),alpha=0,size=1)+
  scale_y_continuous(trans = "log2")+
  labs(y="log2(RPKM)",title="The gene expression in different groups")+
  hy.theme

ggplot(melt(k562.pk.gexp[,c(3:6,24)]))+
  geom_jitter(aes(x=State,y=value),width = 0.8,alpha=0.1)+
  geom_boxplot(aes(x=State,y=value),alpha=0,size=1)+
  scale_y_continuous(trans = "log2")+
  labs(y="log2(RPKM)",title="The gene expression in different groups")+
  facet_wrap(~variable)+
  hy.theme

k562.pk.gexp.dep<- k562.pk.gexp[k562.pk.gexp$State!="Control",]
k562.pk.gexp.dep<- rbind(k562.pk.gexp.dep,k562.pk.gexp[k562.pk.gexp$State=="Control",])
k562.pk.gexp.dep<- k562.pk.gexp.dep[!duplicated(k562.pk.gexp.dep$Gene.Symbol),]
ggplot(melt(k562.pk.gexp.dep[,c(3:6,24)]))+
  geom_jitter(aes(x=State,y=value),width = 0.8,alpha=0.1)+
  geom_boxplot(aes(x=State,y=value),alpha=0,size=1)+
  scale_y_continuous(trans = "log2")+
  labs(y="log2(RPKM)",title="The gene expression in different groups(dep)")+
  facet_wrap(~variable)+
  hy.theme

###choose K562-distal only H3K4me3 region
k562.ds.me3.exp<- k562.pk.gexp.dep[k562.pk.gexp.dep$State=="None",]
k562.ds.me3.exp.s<- as.data.frame(sapply(3:6, function(i)summary(k562.ds.me3.exp[,i])))
high.exp.idx<- rowSums(as.matrix(k562.ds.me3.exp[,3:6]>t(replicate(318,as.numeric(k562.ds.me3.exp.s[3,])))))
k562.ds.me3.exp.high<- k562.ds.me3.exp[high.exp.idx==4,]

k562.bp.me3<- readBroadPeak("/mnt/local-disk1/rsgeno2/MAmotif/macs2/broad_test_peaks.broadPeak")[1:20000]
k562.np.me3<- readNarrowPeak("/mnt/local-disk1/rsgeno2/MAmotif/macs2/narrow_test_peaks.narrowPeak")[1:20000]
k562.ds.me3.only<- k562.pk.dnase[k562.ds.me3.exp.high$Index]
k562.ds.me3.only$BP<- countOverlaps(k562.ds.me3.only,k562.bp.me3)
k562.ds.me3.only$NP<- countOverlaps(k562.ds.me3.only,k562.np.me3)
saveRDS(k562.ds.me3.only,file = "data/K562_only_H3K4me3_peaks.rds")


#Cage signal in K562 annotated peaks-----------
cage.seq<- readRDS("/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/data/gro.cage.all.seq.rds")
k562.pk.dnase.1k<- resize(k562.pk.dnase,width = 1000,fix="center")
k562.pk.dnase.1k<- k562.pk.dnase.1k[countOverlaps(k562.pk.dnase.1k,k562.pk.dnase.1k)==1]
k562.pk.dnase.1k<- keepSeqlevels(k562.pk.dnase.1k,seqlevels(cage.seq$K562_cell_rep1.minus))
k562.1k.bs<- region.base.signal(k562.pk.dnase.1k,cage.seq,strand = FALSE,weight.col = "V5")
#saveRDS(k562.1k.bs,file = "/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/data/K562_annotated_peaks_1kp_cage_signal.rds")

##The distribution of CAGE reads in K562----------
k562.5h.bs<- readRDS("/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/data/K562_annotated_peaks_1kp_cage_signal.rds")
k562.5h.sum<- lapply(k562.5h.bs,rowSums)
k562.5h.sum<- lapply(seq(1,length(k562.5h.sum),by = 2),function(i,y)return(y[[i]]+y[[i+1]]),y=k562.5h.sum)
names(k562.5h.sum)<- names(k562.5h.bs)[seq(1,length(k562.5h.bs),by=2)]
names(k562.5h.sum)<- substring(names(k562.5h.sum),first = 1,last = nchar(names(k562.5h.sum))-6)

k562.5h.gd<- as.data.frame(sapply(k562.5h.sum,function(x)return(x)))
colnames(k562.5h.gd)<- names(k562.5h.sum)
k562.5h.gd[k562.5h.gd>30]<-30
k562.5h.gd$State<- k562.pk.dnase.5h$State
k562.5h.gd$State<- factor(k562.5h.gd$State,levels = c("None","H3K4me1","H3K27ac","Both","Control","Enhancer","poiProm","actProm"))
k562.5h.gd$Type<- "H3K4me3"
k562.5h.gd$Type[k562.5h.gd$State=="Control"]<- "Other"
k562.5h.gd$Type[k562.5h.gd$State=="Enhancer"]<- "Other"
k562.5h.gd$Type[k562.5h.gd$State=="poiProm"]<- "Promoter"
k562.5h.gd$Type[k562.5h.gd$State=="actProm"]<- "Promoter"

ggplot(melt(k562.5h.gd[,c(12:19)]))+geom_histogram(aes(x=value,fill=Type))+
  facet_grid(variable ~ State)+
  labs(x="read",title="The distribution of CAGE signal",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size = 16)
  )

ggplot(melt(k562.5h.gd[,c(13:19)]))+geom_histogram(aes(x=value,y=..density..,fill=Type))+
  facet_grid(variable ~ State)+
  labs(x="read",y="density",title="The distribution of CAGE signal",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size = 16)
  )

