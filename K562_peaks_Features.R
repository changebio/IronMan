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
