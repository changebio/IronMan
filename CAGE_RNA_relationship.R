k562.rpkm.rep1<- read.table("/mnt/local-disk1/rsgeno2/MAmotif/2.Processing/5.RNAseq_hg19_DEseq/rpkms/K562_1x75_rpkms.txt",header = TRUE)
k562.rpkm.rep1<- k562.rpkm.rep1[!duplicated(k562.rpkm.rep1$Gene.Symbol),]
rownames(k562.rpkm.rep1)<- k562.rpkm.rep1$Gene.Symbol
k562.rpkm.rep1<- k562.rpkm.rep1[hmg.pro$symbol,]

k562.rpkm.rep2<- read.table("/mnt/local-disk1/rsgeno2/MAmotif/2.Processing/5.RNAseq_hg19_DEseq/rpkms/K562_2x75_rpkms.txt",header = TRUE)
k562.rpkm.rep2<- k562.rpkm.rep2[!duplicated(k562.rpkm.rep2$Gene.Symbol),]
rownames(k562.rpkm.rep2)<- k562.rpkm.rep2$Gene.Symbol
k562.rpkm.rep2<- k562.rpkm.rep2[hmg.pro$symbol,]

cor(k562.rpkm.rep1$wgEncodeCaltechRnaSeqK562R1x75dAlignsRep1V2.bed,k562.rpkm.rep1$wgEncodeCaltechRnaSeqK562R1x75dAlignsRep2V2.bed)
plot(k562.rpkm.rep1$wgEncodeCaltechRnaSeqK562R1x75dAlignsRep1V2.bed,k562.rpkm.rep1$wgEncodeCaltechRnaSeqK562R1x75dAlignsRep2V2.bed)

require(ggplot2)
require(reshape2)
ggplot(k562.rpkm.rep1)+geom_point(aes(x=wgEncodeCaltechRnaSeqK562R1x75dAlignsRep1V2.bed,y=wgEncodeCaltechRnaSeqK562R1x75dAlignsRep2V2.bed),alpha = 1/10)+
  scale_x_continuous(trans = "log2")+
  scale_y_continuous(trans = "log2")
  
require(GenomicAlignments)

hmg<- genes(txdb,columns="gene_id")
gene.df <- bitr(names(hmg), fromType = "ENTREZID", toType = c("SYMBOL"), annoDb = "org.Hs.eg.db")

hmg<- hmg[gene.df$ENTREZID]
hmg$symbol<- gene.df$SYMBOL
hmg.pro<- promoters(hmg,upstream = 250,downstream = 250)
hmg.pro<- keepSeqlevels(hmg.pro,seqlevels(cage.seq$K562_cell_rep1.minus))
hmg.pro.bs<- region.base.signal(hmg.pro,gro.cage.seq,strand = FALSE,weight.col = "V5")
#saveRDS(hmg.pro.bs,file = "/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/data/human_Promoter_cage_bsignl.rds")
hmg.pro.bs<- readRDS("/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/data/human_Promoter_cage_bsignl.rds")

hmg.pro.sum<- lapply(hmg.pro.bs,rowSums)
hmg.pro.sum<- lapply(seq(1,length(hmg.pro.sum),by = 2),function(i,y)return(y[[i]]+y[[i+1]]),y=hmg.pro.sum)
names(hmg.pro.sum)<- names(hmg.pro.bs)[seq(1,length(hmg.pro.bs),by=2)]
names(hmg.pro.sum)<- substring(names(hmg.pro.sum),first = 1,last = nchar(names(hmg.pro.sum))-6)

hmg.pro.rpkm<- as.data.frame(sapply(hmg.pro.sum,function(x)return(x*10^9/500/sum(x))))
colnames(hmg.pro.rpkm)<- names(hmg.pro.sum)
ggplot(hmg.pro.rpkm)+geom_point(aes(x=K562_cell_rep1,y=K562_cell_rep2),alpha = 1/10)+
  scale_x_continuous(trans = "log2")+
  scale_y_continuous(trans = "log2")

temp<- cbind(k562.rpkm.rep1,k562.rpkm.rep2[,3:4],hmg.pro.rpkm)
temp<- na.omit(temp)
#saveRDS(temp,file = "data/CAGE_RNA_relationship.rds")
temp<- readRDS("data/CAGE_RNA_relationship.rds")
##heatmap of the association between CAGE and RNA
#heatmap(cor(temp[,3:23]),labCol = FALSE,main="The association between CAGE signal and RNA signal in genes")
data <- cor(temp[,3:23])
ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order
ord
pd <- as.data.frame( data[ord,ord] )
pd$Sample<- rownames
pd.m<- melt(pd)
pd.m$Sample<- factor(pd.m$Sample,levels = rownames(pd))
pd.m$variable<- factor(pd.m$variable,levels = rownames(pd))
ggplot(pd.m, aes(Sample,variable) ) +
  geom_tile(aes(fill = value))+
  labs(x="",y="",title="The association between RNA and CAGE") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "white"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_blank(),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size = 16)
  )


###RNA production as measured by CAGE tag density at TSSs in K562 cells-----------
require(ggplot2)
require(reshape2)

temp.gd<-melt(temp,id.vars = c("Transcipt.ID","Gene.Symbol","wgEncodeCaltechRnaSeqK562R1x75dAlignsRep1V2.bed"))
ggplot(temp.gd)+geom_point(aes(x=wgEncodeCaltechRnaSeqK562R1x75dAlignsRep1V2.bed,y=value),alpha = 1/10)+
  scale_x_continuous(trans = "log2")+
  scale_y_continuous(trans = "log2")+
  facet_wrap( ~ variable)+
  labs(title="The association between RNA and CAGE") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size = 16)
  )


#region.base.signal---------
region.base.signal<- function(region,cage.seq,strand=TRUE,weight.col=NULL,...){
  if(strand){
    cage.sg.p<-lapply(cage.seq,function(x)ScoreMatrix(x,region[strand(region)=="+"],weight.col = weight.col))
    cage.sg.m<-lapply(cage.seq,function(x)ScoreMatrix(x,region[strand(region)=="-"],weight.col = weight.col))
    cage.sg<- list(plus=cage.sg.p,minus=cage.sg.m)
    
  }else{
    cage.sg<-lapply(cage.seq,function(x)ScoreMatrix(x,region,weight.col = weight.col))
    
  }
  return(cage.sg)
}


