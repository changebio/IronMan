k562.rpkm.rep1<- read.table("/mnt/local-disk1/rsgeno2/MAmotif/2.Processing/5.RNAseq_hg19_DEseq/rpkms/K562_1x75_rpkms.txt",header = TRUE)
k562.rpkm.rep1<- k562.rpkm.rep1[!duplicated(k562.rpkm.rep1$Gene.Symbol),]
rownames(k562.rpkm.rep1)<- k562.rpkm.rep1$Gene.Symbol
k562.rpkm.rep1<- k562.rpkm.rep1[hmg.pro$symbol,]

k562.rpkm.rep2<- read.table("/mnt/local-disk1/rsgeno2/MAmotif/2.Processing/5.RNAseq_hg19_DEseq/rpkms/K562_2x75_rpkms.txt",header = TRUE)
cor(k562.rpkm.rep1$wgEncodeCaltechRnaSeqK562R1x75dAlignsRep1V2.bed,k562.rpkm.rep1$wgEncodeCaltechRnaSeqK562R1x75dAlignsRep2V2.bed)
plot(k562.rpkm.rep1$wgEncodeCaltechRnaSeqK562R1x75dAlignsRep1V2.bed,k562.rpkm.rep1$wgEncodeCaltechRnaSeqK562R1x75dAlignsRep2V2.bed)

require(ggplot2)
require(reshape2)
ggplot(k562.rpkm.rep1)+geom_point(aes(x=wgEncodeCaltechRnaSeqK562R1x75dAlignsRep1V2.bed,y=wgEncodeCaltechRnaSeqK562R1x75dAlignsRep2V2.bed),alpha = 1/10)+
  scale_x_continuous(trans = "log2")+
  scale_y_continuous(trans = "log2")
  
require(GenomicAlignments)

hmg<- genes(txdb,columns="gene_id")
hmg<- hmg[gene.df$ENTREZID]
hmg$symbol<- gene.df$SYMBOL
hmg.pro<- promoters(hmg,upstream = 250,downstream = 250)
hmg.pro<- keepSeqlevels(hmg.pro,seqlevels(cage.seq$K562_cell_rep1.minus))
hmg.pro.bs<- region.base.signal(hmg.pro,cage.seq,strand = FALSE,weight.col = "V5")
hmg.pro.sum<- lapply(hmg.pro.bs,rowSums)
hmg.pro.sum<- lapply(seq(1,length(hmg.pro.sum),by = 2),function(i,y)return(y[[i]]+y[[i+1]]),y=hmg.pro.sum)
names(hmg.pro.sum)<- names(cage.seq)[seq(1,length(cage.seq),by=2)]
names(hmg.pro.sum)<- substring(names(hmg.pro.sum),first = 1,last = nchar(names(hmg.pro.sum))-6)
cage.sum<- sapply(cage.seq,function(x)sum(x$V5))
cage.sum<- sapply(seq(1,length(cage.sum),by = 2),function(i,y)return(y[[i]]+y[[i+1]]),y=cage.sum)
hmg.pro.rpkm<- as.data.frame(sapply(1:13,function(i,x,y)return(x[[i]]*10^9/500/y[i]),x=hmg.pro.sum,y=cage.sum))
colnames(hmg.pro.rpkm)<- names(hmg.pro.sum)
ggplot(hmg.pro.rpkm)+geom_point(aes(x=K562_cell_rep1,y=K562_cell_rep2),alpha = 1/10)+
  scale_x_continuous(trans = "log2")+
  scale_y_continuous(trans = "log2")


ind<-seqnames(hmg.pro)!="chrY" & !is.na(k562.rpkm.rep1$Gene.Symbol)

temp<- cbind(k562.rpkm.rep1,hmg.pro.rpkm)
temp<- na.omit(temp)
ggplot(temp)+geom_point(aes(x=wgEncodeCaltechRnaSeqK562R1x75dAlignsRep1V2.bed,y=K562_cell_rep2),alpha = 1/10)+
  scale_x_continuous(trans = "log2")+
  scale_y_continuous(trans = "log2")


gene.df <- bitr(names(hmg), fromType = "ENTREZID", toType = c("SYMBOL"), annoDb = "org.Hs.eg.db")
gene.df


