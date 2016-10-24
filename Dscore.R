#################################################################
# Direction by CAGE and GRO 
# There are two parts
# First: average signal
# Second: Dscore
#################################################################


#load packages-----
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
require(GenomicFeatures)
require(genomation)

# Part I: average signal
##functions
###---------------------------------
region.cage.signal<- function(region,cage.seq,strand=TRUE,bin.num = 50,weight.col=NULL){
  if(strand){
    cage.sg.p<-lapply(cage.seq,function(x)ScoreMatrixBin(x,region[strand(region)=="+"],bin.num = bin.num,weight.col = weight.col))
    cage.sg.m<-lapply(cage.seq,function(x)ScoreMatrixBin(x,region[strand(region)=="-"],bin.num = bin.num,weight.col = weight.col))
    cage.sg<- melt(rbind(sapply(cage.sg.p,colMeans),sapply(cage.sg.m,colMeans)))
    cage.sg$gene<- factor(c(rep("+",bin.num),rep("-",bin.num)),levels = c("+","-"))
    cage.sg$read<- factor(c(rep("minus",2*bin.num),rep("plus",2*bin.num)),levels = c("plus","minus"))
    cage.sg$type <- substring(cage.sg$Var2,first = 1,last = nchar(as.character(cage.sg$Var2))-c(rep(6,2*bin.num),rep(5,2*bin.num)))
    cage.sg$Var1<- c(-24.5:24.5)*10
  }else{
    cage.sg<-lapply(cage.seq,function(x)ScoreMatrixBin(x,region,bin.num = bin.num,weight.col = weight.col))
    cage.sg<- melt(sapply(cage.sg,colMeans))
    cage.sg$gene<- factor(rep("*",bin.num),"*")
    cage.sg$read<- factor(c(rep("minus",bin.num),rep("plus",bin.num)),levels = c("plus","minus"))
    cage.sg$type <- substring(cage.sg$Var2,first = 1,last = nchar(as.character(cage.sg$Var2))-c(rep(6,bin.num),rep(5,bin.num)))
    cage.sg$Var1<- c(-24.5:24.5)*10
  }
  return(cage.sg)
}


region.cage.signal.geom_line<- function(cage.sg,title="The average signal of CAGE",...){
  require(ggplot2)
  p<- ggplot(cage.sg)
  if(length(levels(cage.sg$gene))==1){
    p<- p + geom_line(aes(x=Var1,y=value,linetype= read))
  }else{
    p<- p + geom_line(data=cage.sg[cage.sg$gene=="+",],aes(x=Var1,y=value,colour=gene,linetype=read)) + geom_line(data=cage.sg[cage.sg$gene=="-",],aes(x=Var1,y=-value,colour=gene,linetype=read))
  }
  p<- p + scale_color_manual(values=c("#0099FF","#FF6600"))+
    facet_wrap( ~ type)+
    labs(x = "",y = " ",title=title) +
    theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
          axis.title.x = element_text( face="bold",size=14),
          axis.title.y = element_text(color="black", size=14, face="bold"),
          legend.title =element_text(face = "bold", size = 14, color = "black"),
          legend.text = element_text(face = "bold", size = 12),
          axis.text.x = element_text(face="bold",size=14),
          axis.text.y = element_text(face="bold", size=14),
          strip.text.x = element_text(face = "bold",size = 16)
    )
  p
}

###================================



#gene signal calculation
##gro and cage signal distribution in promoters
ucsc.gene<- genes(txdb)
ucsc.gene.pt<- promoters(ucsc.gene,upstream = 250,downstream = 250)
gro.sg.p<-lapply(gro.seq,function(x)ScoreMatrixBin(x,ucsc.gene.pt[strand(ucsc.gene.pt)=="+"],bin.num = 50,weight.col = "V5"))
gro.sg.m<-lapply(gro.seq,function(x)ScoreMatrixBin(x,ucsc.gene.pt[strand(ucsc.gene.pt)=="-"],bin.num = 50,weight.col = "V5"))

gro.sg<- melt(rbind(sapply(gro.sg.p,colMeans),sapply(gro.sg.m,colMeans)))
gro.sg$gene<- c(rep("+",50),rep("-",50))
gro.sg$read<- factor(c(rep("minus",100),rep("plus",100)),levels = c("plus","minus"))
gro.sg$type <- c(rep("wTAP",200),rep("nTAP",200),rep("CAGE",200))
gro.sg$Var1<- c(-24.5:24.5)*10


ggplot(gro.sg)+geom_line(data=gro.sg[gro.sg$gene=="+",],aes(x=Var1,y=value,colour=gene,linetype=as.factor(read)))+
  geom_line(data=gro.sg[gro.sg$gene=="-",],aes(x=Var1,y=-value,colour=gene,linetype=as.factor(read)))+
  scale_color_manual(values=c("#0099FF","#FF6600"))+
  facet_grid(. ~ type)+
  labs(x = "",y = " ",title="The average signal of Promoters with 500bp width in K562",linetype="read") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size = 16)
  )

###cage from ENCODE
cage.seq<- readRDS("/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/data/cage_gr.rds")
cage.sg.p<-lapply(cage.seq,function(x)ScoreMatrixBin(x,ucsc.gene.pt[strand(ucsc.gene.pt)=="+"],bin.num = 50,weight.col = "V5"))
cage.sg.m<-lapply(cage.seq,function(x)ScoreMatrixBin(x,ucsc.gene.pt[strand(ucsc.gene.pt)=="-"],bin.num = 50,weight.col = "V5"))

cage.sg<- melt(rbind(sapply(cage.sg.p,colMeans),sapply(cage.sg.m,colMeans)))
cage.sg$gene<- c(rep("+",50),rep("-",50))
cage.sg$read<- factor(c(rep("minus",100),rep("plus",100)),levels = c("plus","minus"))
cage.sg$type <- substring(cage.sg$Var2,first = 1,last = nchar(as.character(cage.sg$Var2))-c(rep(6,100),rep(5,100)))
#cage.sg$region<- c(rep("cell",400),rep("chromatin",200),rep("cytosol",800),rep("nucleolus",200),rep("nucleoplasm",200),rep("nucleus",600),rep("polysome",200))
cage.sg$Var1<- c(-24.5:24.5)*10


gro.cage.sg<- rbind(gro.sg,cage.sg)
ggplot(gro.cage.sg)+geom_line(data=gro.cage.sg[gro.cage.sg$gene=="+",],aes(x=Var1,y=value,colour=gene,linetype=as.factor(read)))+
  geom_line(data=gro.cage.sg[gro.cage.sg$gene=="-",],aes(x=Var1,y=-value,colour=gene,linetype=as.factor(read)))+
  scale_color_manual(values=c("#0099FF","#FF6600"))+
  facet_wrap( ~ type)+
  labs(x = "",y = " ",title="The average signal of Promoters with 500bp width in K562",linetype="read") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size = 16)
  )

### the average signal of CAGE bigWig from ENCODE
cage.bw<- readRDS("/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/data/cage_bw.rds")

### the average signal for all merged cage
ucsc.gene<- genes(txdb)
ucsc.gene.pt<- promoters(ucsc.gene,upstream = 250,downstream = 250)
all.gr.m<- readRDS("/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/data/hg19.ctss_all_minus_pool.rds")
all.gr.p<- readRDS("/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/data/hg19.ctss_all_plus_pool.rds")

all.sg.p<-lapply(list(all.gr.m,all.gr.p),function(x)ScoreMatrixBin(x,ucsc.gene.pt[strand(ucsc.gene.pt)=="+"],bin.num = 50,weight.col = "V5"))
all.sg.m<-lapply(list(all.gr.m,all.gr.p),function(x)ScoreMatrixBin(x,ucsc.gene.pt[strand(ucsc.gene.pt)=="-"],bin.num = 50,weight.col = "V5"))

all.sg<- melt(rbind(sapply(all.sg.p,colMeans),sapply(all.sg.m,colMeans)))
all.sg$gene<- c(rep("+",50),rep("-",50))
all.sg$read<- factor(c(rep("minus",100),rep("plus",100)),levels = c("plus","minus"))
all.sg$Var1<- c(-24.5:24.5)*10
saveRDS(all.sg,file = "data/all.sg.rds")
ggplot(all.sg)+geom_line(data=all.sg[all.sg$gene=="+",],aes(x=Var1,y=value,colour=gene,linetype=as.factor(read)))+
  geom_line(data=all.sg[all.sg$gene=="-",],aes(x=Var1,y=-value,colour=gene,linetype=as.factor(read)))+
  scale_color_manual(values=c("#0099FF","#FF6600"))+
  labs(x = "",y = " ",title="The average signal of Promoters with 500bp width in all cells",linetype="read") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size = 16)
  )


##other signal distribution (H3K4me3,H3K4me1,....)

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#part 2 Dscore =========================

##Dscore function
Dscore<- function(plus,minus,windows,wt = NULL,cutoff=0){
  plus<- ScoreMatrixBin(plus,windows,weight.col = wt,bin.num = 2 )
  minus<- ScoreMatrixBin(minus,windows,weight.col = wt,bin.num = 2)
  index<- intersect(rownames(plus),rownames(minus))
  windows <- windows[as.numeric(index)]
  plus<- rowSums(plus[index])*width(windows)
  minus<- rowSums(minus[index])*width(windows)
  Dscore<- (plus-minus)/(plus+minus)
  mcols(windows)<- DataFrame(plus,minus,Dscore)
  return(windows[plus+minus > cutoff])
}

region.cage.Dscore<- function(region,cage.seq,wt = "V5",cutoff=6){
  names(region)<- 1:length(region)
  gc.ds<- lapply(seq(from=1,to=length(cage.seq),by=2),function(i,y)Dscore(y[[i+1]],y[[i]],region,wt = "V5",cutoff = cutoff),y=cage.seq)
  names(gc.ds)<- substring(names(cage.seq),first = 1,last = nchar(names(cage.seq))-c(6,5))[seq(from=1,to=length(cage.seq),by=2)]
  gc.ds.idx<- Reduce(intersect,lapply(gc.ds,names))
  gene.ds<- sapply(gc.ds,function(x)x[gc.ds.idx]$Dscore)
  gene.ds<- as.data.frame(gene.ds)
  gene.ds$strand<- as.character(strand(region[gc.ds.idx]))
  gene.ds.gp<-melt(gene.ds)
  return(gene.ds.gp)
}



###Dscore distribution in region 500bp
ucsc.gene<- genes(txdb)
ucsc.gene.pt<- promoters(ucsc.gene,upstream = 250,downstream = 250)
names(gro.seq)<- substring(names(gro.seq),first = 1,last = nchar(names(gro.seq))-9)
gro.cage.seq<- c(gro.seq,cage.seq)
gc.ds<- lapply(seq(from=1,to=length(gro.cage.seq),by=2),function(i,y)Dscore(y[[i+1]],y[[i]],ucsc.gene.pt,wt = "V5"),y=gro.cage.seq)
names(gc.ds)<- substring(names(gro.cage.seq),first = 1,last = nchar(names(gro.cage.seq))-c(6,5))[seq(from=1,to=length(gro.cage.seq),by=2)]
gc.ds.idx<- Reduce(intersect,lapply(gc.ds,names))
gene.ds<- sapply(gc.ds,function(x)x[gc.ds.idx]$Dscore)
gene.ds<- as.data.frame(gene.ds)
gene.ds$strand<- factor(as.character(strand(ucsc.gene.pt[gc.ds.idx])),levels = c("+","-"))
saveRDS(gene.ds,file = "/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/data/gene.ds.rds")

gene.ds<- readRDS("/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/data/gene.ds.rds")
gene.ds.gp<-melt(gene.ds)

ggplot(gene.ds.gp)+geom_histogram(aes(x=value,y=..count..,fill=strand))+
  facet_wrap( ~ variable)+
  labs(x = "",y = " ",title="Histogram of D score of Promoters with 500bp width in K562") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size=16)
  )


ggplot(gene.ds.gp)+geom_histogram(data=gene.ds.gp[gene.ds.gp$strand=="+",],aes(x=value,y=..count..,fill=strand))+
  geom_histogram(data=gene.ds.gp[gene.ds.gp$strand=="-",],aes(x=value,y=-..count..,fill=strand))+
  scale_fill_manual(values=c("#0099FF","#FF6600"))+
  facet_wrap( ~ variable)+
  labs(x = "",y = " ",title="Histogram of D score of Promoters with 500bp width in K562") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size=16)
  )


### the Dscore of all merged cage
all.cage.ds<- Dscore(all.gr.p,all.gr.m,ucsc.gene.pt,wt = "score")
saveRDS(all.cage.ds,file = "data/all.cage.ds.rds")

all.cage.ds<- readRDS("data/all.cage.ds.rds")
all.cage.ds<- as.data.frame(all.cage.ds)

ggplot(all.cage.ds)+geom_histogram(aes(x=Dscore,y=..count..,fill=strand))+
  labs(x = "",y = " ",title="Histogram of D score of Promoters with 500bp width") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size=16)
  )

ggplot(all.cage.ds)+geom_histogram(data=all.cage.ds[all.cage.ds$strand=="+",],aes(x=Dscore,y=..count..,fill=strand))+
  geom_histogram(data=all.cage.ds[all.cage.ds$strand=="-",],aes(x=Dscore,y=-..count..,fill=strand))+
  scale_fill_manual(values=c("#0099FF","#FF6600"))+
  labs(x = "",y = " ",title="Histogram of D score of Promoters with 500bp width") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size=16)
  )

### the Dscore of wTAP - noTAP
tap.p<- intersectScoreMatrixList(ScoreMatrixList(gro.sg.p[1:4]))
tap.p.rs<- lapply(tap.p,rowSums)
tap.p.rs$plus<- tap.p.rs$GSM1480321_K562_GROcap_wTAP_plus - tap.p.rs$GSM1480322_K562_GROcap_noTAP_plus
tap.p.rs$minus<- tap.p.rs$GSM1480321_K562_GROcap_wTAP_minus - tap.p.rs$GSM1480322_K562_GROcap_noTAP_minus
tap.p.ds<- (tap.p.rs$plus-tap.p.rs$minus)/(tap.p.rs$plus+tap.p.rs$minus)

tap.m<- intersectScoreMatrixList(ScoreMatrixList(gro.sg.m[1:4]))
tap.m.rs <- lapply(tap.m,rowSums)
tap.m.rs$plus<- tap.m.rs$GSM1480321_K562_GROcap_wTAP_plus - tap.m.rs$GSM1480322_K562_GROcap_noTAP_plus
tap.m.rs$minus<- tap.m.rs$GSM1480321_K562_GROcap_wTAP_minus - tap.m.rs$GSM1480322_K562_GROcap_noTAP_minus
tap.m.ds<- (tap.m.rs$plus - tap.m.rs$minus)/(tap.m.rs$plus+tap.m.rs$minus)

tap.ds<- data.frame(Dscore=c(tap.p.ds,tap.m.ds),strand=c(rep("+",length(tap.p.ds)),rep("-",length(tap.m.ds))))



ggplot(tap.ds)+geom_histogram(data=tap.ds[tap.ds$strand=="+",],aes(x=Dscore,y=..count..,fill=strand))+
  geom_histogram(data=tap.ds[tap.ds$strand=="-",],aes(x=Dscore,y=-..count..,fill=strand))+
  xlim(-1,1)+
  scale_fill_manual(values=c("#0099FF","#FF6600"))+
  labs(x = "",y = " ",title="Histogram of D score of Promoters with 500bp width(wTAP - noTAP)") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size=16)
  )



##

grl.k4me3.ma.c <- lapply(k4me3.smt,function(x)resize(x,width = 240,fix="center"))
grl.k4me3.ma.ds<- lapply(grl.k4me3.ma.c,Dscore)
lapply(grl.k4me3.ma.ds,function(x)print(hist(as.numeric(mcols(x)[,2]),main = "Histogram of D score of genes in K562 GROcap wTAP")))
temp<-sapply(grl.k4me3.ma.ds,function(x)as.numeric(mcols(x)[,3]))
temp<- melt(temp)
ggplot(temp,aes(value))+geom_histogram()+facet_grid(.~L1)+
  labs(x = "",y = " ",title="The Dscore of H3K4me3 in K562 CAGE",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )

ggplot(temp,aes(value))+geom_density()+facet_grid(.~L1)

ucsc.prom<- promoters(txdb,upstream = 2000,downstream = 2000,columns=c("tx_id", "tx_name","gene_id"))
ucsc.prom<- keepSeqlevels(ucsc.prom,paste0("chr",c(1:22,"X")))
ucsc.prom<- unique(ucsc.prom)
CTSS.prom <- readRDS("data/CTSS_promoters.rds")
CTSS.prom <- unique(CTSS.prom)
only.prom<- ucsc.prom[countOverlaps(ucsc.prom,CTSS.prom,ignore.strand=TRUE)==1]
only.prom<- only.prom[countOverlaps(only.prom,ucsc.gene,ignore.strand=TRUE)<2]




only.prom<- only.prom[countOverlaps(only.prom,ucsc.gene,ignore.strand=TRUE)<2]
