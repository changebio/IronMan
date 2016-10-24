#here analize genomic feature of K562.pk.dnase which was annotated in K562_peaks_annotation.R

##the relationship between K562.pk.dnase and the expression of their target gene


#================================================

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
k562.bp.me3.sub<- subsetByOverlaps(k562.bp.me3,k562.ds.me3.only)

#Cage signal in K562 annotated peaks-----------
cage.seq<- readRDS("/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/data/gro.cage.all.seq.rds")
cage.seq.dp<- sapply(cage.seq,function(x)sum(x$V5))
#saveRDS(cage.seq.dp,file = "data/Cage_seq_deep.rds")
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

#Dscore
Dsuper<- function(base.sg,width){
  cage.rd.sum<- lapply(base.sg,function(x)rowSums(x[,(ncol(x)/2-width/2+1):(ncol(x)/2+width/2)]))
  cage.rd.sum<- lapply(seq(1,length(cage.rd.sum),by = 2),function(i,y)return((y[[i+1]]-y[[i]])/(y[[i]]+y[[i+1]])),y=cage.rd.sum)
  names(cage.rd.sum)<- names(base.sg)[seq(1,length(base.sg),by=2)]
  names(cage.rd.sum)<- substring(names(cage.rd.sum),first = 1,last = nchar(names(cage.rd.sum))-6)
  
  cage.rd.gd<- as.data.frame(sapply(cage.rd.sum,function(x)return(x)))
  colnames(cage.rd.gd)<- names(cage.rd.sum)
  cage.rd.gd$State<- k562.pk.dnase.5h$State
  cage.rd.gd$State<- factor(cage.rd.gd$State,levels = c("None","H3K4me1","H3K27ac","Both","Control","Enhancer","poiProm","actProm"))
  cage.rd.gd$Type<- "H3K4me3"
  cage.rd.gd$Type[cage.rd.gd$State=="Control"]<- "Other"
  cage.rd.gd$Type[cage.rd.gd$State=="Enhancer"]<- "Other"
  cage.rd.gd$Type[cage.rd.gd$State=="poiProm"]<- "Promoter"
  cage.rd.gd$Type[cage.rd.gd$State=="actProm"]<- "Promoter"
  return(cage.rd.gd)
}
cage.rd.gd<- Dsuper(k562.5h.bs,width = 100)
ggplot(melt(cage.rd.gd[,c(1,17:19)]))+geom_histogram(aes(x=value,y = ..density..,fill=Type))+
  facet_grid(variable ~ State)+
  labs(x="read",title="The distribution of D score(100bp)",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size = 16)
  )

Dscore<- function(base.sg,width){
  cage.rd.sum<- lapply(1:length(base.sg),function(i,y){x=y[[i]];if(i%%2){
    return(rowSums(x[,(ncol(x)/2-width/2+1):(ncol(x)/2)]))
    }else{
        return(rowSums(x[,(ncol(x)/2):(ncol(x)/2+width/2)]))}},y=base.sg)
  cage.rd.sum<- lapply(seq(1,length(cage.rd.sum),by = 2),function(i,y)return((y[[i+1]]-y[[i]])/(y[[i]]+y[[i+1]])),y=cage.rd.sum)
  names(cage.rd.sum)<- names(base.sg)[seq(1,length(base.sg),by=2)]
  names(cage.rd.sum)<- substring(names(cage.rd.sum),first = 1,last = nchar(names(cage.rd.sum))-6)
  
  cage.rd.gd<- as.data.frame(sapply(cage.rd.sum,function(x)return(x)))
  colnames(cage.rd.gd)<- names(cage.rd.sum)
  cage.rd.gd$State<- k562.pk.dnase.5h$State
  cage.rd.gd$State<- factor(cage.rd.gd$State,levels = c("None","H3K4me1","H3K27ac","Both","Control","Enhancer","poiProm","actProm"))
  cage.rd.gd$Type<- "H3K4me3"
  cage.rd.gd$Type[cage.rd.gd$State=="Control"]<- "Other"
  cage.rd.gd$Type[cage.rd.gd$State=="Enhancer"]<- "Other"
  cage.rd.gd$Type[cage.rd.gd$State=="poiProm"]<- "Promoter"
  cage.rd.gd$Type[cage.rd.gd$State=="actProm"]<- "Promoter"
  return(cage.rd.gd)
}
cage.sc.gd<- Dscore(k562.5h.bs,width = 400)
ggplot(melt(cage.sc.gd[,c(1,17:19)]))+geom_histogram(aes(x=value,y = ..density..,fill=Type))+
  facet_grid(variable ~ State)+
  labs(x="read",title="The distribution of normal D score(400bp)",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size = 16)
  )
##average signal

cage.ave.sg<- lapply(1:length(k562.5h.bs),function(i)sapply(split(as.data.frame(k562.5h.bs[[i]]@.Data),k562.pk.dnase.1k$State),function(x)colMeans(x)*10^9/cage.seq.dp[i]))
k562.cage.sg<- as.data.frame(rbind(cage.ave.sg$hg19.ctss_all_plus,cage.ave.sg$hg19.ctss_all_minus))
k562.cage.sg$strand<- c(rep("+",1000),rep("-",1000))
k562.cage.sg$name<- -500:499

ggplot(melt(k562.cage.sg,id.var=c("strand","name")))+geom_line(aes(x=name,y=value,color=strand))+
  facet_grid(variable~.,scales = "free")+
  labs(x = "",y = " ",title="The average of CAGE signal in 1kp regions",linetype="read") +
  hy.theme

##trimed
cage.ave.tm<- lapply(1:length(k562.5h.bs),function(i)sapply(split(as.data.frame(k562.5h.bs[[i]]@.Data),k562.pk.dnase.1k$State),function(x)apply(x,2,mean,trim=.05)*10^9/cage.seq.dp[i]))
names(cage.ave.tm)<- names(k562.5h.bs)
k562.cage.tm<- as.data.frame(rbind(cage.ave.tm$hg19.ctss_all_plus,cage.ave.tm$hg19.ctss_all_minus))
k562.cage.tm$strand<- c(rep("+",1000),rep("-",1000))
k562.cage.tm$name<- -500:499

ggplot(melt(k562.cage.tm,id.var=c("strand","name")))+geom_line(aes(x=name,y=value,color=strand))+
  stat_smooth(method = "auto",level = 0.95)+
  facet_grid(variable~.,scales = "free")+
  labs(x = "",y = " ",title="The average of CAGE signal in 1kp regions",linetype="read") +
  hy.theme


##H3K4me3 
k562.me3.f<- list.files("/mnt/local-disk1/rsgeno2/MAmotif/3.Histone_Broad_hg19/H3K4me3/K562/")[2:3]
k562.me3.bed<- lapply(k562.me3.f, function(x)import.bed(paste0("/mnt/local-disk1/rsgeno2/MAmotif/3.Histone_Broad_hg19/H3K4me3/K562/",x)))
k562.me3.bed<- lapply(k562.me3.bed,function(x){seqlengths(x)<- seqlengths(Hsapiens)[as.character(seqlevels(x))];return(x)})
names(k562.me3.bed)<- k562.me3.f
k562.me3.cov<- lapply(k562.me3.bed,function(x)coverage(x,ifelse(strand(x)=="+",115,-115)))
me3.seq.dp<- lapply(k562.me3.cov,function(x)sum(sum(x)))
k562.me3.bs<- region.base.signal(k562.pk.dnase.1k,k562.me3.cov,strand=FALSE)
me3.seq.dp<- lapply(k562.me3.bs,sum)
me3.ave.sg<- lapply(1:length(k562.me3.bs),function(i)sapply(split(as.data.frame(k562.me3.bs[[i]]@.Data),k562.pk.dnase.1k$State),function(x)colMeans(x)*10^9/me3.seq.dp[[i]]))
names(me3.ave.sg)<- k562.me3.f
k562.me3.sg<- as.data.frame(rbind(me3.ave.sg$wgEncodeBroadHistoneK562H3k4me3StdAlnRep1.bed,me3.ave.sg$wgEncodeBroadHistoneK562H3k4me3StdAlnRep2.bed))
k562.me3.sg$strand<- c(rep("Rep1",1000),rep("Rep2",1000))
k562.me3.sg$name<- -500:499
saveRDS(k562.me3.sg,file = "data/K562_me3_sg.rds")
ggplot(melt(k562.me3.sg,id.var=c("strand","name")))+geom_line(aes(x=name,y=value,color=strand))+
  facet_grid(variable~.,scales = "free")+
  labs(x = "",y = " ",title="The average of H3K4me3 signal in 1kp regions",linetype="read") +
  hy.theme

#DNase
k562.ase.f<- list.files("/mnt/local-disk1/rsgeno2/MAmotif/ENCODE/2.DNase_Duke_hg19/",pattern = "bed")[3:4]
k562.ase.bed<- lapply(k562.ase.f, function(x)import.bed(paste0("/mnt/local-disk1/rsgeno2/MAmotif/ENCODE/2.DNase_Duke_hg19/",x)))
k562.ase.bed<- lapply(k562.ase.bed,function(x){seqlengths(x)<- seqlengths(Hsapiens)[as.character(seqlevels(x))];return(x)})
names(k562.ase.bed)<- k562.ase.f
k562.ase.cov<- lapply(k562.ase.bed,coverage)
k562.dnase.bw<- import.bw("/mnt/local-disk1/rsgeno2/MAmotif/ENCODE/2.DNase_Duke_hg19/wgEncodeOpenChromDnaseK562BaseOverlapSignalV2.bigWig",asRle=TRUE)
k562.dnase.sig<- import.bw("/mnt/local-disk1/rsgeno2/MAmotif/ENCODE/2.DNase_Duke_hg19/wgEncodeOpenChromDnaseK562SigV2.bigWig",asRle=TRUE)
k562.ase.bs<- region.base.signal(k562.pk.dnase.1k,c(k562.ase.cov,list(K562Sig=k562.dnase.sig,K562BaseSig=k562.dnase.bw)),strand=FALSE)
saveRDS(k562.ase.bs,file = "data/K562_Dnase_base_signal.rds")

ase.seq.dp<- lapply(k562.ase.bs,sum)
ase.ave.sg<- lapply(1:length(k562.ase.bs),function(i)sapply(split(as.data.frame(k562.ase.bs[[i]]@.Data),k562.pk.dnase.1k$State),function(x)colMeans(x)*10^9/ase.seq.dp[[i]]))
names(ase.ave.sg)<- k562.ase.f
k562.ase.sg<- as.data.frame(rbind(ase.ave.sg$wgEncodeOpenChromDnaseK562AlnRep1.bed,ase.ave.sg$wgEncodeOpenChromDnaseK562AlnRep2.bed))
k562.ase.sg$strand<- c(rep("Rep1",1000),rep("Rep2",1000))
k562.ase.sg$name<- -500:499
saveRDS(k562.ase.sg,file = "data/K562_ase_sg.rds")
ggplot(melt(k562.ase.sg,id.var=c("strand","name")))+geom_line(aes(x=name,y=value,color=strand))+
  facet_grid(variable~.,scales = "free")+
  labs(x = "",y = " ",title="The average of DNase signal in 1kp regions") +
  hy.theme

#summit position
ase.smt.pos<- lapply(1:length(k562.ase.bs),function(i)sapply(split(as.data.frame(k562.ase.bs[[i]]@.Data[,425:575]),k562.pk.dnase.1k$State),function(x)apply(x,1,which.max)))
names(ase.smt.pos)<- names(k562.ase.bs)


#Pol2
k562.pol.f<- list.files("/mnt/local-disk1/rsgeno2/MAmotif/ENCODE/Tfbs/",pattern = "Pol2")
k562.pol.bed<- lapply(k562.pol.f, function(x)import.bw(paste0("/mnt/local-disk1/rsgeno2/MAmotif/ENCODE/Tfbs/",x),asRle=TRUE))
k562.pol.bed<- lapply(k562.pol.bed,function(x){seqlengths(x)<- seqlengths(Hsapiens)[as.character(seqlevels(x))];return(x)})
names(k562.pol.bed)<- k562.pol.f
k562.pol.cov<- lapply(k562.pol.bed,function(x)coverage(x,ifelse(strand(x)=="+",100,-100)))

k562.pol.bs<- region.base.signal(k562.pk.dnase.1k,k562.pol.bed,strand=FALSE)
pol.seq.dp<- lapply(k562.pol.bs,sum)
pol.ave.sg<- lapply(1:length(k562.pol.bs),function(i)sapply(split(as.data.frame(k562.pol.bs[[i]]@.Data),k562.pk.dnase.1k$State),function(x)colMeans(x)*10^9/pol.seq.dp[[i]]))
names(pol.ave.sg)<- k562.pol.f
k562.pol.sg<- do.call(rbind.data.frame,pol.ave.sg)
k562.pol.sg$strand<- rep(names(pol.ave.sg),each=1000)
k562.pol.sg$name<- -500:499
saveRDS(k562.pol.sg,file = "data/K562_pol_sg.rds")
ggplot(melt(k562.pol.sg,id.var=c("strand","name")))+geom_line(aes(x=name,y=value,color=strand))+
  facet_grid(variable~strand,scales = "free")+
  labs(x = "",y = " ",title="The average of Pol2 signal in 1kp regions") +
  hy.theme



#K562 DNA Methylation
sl725<- import("/mnt/local-disk1/rsgeno2/MAmotif/ENCODE/Methyl/SL725.DCC.CGs.bed")
mcols(sl725)<- as.numeric(sl725$blockSizes)
sl726<- import("/mnt/local-disk1/rsgeno2/MAmotif/ENCODE/Methyl/SL726.DCC.CGs.bed")
mcols(sl726)<- as.numeric(sl726$blockSizes)
k562.dna.met<- list(SL725=sl725,SL726=sl726)
k562.dna.met<- lapply(k562.dna.met, function(x){seqlengths(x)<-seqlengths(Hsapiens)[as.character(seqlevels(x))];return(x)})
k562.met.bs<- region.base.signal(k562.pk.dnase.1k,k562.dna.met,strand=FALSE,weight.col = "X")
met.ave.sg<- lapply(1:length(k562.met.bs),function(i)sapply(split(as.data.frame(k562.met.bs[[i]]@.Data),k562.pk.dnase.1k$State),function(x)colMeans(x)))
names(met.ave.sg)<- names(k562.met.bs)
met.ave.sg<- lapply(met.ave.sg,function(x)t(sapply(seq(1,1000,50),function(i)colMeans(x[i:(i+49),]))))
k562.met.sg<- do.call(rbind.data.frame,met.ave.sg)
k562.met.sg$strand<- rep(names(met.ave.sg),each=20)
k562.met.sg$name<- seq(-500,499,50)
saveRDS(k562.met.sg,file = "data/K562_met_sg.rds")
ggplot(melt(k562.met.sg,id.var=c("strand","name")))+geom_line(aes(x=name,y=value,color=strand))+
  facet_grid(variable~.,scales = "free")+
  labs(x = "",y = " ",title="The average of methylation CpG signal in 1kp regions") +
  hy.theme
