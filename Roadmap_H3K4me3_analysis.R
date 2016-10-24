
k562.me3.np.f<- list.files("/mnt/local-disk1/rsgeno2/MAmotif/Roadmap/cons_np/narrow/")
k562.me3.np<- lapply(k562.me3.np.f,function(x)readNarrowPeak(paste0("/mnt/local-disk1/rsgeno2/MAmotif/Roadmap/cons_np/narrow/",x))[1:20000])
names(k562.me3.np)<- substr(k562.me3.np.f,start = 1,stop = 4)



k562.me3.bp.f<- list.files("/mnt/local-disk1/rsgeno2/MAmotif/Roadmap/cons_np/broadp/")
k562.me3.bp<- lapply(k562.me3.bp.f,function(x)readBroadPeak(paste0("/mnt/local-disk1/rsgeno2/MAmotif/Roadmap/cons_np/broadp/",x))[1:20000])
names(k562.me3.bp)<- substr(k562.me3.bp.f,start = 1,stop = 4)

k562.me1.np.f<- list.files("/mnt/local-disk1/rsgeno2/MAmotif/Roadmap/cons_np/H3K4me1/narrow/")
k562.me1.np<- lapply(k562.me1.np.f,function(x)readNarrowPeak(paste0("/mnt/local-disk1/rsgeno2/MAmotif/Roadmap/cons_np/H3K4me1/narrow/",x))[1:20000])
names(k562.me1.np)<- substr(k562.me1.np.f,start = 1,stop = 4)

k562.me1.bp.f<- list.files("/mnt/local-disk1/rsgeno2/MAmotif/Roadmap/cons_np/H3K4me1/broadp/")
k562.me1.bp<- lapply(k562.me1.bp.f,function(x)readBroadPeak(paste0("/mnt/local-disk1/rsgeno2/MAmotif/Roadmap/cons_np/H3K4me1/broadp/",x))[1:20000])
names(k562.me1.bp)<- substr(k562.me1.bp.f,start = 1,stop = 4)

k562.ac.np.f<- list.files("/mnt/local-disk1/rsgeno2/MAmotif/Roadmap/cons_np/H3K27ac/narrow/")
k562.ac.np<- lapply(k562.ac.np.f,function(x)readNarrowPeak(paste0("/mnt/local-disk1/rsgeno2/MAmotif/Roadmap/cons_np/H3K27ac/narrow/",x))[1:20000])
names(k562.ac.np)<- substr(k562.ac.np.f,start = 1,stop = 4)

H3K4me3.annotate<- function(x){
  me3.anno<- annotatePeak(k562.me3.np[[x]], tssRegion=c(-2000, 2000), TxDb =txdb, annoDb="org.Hs.eg.db",genomicAnnotationPriority = c("Promoter", "Exon", "Intron","Intergenic"))
  me3.anno<- as.GRanges(me3.anno)
  me3.anno$annotation[grep("Intron",me3.anno$annotation)]<- "Intron"
  me3.anno$annotation[grep("Exon",me3.anno$annotation)]<- "Exon"
  me3.anno$annotation[grep("Downstream",me3.anno$annotation)]<- "Intergenic"
  me3.anno$annotation[me3.anno$annotation=="Distal Intergenic"]<- "Intergenic"
  me3.anno$annotation[me3.anno$annotation=="Promoter (<=1kb)"]<- "Promoter"
  me3.anno$annotation[me3.anno$annotation=="Promoter (1-2kb)"]<- "Promoter"
  return(me3.anno)
}

anno.me3.np<- lapply(names(k562.me3.np),H3K4me3.annotate)
saveRDS(anno.me3.np,file = "data/Roadmap_annotated_H3K4me3_np.rds")

Histone.annotate<- function(x){
  me3.anno<- anno.me3.np[[x]]
  if(x %in% names(k562.ac.np)){
    me3.anno$State<- factor("None",levels = c("Both","H3K4me1","H3K27ac","None"))
    me1.idx<- countOverlaps(me3.anno,k562.me1.np[[x]])>0
    ac.idx<-  countOverlaps(me3.anno,k562.ac.np[[x]])>0
    me3.anno$State[me1.idx]<- "H3K4me1"
    me3.anno$State[ac.idx]<- "H3K27ac"
    me3.anno$State[me1.idx & ac.idx]<- "Both"
  }else{
    me3.anno$State<- factor("None",levels = c("Both","H3K4me1","H3K27ac","None"))
    me1.idx<- countOverlaps(me3.anno,k562.me1.np[[x]])>0
    me3.anno$State[me1.idx]<- "H3K4me1"
  }
  return(me3.anno)
}

anno.me3.np<- lapply(names(anno.me3.np),Histone.annotate)
names(anno.me3.np)<- names(k562.me3.np)

countOverlaps(a,k562.me1.np$E001)
names(anno.me3.np)<- names(k562.me3.np)
anno.me3.np.d<- sapply(anno.me3.np, function(x)table(x$annotation))
anno.me3.np.d<- t(anno.me3.np.d)

anno.me3.np.s<- t(sapply(anno.me3.np,function(x)table(x$State[x$annotation!="Promoter"])))
gd<- as.data.frame(cbind(anno.me3.np.d,anno.me3.np.s))

roadmap.meta<- read.csv("/mnt/local-disk1/rsgeno2/MAmotif/Roadmap/jul2013.roadmapData.qc - Consolidated_EpigenomeIDs_summary_Table.csv")
roadmap.meta<- roadmap.meta[-(1:2),]
rownames(roadmap.meta)<- roadmap.meta$Epigenome.ID..EID.
gd$name<- roadmap.meta[names(anno.me3.np),"Epigenome.Mnemonic"]
gd<- gd[order(gd$Promoter),]
gd$name<- factor(gd$name,levels = gd$name)
ggplot(melt(gd[,5:9]))+geom_bar(aes(x=name,y=value,fill=variable),stat = "identity",position = "fill")+coord_flip()
