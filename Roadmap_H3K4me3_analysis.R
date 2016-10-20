
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