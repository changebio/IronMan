require(rtracklayer)
cpgm<- read.table("/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/cpgIslandExt.txt",sep = "\t")
colnames(cpgm)<- c("bin","chrom",	"chromStart","chromEnd","name","length","cpgNum","gcNum","perCpg","perGc","obsExp")
cpgm<- makeGRangesFromDataFrame(cpgm,keep.extra.columns = TRUE,starts.in.df.are.0based = TRUE)

cpgu<- read.table("/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/cpgIslandExtUnmasked.txt",sep = "\t")
colnames(cpgu)<- c("bin","chrom",	"chromStart","chromEnd","name","length","cpgNum","gcNum","perCpg","perGc","obsExp")
cpgu<- makeGRangesFromDataFrame(cpgu,keep.extra.columns = TRUE,starts.in.df.are.0based = TRUE)

temp<- countOverlaps(k562.pk.dnase,cpgu)
table(temp,k562.pk.dnase$State)

#chromHMM
k562.hmm<- import.bed("/mnt/local-disk1/rsgeno2/MAmotif/ENCODE/wgEncodeAwgSegmentationSegwayK562.bed")

a<- findOverlaps(k562.pk.dnase,k562.hmm,select = "arbitrary")
table(k562.hmm$name[a],k562.pk.dnase$State)

##
k562.tf<- import.bed("/mnt/local-disk1/rsgeno2/MAmotif/ENCODE/ENCODE_AwgTfbs.hg19.txt")
table(countOverlaps(k562.pk.dnase,k562.tf),k562.pk.dnase$)
