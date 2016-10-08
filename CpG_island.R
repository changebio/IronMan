require(rtracklayer)
cpgm<- read.table("/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/cpgIslandExt.txt",sep = "\t")
colnames(cpgm)<- c("bin","chrom",	"chromStart","chromEnd","name","length","cpgNum","gcNum","perCpg","perGc","obsExp")
cpgm<- makeGRangesFromDataFrame(cpgm,keep.extra.columns = TRUE,starts.in.df.are.0based = TRUE)

cpgu<- read.table("/mnt/local-disk1/rsgeno2/huangyin/Rstudio/Iranman/cpgIslandExtUnmasked.txt",sep = "\t")
colnames(cpgu)<- c("bin","chrom",	"chromStart","chromEnd","name","length","cpgNum","gcNum","perCpg","perGc","obsExp")
cpgu<- makeGRangesFromDataFrame(cpgu,keep.extra.columns = TRUE,starts.in.df.are.0based = TRUE)

temp<- countOverlaps(k562.pk.dnase,cpgu)
table(temp,k562.pk.dnase$State)
