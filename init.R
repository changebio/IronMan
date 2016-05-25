# setting color -----------------------------------------------------------

tropical=  c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)
par(pch=19)
# load refgene hg19 -------------------------------------------------------

load("/opt/rstudio/rstudio1/rstudio1/proj/refdata/refhg19.RData")
gtf<- read.table("/mnt/local-disk1/rsgeno2/huangyin/PRC2/Refgenes/genes.gtf",sep = "\t")

# load hg19 cpg island  ---------------------------------------------------

cpghg19<- read.table("/mnt/local-disk1/rsgeno2/huangyin/PRC2/hg19.cpgIslandExt.txt",header=TRUE,row.names = NULL,sep = "\t")
cpghg19<- cpghg19[,-1]

# load two H3K4me3 ChIP Seq data from ENCODE ------------------------------
h3k4me3.files <- list.files("/mnt/local-disk1/rsgeno2/MAmotif/macs.20151224/H3K4me3", pattern = "xls",include.dirs = TRUE)
h3k4me3<- lapply(paste0("/mnt/local-disk1/rsgeno2/MAmotif/macs.20151224/H3K4me3/",h3k4me3.files),function(x)read.table(x,header = TRUE))
names(h3k4me3)<- substr(h3k4me3.files,1,nchar(h3k4me3.files)-24)
h3k4me3<- lapply(h3k4me3,function(x)Peaks.Distribution(x,refhg19))

# load two H3K27ac ChIP Seq data from ENCODE ------------------------------
h3k27ac.files <- list.files("/mnt/local-disk1/rsgeno2/MAmotif/macs.20151224/H3K27ac", pattern = "xls",include.dirs = TRUE)
h3k27ac<- lapply(paste0("/mnt/local-disk1/rsgeno2/MAmotif/macs.20151224/H3K27ac/",h3k27ac.files),function(x)read.table(x,header = TRUE))
names(h3k27ac)<- substr(h3k27ac.files,1,nchar(h3k27ac.files)-24)
h3k27ac<- lapply(h3k27ac,function(x)Peaks.Distribution(x,refhg19))

# load two H3K4me1 ChIP Seq data from ENCODE ------------------------------
h3k4me1.files <- list.files("/mnt/local-disk1/rsgeno2/MAmotif/macs.20151224/H3K4me1", pattern = "xls",include.dirs = TRUE)
h3k4me1<- lapply(paste0("/mnt/local-disk1/rsgeno2/MAmotif/macs.20151224/h3K4me1/",h3k4me1.files),function(x)read.table(x,header = TRUE))
names(h3k4me1)<- substr(h3k4me1.files,1,nchar(h3k4me1.files)-24)
h3k4me1<- lapply(h3k4me1,function(x)Peaks.Distribution(x,refhg19))


