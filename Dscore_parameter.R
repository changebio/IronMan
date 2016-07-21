#find the best region size and the best cutoff of TSS count

#load package-------
require(reshape2)
require(ggplot2)

#load data-------


ucsc.gene.pt<- promoters(ucsc.gene,upstream = 2048,downstream = 2048)

ucsc.gene.pt.s<- ScoreMatrixList(gro.seq[5:6],ucsc.gene.pt,weight.col = "V5")
ucsc.gene.pt <- ucsc.gene.pt[as.numeric(Reduce(intersect,lapply(ucsc.gene.pt.s,rownames)))]
ucsc.gene.pt.ml<- intersectScoreMatrixList(ucsc.gene.pt.s)

# a<-sapply(0:11,function(n){
#   score<- lapply(ucsc.gene.pt.ml,function(x)rowSums(x[,(2049-2^n):(2048+2^n)]))
#   temp<- (score[[2]]-score[[1]])/(score[[2]]+score[[1]])
# })

a<-sapply(0:11,function(n){
  score<- lapply(ucsc.gene.pt.ml,function(x)rowSums(x[,(2049-2^n):(2048+2^n)]))
  temp<- rep(NA,length(ucsc.gene.pt))
  ind<- (score[[2]]+score[[1]])>=6
  score <- lapply(score, function(x)x[ind])
  temp[ind]<- (score[[2]]-score[[1]])/(score[[2]]+score[[1]])
  temp
})

Dscore<- as.data.frame(t(apply(a,2,function(x)table(x>0,as.data.frame(ucsc.gene.pt)$strand))))
Dscore$number<- apply(a,2,function(x)(length(ucsc.gene.pt)-sum(is.na(x)))/length(ucsc.gene.pt))
Dscore$plus<- Dscore$V2/(Dscore$V1+Dscore$V2)
Dscore$minus<- Dscore$V3/(Dscore$V3+Dscore$V4)
Dscore$rate <- (Dscore$V2+Dscore$V3)/rowSums(Dscore[,1:4])
Dscore$region<- as.factor(1:12)


ggplot(melt(Dscore[7:11]))+geom_line(aes(x=2^as.numeric(region),y=value,colour=variable))+
  labs(x = "",y = " ",title="ratio",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )


#the different cutoff of cage read----------
a<-sapply(1:20,function(n){
  score<- lapply(ucsc.gene.pt.ml,function(x)rowSums(x[,(2049-50):(2048+50)]))
  temp<- rep(NA,length(ucsc.gene.pt))
  ind<- (score[[2]]+score[[1]])>=n
  score <- lapply(score, function(x)x[ind])
  temp[ind]<- (score[[2]]-score[[1]])/(score[[2]]+score[[1]])
  temp
})

Dscore<- as.data.frame(t(apply(a,2,function(x)table(x>0,as.data.frame(ucsc.gene.pt)$strand))))
Dscore$number<- apply(a,2,function(x)(length(ucsc.gene.pt)-sum(is.na(x)))/length(ucsc.gene.pt))
Dscore$plus<- Dscore$V2/(Dscore$V1+Dscore$V2)
Dscore$minus<- Dscore$V3/(Dscore$V3+Dscore$V4)
Dscore$rate <- (Dscore$V2+Dscore$V3)/rowSums(Dscore[,1:4])
Dscore$region<- as.factor(1:20)

require(reshape2)
require(ggplot2)
ggplot(melt(Dscore[7:11]))+geom_line(aes(x=as.numeric(region),y=value,colour=variable))+
  labs(x = "",y = " ",title="ratio",fill="") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14)
  )




