require(DOSE)
require(clusterProfiler)
library(org.Hs.eg.db)

ego<- lapply(1:44,function(i,x)enrichGO(gene = x[[i]]$geneId,organism = "human",ont = "CC",pvalueCutoff = 0.01,pAdjustMethod = "BH",qvalueCutoff = 0.05,readable = TRUE),x = dnase.K562.anno.sl)
saveRDS(ego,file = "data/Dnase_K562_CC_ego.rds")
ego1<- lapply(1:44,function(i,x)enrichGO(gene = x[[i]]$geneId,organism = "human",ont = "BP",pvalueCutoff = 0.01,pAdjustMethod = "BH",qvalueCutoff = 0.05,readable = TRUE),x = dnase.K562.anno.sl)
saveRDS(ego1,file = "data/Dnase_K562_BP_ego.rds")
ego2<- lapply(1:44,function(i,x)enrichGO(gene = x[[i]]$geneId,organism = "human",ont = "MF",pvalueCutoff = 0.01,pAdjustMethod = "BH",qvalueCutoff = 0.05,readable = TRUE),x = dnase.K562.anno.sl)
saveRDS(ego2,file = "data/Dnase_K562_MF_ego.rds")

a<- Reduce(rbind,lapply(ego2,function(x)x@result))
a$group<- rep(names(ego2),as.numeric(sapply(ego2,function(x)nrow(x@result))))
DF<- strsplit(as.character(a$GeneRatio),"/")
a$GeneRatio<- sapply(DF,function(x)as.numeric(x[1])/as.numeric(x[2]))
a$feature<-sapply(strsplit(a$group,"[.]"),function(x)x[1])
a$type<- sapply(strsplit(a$group,"[.]"),function(x)x[5])


ggplot(a[a$p.adjust<10**-5 & a$feature=="Distal Intergenic" & a$type=="K562",])+geom_point(aes(x=Description,y=GeneRatio,size=Count,colour=p.adjust,shape=type))+
  coord_flip()+
  facet_grid(.~type)+
  labs(x = "",y = " ",title="GO analysis for genes targeted by distal H3K4me3(p.adjust<10-5)") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=6),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size=16)
  )