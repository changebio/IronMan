

### CAGE or GRO-Seq signal in FANTOM and ENCODE
dnase.K562.sg<- sapply(dnase.K562.anno.sl,function(x)region.cage.signal(x,cage.seq,strand=FALSE,weight.col = "V5"))

temp<- lapply(seq(1,length(dnase.K562.sg),by=6),function(i){
  a1<- as.data.frame(dnase.K562.sg[i:(i+5)])
  colnames(a1)<- names(gro.K562.sg)[1:6]
  return(a1)
})
names(temp)<-names(dnase.K562.anno.sl)

temp<- Reduce(rbind,temp)
temp$region<- rep(names(dnase.K562.anno.sl),each=1300)
dnase.K562.sg<- temp
saveRDS(dnase.K562.sg,file = "data/Dnase_K562_sg.rds")

gro.K562.sg<- lapply(dnase.K562.anno.sl,function(x)region.cage.signal(x,gro.seq,strand=FALSE,weight.col = "V5"))
temp<- lapply(seq(1,length(gro.K562.sg),by=6),function(i){
a1<- as.data.frame(gro.K562.sg[i:(i+5)])
colnames(a1)<- names(a)
a1$type<- substring(a1$type,first = 1,last = nchar(as.character(a1$type))-9)
return(a1)
})
names(temp)<-names(dnase.K562.anno.sl)

temp<- Reduce(rbind,temp)
temp$region<- rep(names(dnase.K562.anno.sl),each=1300)

gro.K562.sg<- temp
saveRDS(gro.K562.sg,file = "data/GRO_K562_sg.rds")

ggplot(temp)+geom_line(aes(x=Var1,y=value,linetype=as.factor(read),colour=region))+
  facet_grid(.~type)+
  labs(x = "",y = " ",title="The average signal in annotated DNase",linetype="read") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size = 16)
  )



all.K562.sg<- lapply(dnase.K562.anno.sl,function(x)region.cage.signal(x,list(all.gr.m,all.gr.p),strand=FALSE,weight.col = "V5"))
saveRDS(all.K562.sg,file="data/ALL_K562_sg.rds")
temp<- Reduce(rbind,all.K562.sg)
temp$type<- "K562_CTSS_pool"
temp$region<- rep(names(dnase.K562.anno.sl),each=100)


temp<- rbind(temp,dnase.K562.sg)
temp<- rbind(temp,gro.K562.sg)
temp$value[temp$value>10]<- 10

ggplot(temp)+geom_line(aes(x=Var1,y=value,linetype=as.factor(read),colour=type))+
  facet_wrap(~region)+
  labs(x = "",y = " ",title="The average signal in annotated DNase(CTSS pool)",linetype="read") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
        axis.title.x = element_text( face="bold",size=14),
        axis.title.y = element_text(color="black", size=14, face="bold"),
        legend.title =element_text(face = "bold", size = 14, color = "black"),
        legend.text = element_text(face = "bold", size = 12),
        axis.text.x = element_text(face="bold",size=14),
        axis.text.y = element_text(face="bold", size=14),
        strip.text.x = element_text(face = "bold",size = 8)
  )
