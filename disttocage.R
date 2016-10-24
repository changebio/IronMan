#get distance to the nearest cage ------
disttocage<- function(cage.p,one.p){
  mcols(cage.p) <- cbind(mcols(cage.p),DataFrame(TSS=as.data.frame(distanceToNearest(cage.p,one.p,ignore.strand= TRUE))$distance,PLUS=as.data.frame(distanceToNearest(cage.p,one.p))$distance))
  cage.p$PLUS[cage.p$TSS!=cage.p$PLUS]<- -cage.p$TSS[cage.p$TSS!=cage.p$PLUS]
  cage.p <- cage.p[abs(cage.p$PLUS)<251]
  return(cage.p)
}


#cage postions------
cage.p<- gro.seq$K562_biol_rep123.hg19.ctss_plus_pool.bed
cage.p<- cage.p[cage.p$V5>3]
strand(cage.p)<- cage.p$V6

# txdb tss-----------
txdb.tps<-  transcripts(txdb, columns=c("tx_id", "tx_name"))
txdb.tss <- resize(txdb.tps,width = 1,fix = "start")
txdb.tss <- unique(txdb.tss)

txdb.tss.dt<- disttocage(cage.p,txdb.tss)
hist(txdb.tss.dt$PLUS)
##Plus direction genes and minus direction genes-----
txdb.tss.dtp<- disttocage(cage.p,txdb.tss[strand(txdb.tss)=="+"])
hist(txdb.tss.dtp$PLUS)

txdb.tss.dtm<- as.data.frame(distanceToNearest(cage.p,txdb.tss[strand(txdb.tss)=="-"],ignore.strand= TRUE))
hist(txdb.tss.dtm$distance[txdb.tss.dtm$distance<1000])

#combined TSS by merging hg19_Ensembl_Genes, hg19_GENCODE_V14,knownGene, and refGene -----
comb.tss<- resize(CTSS,width = 1,fix="start")
comb.tss<- unique(comb.tss)
comb.tss.dt<- disttocage(cage.p,comb.tss)
hist(comb.tss.dt$PLUS)

#DNase distribution around cage--------
dnase.K562.ct<- unique(resize(dnase.K562,width = 1,fix = "center"))
dnase.K562.dt<- disttocage(cage.p,dnase.K562.ct)
hist(dnase.K562.dt$PLUS)





