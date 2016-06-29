# IronMan
study non-promoter H3K4me3 in K562
1.genomic distribution in different cell lines(so I need their peaks by runing Macs under the same condition)

2.seperate peaks as promoter peaks and non-promoter peaks, look into their overlapping in different cell lines, respectively.

3. motifscan in non-promoter peaks to study the relationship between H3K4me3 and RACK7

4. use K562 and H1 as example, K562-only non-promoter enhancer H3K4me3 peak(1 increase a lot non-promoter H3K4me3, 2 most of the non-promoter H3K4me3 peaks are in enhancer)

#Data
ChIP-Seq
H3K4me3(H1,K562,Helas3,Hmec,Nha,Nhdfad,Nhek,Nhlf)
H3K4me1(H1,K562)
H3K27ac(H1,K562,Helas3,Hmec,Nhek)

RNA-Seq

----------
##20160531
two parts
between cancer and normal
in cancer

##20160602
How to prove overactivation?
use H3K4me3 to define enhancer, enhancer with H3K27ac as active enhancer, enhancer with H3K27me3 as poised enhancer
enhancer with H3K27ac and H3K4me3 as overactive enhancer, enhancer with H3K27me3 and H3K4me3 as what???
* this analysis should be done in one cell line?(I think so, cause many enhancers are cell-type specific. if you do want to compare the enhancer's activity in the same region. first, we should control the variables, which means the H3K4me1 and H3K27ac signal should be the same[M value <1].then, compare gene expressions targeted by enhancer with or without H3K4me3)
    not as expected. Maybe the hypothesis is not right.
    or 
    use chiapet data to define target genes.(use the region interaction to annotate target genes)
    
    
    
* different genes can be comparable?(if you suppose that gene expressions depend on enhancer's activity, they are comparable)
    -in order to do that, we need to choose the genuine peaks. there are two ways, overlaping or M-value cutoff(which is best???how to evaluate???what about their comsistency?)


super-enhancers were found at key oncogenes in many cancer types.(15)
  Do enhancer with H3K4me3 tend to include in super-enhancer?
  
  


Error : package ‘IRanges’ 2.2.9 was found, but >= 2.3.7 is required by ‘Rsamtools’
ERROR: lazy loading failed for package ‘Rsamtools’
* removing ‘/opt/rstudio/rstudio1/R/x86_64-pc-linux-gnu-library/3.2/Rsamtools’
* restoring previous ‘/opt/rstudio/rstudio1/R/x86_64-pc-linux-gnu-library/3.2/Rsamtools’
ERROR: dependency ‘rtracklayer’ is not available for package ‘ChIPseeker’
* removing ‘/opt/rstudio/rstudio1/R/x86_64-pc-linux-gnu-library/3.2/ChIPseeker’
ERROR: dependency ‘SummarizedExperiment’ is not available for package ‘GenomicAlignments’
* removing ‘/opt/rstudio/rstudio1/R/x86_64-pc-linux-gnu-library/3.2/GenomicAlignments’
