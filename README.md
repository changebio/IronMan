# IronMan
study non-promoter H3K4me3 in K562
first part: introduction
there are many study of H3K4me3 in promoter region. but almost didn't find the study of H3K4me3. we found that most H3K4me3 are in Promoter region
in hESC, and there are an increase in non-promoter region in other cell lines. many of these non-promoter H3K4me3 co-occur with H3K4me1 and H3K27ac, which usually defined as active enhancer. they are short than the H3K4me3 peaks in promoter, which have studied that the genes with long H3K4me3 have stable expression. And they are cell-type specific.
second part: ~~~~~~~~~~~~~~~~~~


1.genomic distribution in different cell lines(so I need their peaks by runing Macs under the same condition)

2.seperate peaks as promoter peaks and non-promoter peaks, look into their overlapping in different cell lines, respectively.

3. motifscan in non-promoter peaks to study the relationship between H3K4me3 and RACK7

4. use K562 and H1 as example, K562-only non-promoter enhancer H3K4me3 peak(1 increase a lot non-promoter H3K4me3, 2 most of the non-promoter H3K4me3 peaks are in enhancer)

#Data~~~~~~read all bed-like file by GRange format(1 base)
*macs output file xls(1 base),bed(0 base)*
*macs2 output file 
ChIP-Seq
H3K4me3(H1,K562,Helas3,Hmec,Nha,Nhdfad,Nhek,Nhlf)
H3K4me1(H1,K562)
H3K27ac(H1,K562,Helas3,Hmec,Nhek)

RNA-Seq

#report
three part
H3K4me3 distribution with the same p-value cutoff among different cell types
H3K4me3 distribution with the same peak numbers among different cell types

H3K4me3 common in promoter and specific in non-promoter





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
  
  

##20160811
prediction using mechine learning 

first, try to use multiple linear regression to predict region character by TFBS info.
  -define region
  -annotate these region
  -predict 

Protocols used in RNA Extraction
================================

Long Cytoplasmic RNA A+ and A-:  Long Cytoplasmic RNA was isolated using
Qiagens RNeasy Protocol.  Qiagens Qligotex kit was used to separate
Poly-A(+) from Poly-A(-).

Long Nuclear RNA A+ and A-: Cell were lysed with Qiagens RLN buffer and the
nuclei were spun out, resuspended in Qiagen's RLT buffer, layered over a
Cesium Chloride bed and spun at 24,000 rpm for a minimum of 20 hours. The
pellet was recovered and cleaned up on Qiagens RNeasy kit. Qiagens Qligotex
kit was used to separate Poly-A(+) from Poly-A(-).

Long Polysomal RNA A+ and A-: Cyclohexamide arrested cell lysate were
layered on a 20-60% sucrose gradient and spun at 27,000 rpm for 4 hours.
The fractions and analyzed by UV-spectroscopy using a programmable Density
Gradient Fractionation System Foxy Jr. (Isco).

Long Nucleoplasm, Nucleoli and Chromatin: We used the protocol found in
Bhorjee and Pederson (1973). 

Small Cytoplasmic RNA: Small Cytoplasmic RNA was isolated using Qiagens
RNA/DNA kit.

Small Nuclear RNA:  Small Nuclear RNA was isolated using Qiagens RNA/DNA
kit.


