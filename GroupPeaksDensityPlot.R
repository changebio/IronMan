### =========================================================================
### plot the density of multiple type of data structure in different groups of peaks 
### -------------------------------------------------------------------------


# ---------------------------------------------------------------------------- #
#' Make ScoreMatrixList from multiple targets
#' 
#' GRanges >>bed,bam
#' Coverage >>bigwig,bedgraph
#' 
#' 
#' The function constructs a list of \code{ScoreMatrix} objects in the form
#' of \code{ScoreMatrixList} object. This object can be visualized using 
#' \code{multiHeatMatrix}, \code{heatMeta} or \code{plotMeta}
#'
#' @param targets can be a list of \code{scoreMatrix} objects, that are coerced 
#'        to the \code{ScoreMatrixList}, a list of \code{RleList} objects, or a 
#'        character vector specifying the locations of mulitple bam files  or
#'        bigWig files that 
#'        are used to construct the \code{scoreMatrixList}. If it is either a 
#'        RleList object or a character vector of files, it is obligatory to 
#'        give a windows argument.
#' @param windows \code{GenomicRanges} containing viewpoints for the scoreMatrix 
#'        or ScoreMatrixList functions
#' @param bin.num an integer telling the number of bins to bin the score matrix
#' @param bin.op an name of the function that will be used for smoothing windows of ranges
#' @param strand.aware a boolean telling the function whether to reverse the 
#'        coverage of ranges that come from - strand (e.g. when plotting 
#'        enrichment around transcription start sites)
#' @param weight.col if the object is \code{GRanges} object a numeric column
#'                 in meta data part can be used as weights. This is particularly
#'                useful when genomic regions have scores other than their
#'                coverage values, such as percent methylation, conservation
#'                scores, GC content, etc. 
#' @param is.noCovNA (Default:FALSE)
#'                  if TRUE,and if 'targets' is a GRanges object with 'weight.col'
#'                   provided, the bases that are uncovered will be preserved as
#'                   NA in the returned object. This useful for situations where
#'                   you can not have coverage all over the genome, such as CpG
#'                    methylation values.
#'                    
#' @param type if \code{targets} is a character vector of file paths, then type 
#'        designates the type of the corresponding files (bam or bigWig)
#' @param rpm boolean telling whether to normalize the coverage to per milion reads. 
#'            FALSE by default. See \code{library.size}.
#' @param unique boolean which tells the function to remove duplicated reads 
#'                       based on chr, start, end and strand
#' @param extend numeric which tells the function to extend the features
#'               ( i.e aligned reads) to total
#'               length ofwidth+extend
#' @param param ScanBamParam object  
#' @param library.size a numeric vector of the same length as \code{targets} 
#'                     indicating total number of mapped reads in BAM files (\code{targets}).
#'                     If is not given (default: NULL) then library sizes for every target
#'                     is calculated using the Rsamtools package functions:
#'                     param = ScanBamParam( flag = scanBamFlag(isUnmappedQuery=FALSE) )
#'                     sum(countBam(BamFile(target), param=param)$records).
#'                     \code{rpm} argument has to be set to TRUE.
#' @param cores the number of cores to use (default: 1)
#'
#' @return returns a \code{ScoreMatrixList} object
#' 

GroupPeaksDensityPlot <- function(target,gpeaks,bin.num=NULL,title="",){
  k4me3_NB.gro<- lapply(k4me3_NB,function(x)ScoreMatrixBin(gro.seq[[4]],x,bin.num = 50))
  k4me3_NB.gro1<- melt(sapply(k4me3_NB.gro,colMeans))
  ggplot(melt(k4me3_NB.gro1))+geom_line(aes(x=Var1,y=value,colour=Var2))+
    labs(x = "base",y = "density signal",title="The average GRO-Seq signal in different enhancer",fill="") +
    theme(plot.title = element_text(color="black", size=20, face="bold.italic"),
          axis.title.x = element_text( face="bold",size=14),
          axis.title.y = element_text(color="black", size=14, face="bold"),
          legend.title =element_text(face = "bold", size = 14, color = "black"),
          legend.text = element_text(face = "bold", size = 12),
          axis.text.x = element_text(face="bold",size=14),
          axis.text.y = element_text(face="bold", size=14)
    )
}


lapply(files,function(x)GroupPeaksDensityPlot(x,grl.k4me3.ma,))
