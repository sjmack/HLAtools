#### Compare Allele Sequences v1.0.1 May 29, 2026 

################
##compareSequences
#'Identify Sequence Differences Between Two Alleles
#'
#'Compares the sequences of two alleles, and identifies the differences between them at shared polymorphic sequence positions.
#'Alleles at different loci can be compared. 
#'
#'@param alignType A character string identifying the type of alignment being searched. Allowed values are "codon","gen", nuc" and "prot". Only one 'alignType' value is permitted.
#'@param alleles A character vector containing two full-length names for alleles at the two different loci.
#'
#'@return A two-row data frame identifying the positions and sequences at which the two alleles differ. E.g., compareLoci(alignType = "gen", alleles = c("B\*07:02:01:01","C\*02:02:02:01")). Positions are nor shared between the two loci, or for which the sequence of either allele is unknown are ignored. 
#'
#'@examples
#'\dontrun{
#'compareSequences("gen",c("DPA1*01:03:38:01","DPA1*01:03:38:02"))
#'compareSequences("prot",c("DRB1*01:01:01:01","DRB3*01:01:02:01"))
#'compareSequences("prot", c("B*57:01:01:01", "A*24:02:01:01"))
#'}
#'
#'@export
#'
compareSequences <- function(alignType,alleles) {
  # Data checks
  if(length(alleles)!=2) {stop(paste("Please include exactly two alleles."))}
  
  if(length(alignType)!=1) {stop(paste("Please specify only one 'alignType'."))}
  
  alignType <- checkAlignType(alignType)
  
  locus <- c(alleles[1],alleles[2])
  
  locus[1] <- strsplit(alleles[1],"*",fixed=TRUE)[[1]][1]
  locus[2] <- strsplit(alleles[2],"*",fixed=TRUE)[[1]][1]
  
  if(!validateLocus(locus[1],typeToSource(alignType))) { stop()}
  if(!validateLocus(locus[2],typeToSource(alignType))) { stop()}
  
  if(!validateAllele(alleles[1])) {return(paste(alleles[1], "is not found in the", alignType, "alignment for",paste(locus[1],".",sep=""),sep=" "))}
  if(!validateAllele(alleles[2])) {return(paste(alleles[2], "is not found in the", alignType, "alignment for",paste(locus[2],".",sep=""),sep=" "))}
  
  stopMessage <- ""
  stopErr <- FALSE
  
  for(i in 1:2) {
    if(!alleles[i] %in% HLAalignments[[alignType]][[locus[i]]]$allele_name) {
      stopMessage <- paste(stopMessage,alleles[i]," is not a full-length allele name.\n",sep="")
      stopErr <- TRUE
    }
  }
  
  if(stopErr) {stop(paste("\n",stopMessage,"Please use full-length allele names.",sep="\n"))}
  
  if(alleles[1] == alleles[2]) { return(paste(alleles[1],"and",alleles[2],"are identical.",sep=" ")) }
  
  # Function body
  coor1 <- HLAalignments[[alignType]][[locus[1]]]$allele_name == alleles[1]
  coor2 <- HLAalignments[[alignType]][[locus[2]]]$allele_name == alleles[2]
  
  align1 <- HLAalignments[[alignType]][[locus[1]]][coor1,]
  align2 <- HLAalignments[[alignType]][[locus[2]]][coor2,]
  
  align1 <- align1[4:length(align1)]
  align2 <- align2[4:length(align2)]
  
  align1 <- align1[colnames(align1) %in% colnames(align2)]
  align2 <- align2[colnames(align2) %in% colnames(align1)]
  
  noSeq <- (align1[1,] == "*") | (align2[1,] == "*")
  
  combiPair <- rbind(align1,align2)
  rownames(combiPair) <- c(1,2)
  
  combiPair <- combiPair[,noSeq == F]
  combiPair <- combiPair[,combiPair[1,] != combiPair[2,]]
  
  if(length(combiPair[,combiPair[1,] != combiPair[2,]]) == 2) {return(paste("There are no differences between", alleles[1], "and", alleles[2], "in the", switch(alignType, "codon" = "codon","nuc" = "nucleotide","prot"="protein","gen"="genomic"), "alignment.",sep=" "))}
  
  return(combiPair[,combiPair[1,] != combiPair[2,]])
  
}
