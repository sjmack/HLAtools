#### Compare Allele Sequences v0.8.0 29FEB2024 

################
##compareSequences
#'Identify Sequence Differences Between Two Alleles at a Locus
#'
#'Compares the sequences of two alleles at a locus, and identifies the differences between them at specific positions
#'
#'@param alignType A character string identifying the type of alignment being searched. Allowed values are "codon","gen", nuc" and "prot". Only one 'alignType' value is permitted.
#'@param alleles A character vector containing two full-length names for alleles at the same locus.
#'
#'@return A two-row data frame identifying the positions and sequences at which the two alleles differ. E.g., compareSequences(alignType = "gen", alleles = c("DPA1*01:03:38:01","DPA1*01:03:38:02"). Positions for which the sequence of either allele is unknown are ignored. 
#'
#'@export
#'
compareSequences <- function(alignType,alleles) {
# Data checks
  if(length(alleles)!=2) {stop(paste("Please include exactly two alleles."))}
  
  if(length(alignType)!=1) {stop(paste("Please specify only one 'alignType'."))}
  
  alignType <- checkAlignType(alignType)
  
  if(strsplit(alleles[1],"*",fixed=TRUE)[[1]][1] != strsplit(alleles[2],"*",fixed=TRUE)[[1]][1]) {
    stop(paste(alleles[1],"and",alleles[2],"are alleles at different loci.\n",sep=" ")) }
  
  locus <- strsplit(alleles[1],"*",fixed=TRUE)[[1]][1]
  
  if(!validateLocus(locus,typeToSource(alignType))) { stop()}
    
  if(!validateAllele(alleles[1])) {stop(paste(alleles[1], "is not found in the genomic alignment for",locus,".",sep=" "))}
  if(!validateAllele(alleles[2])) {stop(paste(alleles[2], "is not found in the genomic alignment for",locus,".",sep=" "))}
  
  stopMessage <- ""
  stopErr <- FALSE
  for(i in 1:2) {
    if(!alleles[i] %in% HLAalignments[[alignType]][[locus]]$allele_name) {
      stopMessage <- paste(stopMessage,alleles[i]," is not a full-length allele name.\n",sep="")
      stopErr <- TRUE
    }
  }
  if(stopErr) {stop(paste("\n",stopMessage,"Please use full-length allele names.",sep="\n"))}
  
  if(alleles[1] == alleles[2]) { return(paste(alleles[1],"and",alleles[2],"are identical.",sep=" ")) }

# Function body
    coor1 <- HLAalignments[[alignType]][[locus]]$allele_name == alleles[1]
    coor2 <- HLAalignments[[alignType]][[locus]]$allele_name == alleles[2]

    align1 <- HLAalignments[[alignType]][[locus]][coor1,]
    align2 <- HLAalignments[[alignType]][[locus]][coor2,]

    align1 <- align1[4:length(align1)]
    align2 <- align2[4:length(align2)]

    noSeq <- (align1[1,] == "*") | (align2[1,] == "*")

    pairAlign <- rbind(align1[,colnames(align1[,noSeq == FALSE])],align2[,colnames(align1[,noSeq == FALSE])])

    variants <- pairAlign[1,] != pairAlign[2,] 

    diffs <- pairAlign[,variants]
    
    if(is.null(nrow(diffs))) { 
      
      diffs <- paste("There are no differences between", alleles[1], "and", alleles[2], "in the", switch(alignType, "codon" = "codon","nuc" = "nucleotide","prot"="protein","gen"="genomic"), "alignment.",sep=" ")
      
      } else {  rownames(diffs) <- 1:2 }

   diffs
}
