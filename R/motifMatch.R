## motifMatch v0.2 26May2024 

####################
##MotifMatch
#'Identify Alleles that Share a Sequence Motif
#'
#'Searches alignments for alleles that share a specific sequence motif. 
#'
#'@param motif A character string identifying a sequence variant motif in the following format: Locus*#$~#$~#$, where ## identifies a variant position, and $ identifies the sequence variant. Both nucleotide and peptide motifs can be provided, and any number of variants can be specified.
#'@param alignType A character string identifying the type of alignment being searched. Allowed values are "codon","gen", nuc" and "prot". Only one 'alignType' value is allowed.
#'@param full A logical value that specifies if full (TRUE) or truncated (FALSE) allele names should be returned.
#'
#'@return A character vector of allele names that contain the motif, or NA when no alleles contain the motif, NULL when no alignment is available for specified locus, and FALSE when the locus or motif is invalid.
#'
#'@note This function requires an HLAalignments object that has been populated with alignments via alignmentFull().
#'
#'@export
#'
motifMatch <- function(motif,alignType,full=TRUE){  
  
        aName <- 3 # return truncated names
            if(full == TRUE) {aName <- 4 } # return full names 
        
        if(!length(alignType) == 1) {stop("Please provide a single 'alignType' value.")}
        
        if(!validateMotif(motif,alignType)) {return(FALSE)}
  
        splitMotif <- strsplit(motif,"*",fixed = TRUE)[[1]]

        loc <- splitMotif[1]
 
        motifs <- strsplit(splitMotif[2],"~",fixed=TRUE)[[1]]
        numMotif <- length(motifs)
        splitMotifs <- strsplit(motifs,"",fixed=FALSE)
 
        posList <- as.list(rep(NA,numMotif)) ## 
        seqList <- as.list(rep(NA,numMotif)) ## may not need these
 
            for(i in 1:numMotif){
 
                  seqVar <- splitMotifs[[i]][splitMotifs[[i]] %in% c(letters,LETTERS) == TRUE]
                  seqPos <- paste(splitMotifs[[i]][!splitMotifs[[i]] %in% c(letters,LETTERS) == TRUE],collapse="")
   
                  posList[[i]] <- HLAalignments[[alignType]][[loc]][HLAalignments[[alignType]][[loc]][,colnames(HLAalignments[[alignType]][[loc]]) %in% seqPos] == seqVar,aName]
             }
  
            matchList <- Reduce(intersect,posList)
  
            if(identical(matchList,character(0))) { 
                  
                  message(paste("The",motif,"motif was not found in any",alignType,"alleles.",sep=" "))
                  return(NA)
            } else {
                matchList
        }
  
}

####################
### Validate Motif
#' Determine if a Motif is Properly Formatted
#'
#'Evaluates a motif to determine if the locus is valid and variants are valid. 
#'
#'@param motif A character string identifying a sequence variant motif in the following format: Locus*#$~#$~#$, where ## identifies a variant position, and $ identifies the sequence variant. Both nucleotide and peptide motifs can be provided, and any number of variants can be specified.
#'@param alignType A character string identifying the type of alignment being searched. Allowed values are "codon","gen", nuc" and "prot". Only one 'alignType' value is allowed.
#'
#'@return If the motif is valid, TRUE is returned. If the locus or body of the motif is invalid, FALSE is returned along with a brief message.
#'
#'@examples
#' validateMotif("A*-21M~2P","prot")
#' validateMotif("A*196G~301A~3046T","gen")
#'
#'@export
#'
validateMotif <- function(motif,alignType){

    alignType <- checkAlignType(alignType)
  
    splitMotif <- strsplit(motif,"*",fixed = TRUE)[[1]]
  
    if(!validateLocus(source = typeToSource(alignType),splitMotif[1])) {
      return(FALSE)
    }
  
    if(alignType == "prot") {
          if(any(strsplit(splitMotif[2],"")[[1]] %in% c(letters,":","*","B","J","O","U","Z","+","!","@","#","$","%","^","&","(",")","[","]","{","}",",","|","<",">","?"))) {
                message(paste(motif,"contains invalid characters.",sep = " "))
                return(FALSE)
                }
            } else { ## for "nuc" and "gen" alignments
              if(any(strsplit(splitMotif[2],"")[[1]] %in% c(letters,":","*","B","D","E","F","H","I","J","K","L","M","N","O","P","Q","R","S","U","V","W","X","Y","Z","+","!","@","#","$","%","^","&","(",")","[","]","{","}",",","|","<",">","?"))) {
                    message(paste(motif,"contains invalid characters.",sep = " "))
                    return(FALSE)
              }
      
      
        }
    
    return(TRUE)
} 
