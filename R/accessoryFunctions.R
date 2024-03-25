##Accessory Functions v4.0.0 20MAR24

################
##posSort
#'Numerical sort for sequence alignment positions that contain indels.
#'
#'@description
#'Sorts sequence alignment positions that contain indels in numerical order; e.g., three indels following position X are are identified as X.1, X.2 and X.3.
#'
#'Positions are validated against a specified locus for a specified alignment type. Positions that do not exist for that locus-alignment combination are not returned.
#'
#'@param posVec A character vector of variant positions
#'@param alignType The type of alignment the positions are in. The values 'prot', 'codon', 'nuc' and 'gen' are valid options. 
#'@param locus A locus in the HLAalignments object for the specified alignType.
#'
#'@return A correctly sorted sequence.
#'
#'@note Indel positions must be text-formatted (e.g. "607.12"), as illustrated in the the examples.
#'
#'@export
#'
#'@examples
#'\dontrun{
#'posSort(c(2,4,3,1,5), "nuc","DRB1")
#'posSort(c("607.23","607.10","607.3","607.4"),"nuc","DRB1")
#'}
posSort <- function(posVec,alignType, locus){
  if(!alignType %in% c("prot","codon", "nuc","gen")) {stop(paste(alignType,"is not a valid 'alignType' value. Please chose from 'prot', 'codon', nuc' and 'gen'.\n",sep=" "))}
  if(!locus %in% names(HLAalignments[[alignType]])) {stop(paste(locus,"is not included among the",alignType,"alignments.\n",sep=" "))}
  
  tab <- as.data.frame(cbind(match(posVec,names(HLAalignments[[alignType]][[locus]])),posVec))
  tab$V1 <- as.numeric(tab$V1)
  tab <- tab[order(tab$V1),]
  tab$posVec[!is.na(tab$V1)]
}

################
##numFields
#'Identifies the number of fields in a colon-delimited HLA allele name.
#'
#'Returns the number of fields in a colon-delimited HLA allele name. A value of 1 is returned for digit-delimited HLA allele names.
#'
#'@param allele A colon-delimited HLA allele name.
#'
#'@return The number of fields in the allele name.
#'
#'@export
#'
#'@examples
#'numFields("HLA-A*01:01")
#'numFields("DRB1*04:03:01")
#'numFields("13:02:01:01")
numFields <- function(allele) {
  length(strsplit(allele,":",fixed=TRUE)[[1]])
}

#################
##validateAllele
#'Determine if an allele name is properly formed and present in HLAalignments
#'
#'Returns TRUE if an allele name is found in HLAalignments in either the 'allele_name' column of full-length allele names or the 'trimmed_allele' column of two-field allele names in the pertinent genomic alignment. Returns FALSE if the allele name is not properly formed, or if the allele name is not found in HLAalignments.
#'
#'@param allele A colon-delimited HLA allele name.
#'
#'@return A logical identifying if the allele name is present in the alignments (TRUE) or, if it is not in the alignments or is not valid not (FALSE).
#'
#'@note Messages will be returned to the console if the allele name is malformed, or the locus is invalid. 
#'@note The locus being evaluated must be included in HLAalignments.
#'
#'@export
#'
#'@examples
#'\dontrun{
#'validateAllele("A*01:01:01:117")
#'validateAllele("A*01:01:01")
#'validateAllele("A*01:01")
#'}
validateAllele <- function(allele) {
  alleleParts <- strsplit(allele,"*",fixed=TRUE)[[1]]
  if(length(alleleParts) != 2) {
    message(paste("No asterisk ('*') is present in ",allele,".",sep=""))
    return(FALSE)
    }
  if(length(strsplit(alleleParts[2],":",fixed=TRUE)[[1]]) == 1) {
    message(paste("No colon (':') is present in ",allele,".",sep=""))
    return(FALSE)
  }
  locus <- alleleParts[1]
  if(!suppressMessages(validateLocus(alleleParts[1],"gDNA")) && !suppressMessages(validateLocus(alleleParts[1],"cDNA")) && !suppressMessages(validateLocus(alleleParts[1],"prot"))) {
    message(paste(locus, "is not a valid locus.", sep=" "))
    return(FALSE)
  } 
  if((allele %in% HLAalignments$prot[[locus]]$trimmed_allele) || (allele %in% HLAalignments$nuc[[locus]]$trimmed_allele) || (allele %in% HLAalignments$gen[[locus]]$trimmed_allele) || (allele %in% HLAalignments$prot[[locus]]$allele_name) || (allele %in% HLAalignments$nuc[[locus]]$allele_name) || (allele %in% HLAalignments$gen[[locus]]$allele_name)) {
    return(TRUE) } else { message(paste(allele, "is not found in version", HLAalignments$version,"alignments.",sep=" "))
    return(FALSE)
   }
}


#################
##verifyAllele
#'Determine if a specific allele name ever existed, and (if so) the most recent IPD-IMGT/HLA Database release in which it appeared
#'
#'Returns TRUE if an allele name is present in AlleleListHistory or FALSE it is absent, or c(TRUE,version), where 'version' is the most recent IPD-IMGT/HLA Database release in which that name appeared, when version = TRUE.
#'
#'@param allele An HLA allele name. Colon-delimited and field-delimited names are both accepted.
#'@param version A logical that indicates if the most recent nomenclature release version in which that name was valid should be returned. 
#'
#'@return A logical identifying if the allele name is found in AlleleListHistory (TRUE) or not (FALSE), or c(TRUE,version) if version = TRUE.
#'
#'@note
#'
#'@export
#'
#'@examples
#'\dontrun{
#'verifyAllele("A*01:01:01:01")
#'verifyAllele("A*01:01:01:01",TRUE)
#'verifyAllele("A*010101",TRUE)
#'verifyAllele("A*0101",TRUE)
#'}
verifyAllele <- function(allele, version=FALSE){

    resArray <- which(alleleListHistory$AlleleListHistory == allele,arr.ind = TRUE)
    if(length(resArray)==0) {return(FALSE)}
      
    if(!version) {return(TRUE)
        } else {
          rawVersion <- colnames(alleleListHistory$AlleleListHistory)[resArray[1,2][[1]]]
            return(c(TRUE,expandVersion(substr(rawVersion,start=2,stop = nchar(rawVersion)))))
          } 
  }
    
    
    

  
  

