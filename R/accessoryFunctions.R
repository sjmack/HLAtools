##Accessory Functions v1.0.0 5FEB24

################
##posSort v2.0
#'Numerical sort for sequence alignment positions that contain indels.
#'
#'Sorts sequence alignment positions that contain indels in numerical order; e.g., three indels following position X are are identified as X.1, X.2 and X.3.
#'
#'@param posVec A vector of nucleotide positions
#'@param alignType The type of alignment the positions are in. The values 'prot', 'nuc' and 'gen' are valid options. 
#'@param locus A locus in the HLAalignments object for the specified alignType.
#'
#'@return A correctly sorted sequence.
#'
#'@export
#'
#'@examples
#'posSort(c(2,4,3,1,5), "nuc","DRB1")
posSort <- function(posVec,alignType, locus){
  if(!alignType %in% c("prot","nuc","gen")) {stop(paste(alignType,"is not a valid 'alignType' value. Please chose from 'prot', 'nuc' and 'gen'.\n",sep=" "))}
  if(!locus %in% names(HLAtools.data::HLAalignments[[alignType]])) {stop(paste(locus,"is not included among the",alignType,"alignments.\n",sep=" "))}
  
  tab <- as.data.frame(cbind(match(posVec,names(HLAtools.data::HLAalignments[[alignType]][[locus]])),posVec))
  tab <- tab[order(tab$V1),]
  tab$posVec
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
