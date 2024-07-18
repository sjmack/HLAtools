##Accessory Functions v5.0.0 9APR2024

################
##posSort
#'Numerical Sort of Alignment Positions that Contain Indels
#'
#'@description
#'Sorts sequence alignment positions that contain indels in numerical order; e.g., three indels following position X are are identified as X.1, X.2 and X.3.
#'
#'Positions are validated against a specified locus for a specified alignment type. Positions that do not exist for that locus-alignment combination are not returned.
#'
#'@param posVec A character vector of variant positions.
#'@param alignType A character string describing the type of alignment the positions are found in. The values 'prot', 'codon', 'nuc' and 'gen' are valid options. Only one 'alignType' value is allowed.
#'@param locus A character string describing a locus in the HLAalignments object for the specified alignType.
#'
#'@return A character string of the correctly sorted sequence.
#'
#'@note Indel positions must be text-formatted (e.g. "607.12"). C.f., posSort(c(2,4,3,1,5), "nuc","DRB1") vs posSort(c("607.23","607.10","607.3","607.4"),"nuc","DRB1").
#'
#'@export
#'
posSort <- function(posVec,alignType, locus){
  
    alignType <- checkAlignType(alignType)  
    
    if(length(alignType)!=1) {stop(paste("Please specify only one 'alignType'."))}
    
    if(!locus %in% names(HLAalignments[[alignType]])) {stop(paste(locus,"is not included among the",alignType,"alignments.\n",sep=" "))}
  
  tab <- as.data.frame(cbind(match(posVec,names(HLAalignments[[alignType]][[locus]])),posVec))
  tab$V1 <- as.numeric(tab$V1)
  tab <- tab[order(tab$V1),]
  tab$posVec[!is.na(tab$V1)]
}

################
##numFields
#'Identify the Number of Fields in a Colon-Delimited Allele Name
#'
#'Returns the number of fields in a colon-delimited HLA allele name. A value of 1 is returned for digit-delimited HLA allele names.
#'
#'@param allele A character string of a colon-delimited HLA allele name.
#'
#'@return A numeric value describing the number of fields in the allele name.
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
#'Validate Allele-Name Format and Presence in HLAalignments
#'
#'Returns TRUE if an allele name is found in HLAalignments in either the 'allele_name' column of full-length allele names or the 'trimmed_allele' column of two-field allele names in the pertinent genomic alignment. Returns FALSE if the allele name is not properly formed, or if the allele name is not found in HLAalignments.
#'
#'@param allele A character string of the colon-delimited HLA allele name.
#'
#'@return A logical identifying if the allele name is present in the alignments (TRUE) or, if it is not in the alignments or is not valid not (FALSE).
#'
#'@note Messages will be returned to the console if the allele name is malformed, or the locus is invalid; e.g., validateAllele("C*01:01:01:01") or validateAllele("A*01:01:01:01:01").
#'@note The locus being evaluated must be included in HLAalignments.
#'
#'@export
#'
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
#'Determine if an Allele Name Ever Existed, and (if so) its Most Recent IPD-IMGT/HLA Database Release
#'
#'Returns TRUE if an allele name is present in AlleleListHistory or FALSE it is absent, or c(TRUE,version), where 'version' is the most recent IPD-IMGT/HLA Database release in which that name appeared, when version = TRUE.
#'
#'@param allele A character string of an HLA allele name. Colon-delimited and field-delimited names are both accepted.
#'@param version A logical that indicates if the most recent nomenclature release version in which that name was valid should be returned. 
#'
#'@return A logical identifying if the allele name is found in AlleleListHistory (TRUE) or not (FALSE), or c(TRUE,version) if version = TRUE.
#'
#'@export
#'
#'@examples
#'verifyAllele("A*01:01:01:01")
#'verifyAllele("A*01:01:01:01",TRUE)
#'verifyAllele("A*010101",TRUE)
#'verifyAllele("A*0101",TRUE)
#'
verifyAllele <- function(allele, version=FALSE){

    resArray <- which(alleleListHistory$AlleleListHistory == allele,arr.ind = TRUE)
    if(length(resArray)==0) {return(FALSE)}
      
    if(!version) {return(TRUE)
        } else {
          rawVersion <- colnames(alleleListHistory$AlleleListHistory)[resArray[1,2][[1]]]
            return(c(TRUE,expandVersion(substr(rawVersion,start=2,stop = nchar(rawVersion)))))
          } 
  }

#################
##parseAlignmentHead
#'Guides For Parsing the Header Blocks of Alignment Files
#'
#'Returns a vector describing the location of key information in the header blocks of alignment files.
#'
#'@param version A character string of a validated IPD-IMGT/HLA Database release version (e.g., '3.25.0' or '3250').
#'
#'@return Either FALSE if the version is not valid, or two-value numerical vector describing (1) the header line in which alignment version data appears and (2) the length of the character string in that line preceding version data.
#'
#'@export
#'
#'@examples
#'parseAlignmentHead("3.25.0")
#'
parseAlignmentHead <- function(version){
  
  if(validateVersion(version) ) { version <- squashVersion(version)
      } else {
          return(FALSE) }
  
    field1 <- as.numeric(substr(version,1,1))
        if(field1 != 3) { 
                  message(paste("Alignments are not available for version",field1,"releases.",sep=" "))
                     return(FALSE) 
              }
    field2 <- as.numeric(substr(version,2,3))
    
    if(field2 >= 32) {return(c(3,24))} # on line three, skip the first 24 characters 
    
    if(field2 %in% c(25,27:31)) {return(c(2,22))} # on line 2, skip the first 22 characters
    
    if(field2 %in% c(0:24,26)) {return(c(2,18))} # on line 2, skip the first 18 characters
    
    return(FALSE) ## In case something else makes it through.
    
#    field3 <- substr(version,4,4) ## not needed.
  
}

##################
##queryRelease
#'Search Allele Names Across Release Versions 
#'
#'Searches specific release versions in the AlleleListHistory object for user-defined allele variants. 
#'
#'@param rel An IPD-IMGT/HLA Database release version, represented as either a character (e.g. "3.56.0") or a numeric (e.g., 3560) value.
#'
#'@param variant A character string. The value of 'var' can be any part of a locus or allele name (e.g., "DR", "02:01", "DRB1*08:07"). The default ("") specifies all alleles in 'rel'.
#'
#'@param all A logical. When 'all' = TRUE, a vector of all instances of 'variant' in 'rel' is returned.  When 'all' = FALSE, the number of instances of 'var' in 'rel' is returned. 
#
#'@return A character vector of all matches to 'variant' in 'rel' or the number of all matches to 'variant' in 'rel'.
#'
#'@export
#'
#'@examples
#' # Identify the number of DRB9 alleles in releases 3.30.0 and 3.31.0.
#'queryRelease("3.30.0","DRB9",FALSE) 
#'queryRelease("3.31.0","DRB9",FALSE)
#'
#' # Identify the total number of alleles in release 3.56.0.
#' queryRelease(3560)
#'
queryRelease <- function(rel, variant="", all= FALSE){
  
  if(validateVersion(rel)) {rel <- paste("X",squashVersion(rel),sep="")} else {return(paste("Release version",ifelse(is.numeric(rel),expandVersion(rel),rel),"is not in the local copy of 'AlleleListHisory'.",sep=" "))}
  
      alleles <- alleleListHistory$AlleleListHistory[,rel][!is.na(alleleListHistory$AlleleListHistory[,rel][])]
  
      if(variant == "") {
          
        matchAllele <- alleles 
          
            } else { 
            
            matchAllele <- grep(variant,alleles,fixed=TRUE)
            
              }
  
      numAllele <- length(matchAllele)
  
          if(all) { 
            
                  if(variant == "") {
                          return(matchAllele) 
                    } else {
                   return(alleles[matchAllele])
                    }
    
                 } else { return(numAllele)
    
             }
  }

##################
##checkAlignType
#'Ensure that AlignType Values are Valid
#'
#'Evaluates 'alignType' values to ensure that only "prot", "nuc", "codon" and "gen" are passed to downstream functions. If any other values are entered, a message describing excluded values is generated. If no valid 'alignType' values are are present, an error is generated, and any calling function is ended. 
#'
#'@param alignType A vector of character values specifying sequence alignment types to be used for a function.
#'
#'@return A character vector that includes only allowed 'alignType' values.
#'
#'@export
#'
#'@examples
#'checkAlignType(c("nuc","prot","gDNA")) 
#'
checkAlignType <- function(alignType) {
  
  if(!all(alignType %in% c("codon","nuc","prot","gen"))) {
    
    fixedType <- alignType[alignType %in% c("codon","nuc","prot","gen")]
    
    if(length(fixedType) == 0 ){stop(paste("None of the values",paste(alignType,collapse=", "),"are valid 'alignType' values."))}
    
    message(paste("The value '",alignType[!alignType %in% fixedType],"' was removed from 'alignType', as it is not a valid 'alignType' value.\n",sep=""))
    
    return(fixedType)
    }
  
  alignType
}

##################
##checkSource
#'Ensure that Source Values are Valid
#'
#'Evaluates 'source' values to ensure that only "AA", "cDNA", and "gDNA" are passed to downstream functions. If any other values are entered, a message describing excluded values is generated. 
#'
#'@param source A vector of character values specifying sequence alignment file sources to be used for a function.
#'
#'@return A character vector that includes only allowed 'source' values.
#'
#'@export
#'
#'@examples
#'checkSource(c("AA","cDNA","codon")) 
#'
checkSource <- function(source) {
  
  if(!all(source %in% c("AA","cDNA","gDNA"))) {
    
    fixedSource <- source[source %in% c("AA","cDNA","gDNA")]
    
    if(length(fixedSource) == 0 ){stop(paste("None of the values",paste(source,collapse=", "),"are valid 'source' values."))}
    
    message(paste("The value '",source[!source %in% fixedSource],"' was removed from 'source', as it is not a valid 'source' value.\n",sep=""))
    
    return(fixedSource)
  }
  
  source
}

##################
##TypeToSource
#'Convert AlignType Values to Source Values
#'
#'Converts between 'alignType' values, identifying four types of sequence alignments and 'source' values, identifying three kinds of alignment files.
#'
#'@param alignVector A vector of character values specifying kinds of four sequence alignment types ("prot","nuc","codon" and "gen"), or three source alignment files ("AA","cDNA", and "gDNA").
#'
#'@param toSource A logical. If 'toSource' is true, 'alignType' values are converted to 'source' values. If 'toSource' is false, 'source' values are concerted to 'alignType' values
#'
#'@return A character string of the converted 'align' or 'source' value.  
#'
#'@export
#'
#'@examples
#'typeToSource(c("nuc","prot","gen"),TRUE) 
#'
typeToSource <- function(alignVector,toSource=TRUE){
  
  if(toSource) {
    
    alignVector <- checkAlignType(alignVector)
    
    alignVector[alignVector == "prot"] <- "AA"
    alignVector[alignVector == "nuc"] <- "cDNA"
    alignVector[alignVector == "codon"] <- "cDNA"
    alignVector[alignVector == "gen"] <- "gDNA"
    
    unique(alignVector)
  } else {
    
      alignVector <- checkSource(alignVector)
    
      alignVector[alignVector == "AA"] <- "prot"
      alignVector[alignVector == "cDNA"] <- "nuc"
      alignVector[alignVector == "gDNA"] <- "gen"
    
      if("nuc" %in% alignVector) {alignVector <- append(alignVector,"codon",length(alignVector))} ## cDNA is both nuc and codon
    
      unique(alignVector)
      
      }  
}

#################
##AddCodonLine
#'Add an 'AA codon' Line to Alignments When Missing.
#'
#'Modifies cDNA alignment objects that are missing "AA codon" lines to include these lines in the correct location with the correct codon position information.
#'
#'@param cDNAalign A matrix of cDNA alignment lines, generated from an alignment file that is missing "AA codon" lines.
#'@param firstPos A numeric value identifying the position number of the transcript's N-terminal codon, based on a complete cDNA alignment in another release.
#'@param afterLine A numeric value identifying the number of the line below the first "cDNA" line in the alignment. The default value is 8.
#'@param codons A numeric value identifying the number of codons in each line of the nucleotide alignment. The default value is 25.
#'
#'@return A complete cDNA alignment data frame that includes "AA codon" rows.
#'
#'@export
#'
#'@note For internal HLAtools use.
#'
addCodonLine <- function(cDNAalign,firstPos,afterLine = 8, codons = 25){
  
  # Assigns the AA codon positions using the first position in the codon row of another release, and adds the missing lines.
  afterRow <- rev(which(cDNAalign == cDNAalign[afterLine])-1)
  numRow <- length(afterRow)
  startPos <- rep(firstPos,numRow) ## this needs to be set based on a correct alignment 
  for(k in 2:length(startPos)) {startPos[k] <- startPos[k-1]+codons}
  if(any(startPos < 1)) {startPos[startPos > 0] = startPos[startPos > 0 ]+1 }
  startPos <- rev(startPos)
  
  for(f in 1:length(startPos)) {
    cDNAalign <- append(cDNAalign,paste(" AA codon          ",startPos[f],sep=""),after=afterRow[f])
  }      
  
  cDNAalign
}
