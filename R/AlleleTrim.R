#alleleTrim v2.1.0 8/20/2024 -- Steven J. Mack

################
##alleleTrim
#'Trim All Versions of HLA Allele Names
#'
#'Trims an HLA allele name to a specified number of fields or number of digits, depending on the nomenclature epoch. Expression variant suffixes in full-length allele names can be appended to truncated allele names.
#'
#'@param allele A character string of an HLA allele name formatted as locus*allele_name, optionally including the "HLA-" prefix.
#'@param resolution A number identifying the number of fields to trim the allele name to.
#'@param version A number identifying the HLA nomenclature epoch under which the allele was named. Epoch 1 allele names are found in IPD-IMGT/HLA Database releases 1.0.0 to 1.16.0. Epoch 2 allele names are found in IPD-IMGT/HLA Database releases 2.0.0 to 2.28.0. Epoch 3 allele names are found in IPD_IMGT/HLA Database releases 3.0.0 and onward.
#'@param append A logical. When append = TRUE, the expression variant suffix of a full-length allele name is appended to a truncated allele name. The default value is FALSE.
#'
#'@return A character string of the trimmed allele name, shortened according to the input parameters.
#'
#'@note Expression variant suffixes will not be removed from full-length allele names that have fewer than four fields or eight digits.
#'
#'@export
#'
#'@examples
#'alleleTrim(allele = "A*03:01:01", resolution = 2)
#'alleleTrim(allele = "A*030101", resolution = 2,version = 2)
#'alleleTrim(allele = "A*0303N",resolution = 1,version = 1, append = TRUE)
#'alleleTrim("HLA-A*24020102L",resolution = 3,version = 2, append = TRUE)
#'
alleleTrim <- function(allele,resolution,version=3,append = FALSE){

  if(!resolution %in% c(1,2,3,4)) {
    message("Resolution must range from 1 and 4.")
    return(allele)
  }

  if(!version %in% c(1,2,3)) { 
    message("Version must range from 1 and 3.")
    return(allele)
    }
  
  if(version == 3){
  currRes <- lengths(regmatches(allele,gregexpr(":",allele)))+1

  if(currRes > resolution) {allele <- getField(allele,resolution,append)}

  }

  ## For versions 1 and 2, split out the allele prefix and suffix

  if(version %in% c(1,2)) {
    suffix <- strsplit(allele,"*",fixed=TRUE)[[1]][2]
    prefix <- strsplit(allele,"*",fixed=TRUE)[[1]][1]
    
      expVar <- ""
      if(append){
         if(substr(allele,nchar(allele),nchar(allele)) %in% c("N","L","S","C","A","Q")) {
                expVar <- substr(allele,nchar(allele),nchar(allele))
             }
       }
  }

  if(version == 2){
    currRes <- nchar(strsplit(allele,"*",fixed=TRUE)[[1]][2])/2
    if(round(currRes,0) > resolution) {allele <- paste(prefix,substr(suffix,1,2*resolution),sep="*")}
  }
  if(version == 1){
    suffixLen <- nchar(suffix)
     if(suffixLen == 4) {currRes <- 2}
     if(suffixLen == 5 || suffixLen == 6) {currRes <- 3} ## Here and in the next line, the second option is due to the presence of an expression variant (e.g. "N")
     if(suffixLen == 7 || suffixLen == 8) {currRes <- 4}
     if(suffixLen == 9) {
       return(alleleTrim(allele,resolution,version=2, append=append))
        } ## This is for the last version 1 release (1.16.0), which frustratingly used 2 digits for each field

     if(suffixLen < 8 && substr(suffix,nchar(suffix),nchar(suffix)) %in% LETTERS) {currRes <- currRes-1}

     if(currRes > resolution) {
       if(resolution %in% c(1,2) ) {allele <- paste(prefix,substr(suffix,1,resolution*2),sep="*")}
       if(resolution == 3) {allele <- paste(prefix,substr(suffix,1,5),sep="*")}

     }
  }
  
  if(version %in% c(1,2) && append && expVar != "") {
    
      if(!substr(allele,nchar(allele),nchar(allele)) %in% c("N","L","S","C","A","Q")) {
        allele <- paste(allele,expVar,sep="")
            }
      }
  
  allele
}

###############
##getField
#' Trim Colon-Delimited HLA Allele Names by Field
#' 
#'@description
#' Trims a properly formatted colon-delimited HLA allele name to a desired number of fields.
#' 
#' If an allele name with an expression-variant suffix is truncated, the suffix can be appended to the end of the truncated allele name. If a resolution value greater then the number of fields in the submitted field is specified, the original allele is returned.
#'
#' @param allele A character string representing an HLA allele.
#' @param res A numeric describing the resolution desired.
#' @param append A logical. When append = TRUE, the expression variant suffix of a full-length allele name is appended to a truncated allele name. The default value is FALSE.
#'
#' @note For internal HLAtools use.
#' 
#' @return A version of the 'allele' character string that has been trimmed to 'res' resolution.
#'
#' @export
#'
#'@examples
#'getField("HLA-A*01:01:01:01", 3)
#'getField("DRB1*11:01:01:12N", 2,TRUE)
#'
getField <- function(allele,res,append=FALSE) {
  Tmp <- unlist(strsplit(as.character(allele),":"))
  suffix <- ""
  
  if(length(Tmp) < res) return(allele) # added in case the requested resolution is larger than the allele's resolution
  
  if(substr(Tmp[length(Tmp)],nchar(Tmp[length(Tmp)]),nchar(Tmp[length(Tmp)])) %in% c("N","L","S","C","A","Q")) {suffix <- substr(Tmp[length(Tmp)],nchar(Tmp[length(Tmp)]),nchar(Tmp[length(Tmp)])) }
  
  if (length(Tmp)<2) {
    return(allele)
  } else if (res==1) {
    return(paste(Tmp[1],ifelse(append,suffix,""),sep=""))
  } else if (res > 1) {
    Out <- paste(paste(Tmp[1:res],collapse=":"),ifelse(append,suffix,""),sep="")
    return(Out)
  }
}

###############
##multiAlleleTrim
#' Trim Multiple HLA Allele Names
#' 
#'Trim a vector of allele names to a specified number of fields or digits.
#'
#'@param alleles A vector of HLA allele names from a single nomenclature epoch.
#'@param resolution A number identifying the number of fields to which the alleles should be trimmed.
#'@param version A number identifying the HLA nomenclature epoch under which the allele was named. Epoch 1 allele names are found in IPD-IMGT/HLA Database releases 1.0.0 to 1.16.0. Epoch 2 allele names are found in IPD-IMGT/HLA Database releases 2.0.0 to 2.28.0. Epoch 3 allele names are found in IPD_IMGT/HLA Database releases 3.0.0 and onward.
#'@param append A logical. When append = TRUE, the expression variant suffix of a full-length allele name is appended to a truncated allele name. The default value is FALSE.
#'
#'@return A vector of the trimmed allele name, shortened according to the input parameters.
#'
#'@note Expression variant suffixes will not be removed from full-length allele names that have fewer than four fields or eight digits.
#'
#'@export
#'
#'@examples
#'alleles <- c("A*02:01:01:02L","DRB1*08:07", "DQB1*04:02:01:16Q")
#'multiAlleleTrim(alleles,2)
#'multiAlleleTrim(alleles,2,append = TRUE)
#'
#'alleles <- c("A*01010102L","DRB1*1613N","HLA-Cw*0322Q")
#'multiAlleleTrim(alleles,1,2,TRUE)
#'multiAlleleTrim(alleles,2,2)
#'
multiAlleleTrim <- function(alleles,resolution,version=3, append = FALSE){
 
  intTrim <- Vectorize(HLAtools::alleleTrim)
                       
  unname(intTrim(alleles,resolution,version,append))
               
}
