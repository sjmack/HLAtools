#alleleTrim v1.1.1 1/29/2024 -- Steven J. Mack

################
##alleleTrim
#'Trims an HLA allele name to a specified number of fields.
#'
#'@param allele A full HLA allele name formatted as locus*allele_name or HLA-locus*allele_name.
#'@param resolution A number identifying the number of fields to trim the allele down to.
#'@param version the HLA Nomenclature version under which the allele was named. Version 1 allele names are found in IPD-IMGT/HLA Database releases 1.0.0 to 1.16.0. Version 2 allele names are found in IPD-IMGT/HLA Database releases 2.0.0 to 2.28.0. Version 3 allele names are found in IPD_IMGT/HLA Database releases 3.0.0 and onward.
#'
#'@return a trimmed allele name, shortented according to the input parameters.
#'
#'@export
#'
#'@examples
#'alleleTrim(allele = "A*03:01:01", resolution = 2)
#'
alleleTrim <- function(allele,resolution,version=3){

  if(!resolution %in% c(1,2,3,4)) {
    message("Resolution must be between 1 and 4.")
    return(allele)
  }

  if(version == 3){
  currRes <- lengths(regmatches(allele,gregexpr(":",allele)))+1

  if(currRes > resolution) {allele <- GetField(allele,resolution)}

  }

  ## For versions 1 and 2, split out the allele prefix and suffix

  if(version %in% c(1,2)) {
    suffix <- strsplit(allele,"*",fixed=TRUE)[[1]][2]
    prefix <- strsplit(allele,"*",fixed=TRUE)[[1]][1]
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
       return(alleleTrim(allele,resolution,version=2))
        } ## This is for the last version 1 release (1.16.0), which frustratingly used 2 digits for each field

     if(suffixLen < 8 && substr(suffix,nchar(suffix),nchar(suffix)) %in% LETTERS) {currRes <- currRes-1}

     if(currRes > resolution) {
       if(resolution %in% c(1,2) ) {allele <- paste(prefix,substr(suffix,1,resolution*2),sep="*")}
       if(resolution == 3) {allele <- paste(prefix,substr(suffix,1,5),sep="*")}

     }
  }
  allele
}

###############
##GetField
#' An HLA trimming function
#'
#' Trim a properly formatted HLA allele to a desired number of fields.
#'
#' @param x HLA allele.
#' @param Res Resolution desired.
#'
#' @note for internal use only.
#'
#' @export
#'
#'@examples
#'GetField("HLA-A*01:01:01:01", 3)
#'
GetField <- function(x,Res) {
  Tmp <- unlist(strsplit(as.character(x),":"))
  if (length(Tmp)<2) {
    return(x)
  } else if (Res==1) {
    return(Tmp[1])
  } else if (Res > 1) {
    Out <- paste(Tmp[1:Res],collapse=":")
    return(Out)
  }
}
