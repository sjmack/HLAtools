### Tools for Validating IPD-IMGT/HLA Release Version Values v2.2.0 8 March 2024

################
##CheckVersion
#'CheckVersion - Determine if an IPD-IMGT/HLA Release Version is supported in the local HLAtools package.
#'
#'Determines if data for a given IPD-IMGT/HLA Release Version is present in the local HLAtools package.
#'
#'@param version A release version (branch) of the ANHIG/IMGTHLA Github repository (e.g. '3.53.0'). The value 'Latest' refers to the most recent release.
#'
#'@return A logical. TRUE indicates that data for 'version' is available to the local HLAtools package. FALSE means that data for 'version' is not available.
#'
#'@export
#'
#'@examples
#'\dontrun{
#'checkVersion(version = "3.25.0")
#'checkVersion(version = "Latest")
#'}
#'
#'@note
#'For internal HLAtools use.
checkVersion <- function(version){
  if(version == "Latest") { version <- getLatestVersion()}

  paste("X",gsub(".","",version,fixed=TRUE),sep="") %in% colnames(alleleListHistory$AlleleListHistory)

}

################
##ValidateVersion
#'ValidateVersion - Validate a potential IPD-IMGT/HLA Release Version
#'
#'Determines if data for a given IPD-IMGT/HLA Release Version is available in the ANHIG/IMGTHLA GitHub repository.
#'
#'@param version A release version (branch) of the ANHIG/IMGTHLA Github repository (e.g. '3.53.0'). The value 'Latest' refers to the most recent release.
#'
#'@return A logical. TRUE indicates that data for 'version' is available in the ANHIG/IMGTHLA GitHub repository. FALSE means that data for 'version' is not available.
#'
#'@export
#'
#'@examples
#'\dontrun{
#'validateVersion(version = "3.31.1")
#'validateVersion(version = "3.13.0")
#'}
#'
#'@note
#'For internal HLAtools use.
validateVersion <- function(version) {
  if(!checkVersion(version)) {

    releases <- strsplit(readLines(url("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist_history.txt"),n=7,ok=TRUE,skipNul = FALSE)[7],",")[[1]]
    releases <- releases[2:length(releases)]

    squashVersion(version) %in% releases

  } else {TRUE}

}

################
##GetLatestVersion
#'GetLatestVersion - Identifies the latest IPD-IMGT/HLA Database Release
#'
#'Identifies the most recent version of the IPD-IMGT/HLA Database available in the ANHIG/IMGTHLA Github repository.
#'
#'@return A dot-delimited character string identifying the latest release version (branch) of the ANHIG/IMGTHLA Github repository.
#'
#'@export
#'
#'@examples
#'\dontrun{
#'getLatestVersion()
#'}
#'
#'@note
#'For internal HLAtools use.
getLatestVersion <- function() {
  strsplit(readLines(url("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/release_version.txt"),n=3,ok=TRUE,skipNul = FALSE)[3]," ")[[1]][4]
}

################
##SquashVersion
#'SquashVersion -- Reduce an IPD-IMGT/HLA Database Release Version to a sequence of numbers
#'
#'Removes the '.' delimiters from an IPD-IMGT/HLA Database Release Version name
#'
#'@param ver A character vector-formatted IPD-IMGT/HLA Database Release Version name (e.g., '3.54.0').
#'@param num A Logial. A numeric value is returned when num = TRUE. A character vector is returned when num = FALSE. The default value is FALSE.
#'
#'@return A 4-digit value, as either a character string or a number
#'
#'@export
#'
#'@examples
#'squashVersion("3.45.0",TRUE)
#'
#'@note
#'For internal HLAtools use.
squashVersion <- function(ver,num = FALSE){

  ver <- gsub(".","",ver,fixed=TRUE)
  if(num){as.numeric(ver)} else {ver}
}

################
##ExpandVersion
#'ExpandVersion -- Adds 'dot' delimiters to an 'all numeral' IPD-IMGT/HLA Database Release Version
#'
#'Adds the '.' delimiters to a 'squashed' IPD-IMGT/HLA Database Release Version name
#'
#'@param ver A four character/numeral IPD-IMGT/HLA Database Release Version name (e.g., '3540' or 3450).
#'
#'@return A charachter vector of the expanded IPD-IMGT/HLA Database Release Version name
#'
#'@export
#'
#'@examples
#'expandVersion(3450)
#'
#'@note
#'For internal HLAtools use.
#'@note
#'This function assumes single-character-length first and third fields of IPD-IMGT/HLA Database Release Version names.
expandVersion <- function(ver){

  if(is.numeric(ver)) {ver <- as.character(ver)}
  ver <- strsplit(ver,split="")[[1]]
  paste(ver[1],paste(ver[2:(length(ver)-1)],collapse = ""),ver[length(ver)],sep=".")
}

