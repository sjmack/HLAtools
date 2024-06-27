#updateAlleleListHistory v1.2.0 2FEB2024

################
##updateAlleleListHistory
#'Build the AlleleList History R Object
#'
#'Consumes the ANHIG/IMGTHLA GitHub repository's allelelisthistory.txt file to create the AlleleListHistory object.
#'
#'@return A list object containing two data frames. alleleListHistory$Version identifies the date and version of alleleListHistory$AlleleListHistory. alleleListHistory$AlleleListHistory relates each IPD-IMGT/HLA Database Accession ID (HLA_ID) to the corresponding HLA allele name in each release version since 1.05.0.
#'
#'@importFrom utils read.table
#'
#'@export
#'
updateAlleleListHistory <- function() {
  url<-url("https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist_history.txt")
  alleleListHistoryHead <- readLines(url,n = 3)
  alleleListHistory <- read.table(url, skip = 6, header=TRUE, sep=",")
  AlleleListHistory <- list(ALversion <- data.frame(strsplit(alleleListHistoryHead[2]," ")[[1]][3],strsplit(alleleListHistoryHead[3]," ")[[1]][4]), ALtable <- alleleListHistory)
  names(AlleleListHistory) <- c("Version","AlleleListHistory")
  colnames(AlleleListHistory$Version) <- c("Date","Version")

  AlleleListHistory
}

