### buildIMGTAlleleStatus v.1.0.0 25JUL2026

####################
##updateAlleleStatus
#
#'Build the alleleStatus object 
#'
#'Consumes the ANHIG/IMGTHLA GitHub repository's Allele_status.txt file to create the AlleleStatus object.
#'
#'@param version A string identifying the IPD/IMGT-HLA Database release version of the desired Allele_status.txt file. 
#'
#'@return A list-object consisting $AlleleStatus, a data.frame version of the Allele_status.txt file, and $version, a character-string identifying the pertinent IPD-IMGT/HLA Database release version. If an invalid version value is provided, NA will be returned.
#'
#'@note The Allele_status.txt file was introduced in IPD/IMGT-HLA Database release version 3.25.0 and is not available in earlier releases.
##'
#'@export
#' 
updateAlleleStatus <- function(version = "Latest"){
  
            if(!validateVersion(version)) {
                  message(paste("Verson",version,"does not exist.",sep=" "))
                  return(NA)
            }
  
            ## Allele_status.txt was first introduced in release veersion 3.25.0.
            if(version %in% c("3.24.0","3.23.0","3.22.0","3.21.0","3.20.0","3.19.0","3.18.0","3.17.0","3.16.0",
                              "3.15.0","3.14.0","3.13.1","3.12.0","3.11.0","3.10.0","3.09.0","3.08.0","3.07.0",
                              "3.06.0","3.05.0","3.04.0","3.03.0","3.02.0","3.01.0","3.00.0")) {
             message(paste("The Allele_status.txt file does not exist in release verson ",version,".",sep=""))
             return(NA)
                }

            AlleleStatusURL <- paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/refs/heads/",squashVersion(version),"/Allele_status.txt",sep="")
  
            alleleStatusTable <- read.table(AlleleStatusURL,header = T,sep = ",",quote = "",colClasses = "character")
  
            if(version == "Latest"){version <- getLatestVersion()}
            
            alleleStatus  <- list(alleleStatusTable,version)
            names(alleleStatus) <- c("AlleleStatus","version")
    
            alleleStatus
}
