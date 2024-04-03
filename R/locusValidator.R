## Locus Validation Functions v0.3.0 2APR2024

################
##validateLocus
#'Determine if a locus name is present in the HLAgazeteer.
#'
#'Checks a vector of HLA locus names against the HLA gazeteer to determine if the locus name is valid for a specific type of alignment.
#'
#'@param loci A vector of HLA gene names (ex. "DRB1", c("DRB1","DQB1")).
#'@param source A vector of alignment types. "AA", "cDNA", and "gDNA" are allowed types.
#'
#'@return A logical value. TRUE indicates that all of the names and types are valid. FALSE indicates that at least one locus name or alignment type is invalid.
#'
#'@export
#'
#'@note The results of this check should only be considered valid for the IPD-IMGT/HLA Database release version of the current gazeteer.
#'
#'@examples
#'validateLocus(loci = "DRB1", source = "AA")
#'validateLocus(loci = c("V"), source = c("cDNA","gDNA"))
#'validateLocus(loci = c("E","F","G"), source = "aDNA")
validateLocus<-function(loci, source){
valid <- TRUE
for(i in 1:length(source)){
    if(!source[i] %in% c("AA","cDNA","gDNA")) {
      message(paste(source[i], "is not a valid alignment type. Valid types are 'AA', 'cDNA' and 'gDNA'.",paste=" "))
      return(FALSE)
    }
}
for(j in 1:length(loci)){
  if(loci[j]=="DRB1"|loci[j]=="DRB3"|loci[j]=="DRB4"|loci[j]=="DRB5") next
  for(x in 1:length(source)) {
    if(source[x] == "cDNA") {
      if(loci[j]%in% HLAgazeteer$nuc == FALSE) {
        message(paste("The", loci[j], "locus is not present in", source[x],"alignments.",sep=" "))
        valid <- FALSE
      }
    }
    if(source[x] == "gDNA") {
      if(loci[j]%in% HLAgazeteer$gen == FALSE) {
        message(paste("The", loci[j], "locus is not present in", source[x],"alignments.",sep=" "))
        valid <- FALSE
      }
    }
    if(source[x] == "AA") {
      if(loci[j]%in% HLAgazeteer$prot == FALSE) {
        message(paste("The", loci[j], "locus is not present in", source[x],"alignments.",sep=" "))
        valid <- FALSE
        }
      }
    }
  }
valid
}

################
##multiLocusValidation
#'Apply validateLocus() to multiple individual loci
#'
#'Applies validateLocus() to a vector of locus names, validates them against HLAgazeteer$gen, and returns a vector of validated locus names. 
#'
#'@param loci A vector of locus names. 
#'
#'@param verbose A logocal value. If verbose = TRUE, messages describing the invalid names are generated.
#'
#'@return A vector of locus names that are present in HLAgazeteer$gen.
#'
#'@export
#'
#'@note The results of this check should only be considered valid for the IPD-IMGT/HLA Database release version of the current HLAgazeteer.
#'
#'@examples
#'multiLocusValidation(loci = c("DRB1","DPB1","DQB8"))
#'multiLocusValidation(loci = c("A","B","C","D","Q"))
#'
multiLocusValidation <- function(loci) {
    lociTest <- rep(FALSE,length(loci))

      for(i in 1:length(loci)) {
            if(!loci[i] == "all") { 
                if(suppressMessages(validateLocus(loci[i],"gDNA"))) { 
                
                lociTest[i] <- TRUE 
                  } 
                }
              }
        if(any(lociTest == FALSE)) {message(paste("The",loci[lociTest == FALSE],"locus is invalid and has been removed.\n"),sep=" ")}

      loci <- loci[lociTest == TRUE]

    loci
}

