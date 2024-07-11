## Locus Validation Functions v1.0.0 17APR2024

################
##validateLocus
#'Determine if a Locus Name is in the HLAgazeteer
#'
#'Checks a vector of HLA locus names against the HLA gazeteer to determine if the locus name is valid for a specific type of alignment.
#'
#'@param loci A character vector of HLA gene names (ex. "DRB1", c("DRB1","DQB1")).
#'@param source A character vector of alignment source types. "AA", "cDNA", and "gDNA" are allowed types.
#'
#'@return A logical value. TRUE indicates that all of the names and source types are valid. FALSE indicates that at least one locus name or alignment source type is invalid.
#'
#'@export
#'
#'@note The results of this check should only be considered valid for the IPD-IMGT/HLA Database release version of the current gazeteer.
#'
#'@examples
#'validateLocus(loci = "DRB1", source = "AA")
#'validateLocus(loci = c("V"), source = c("cDNA","gDNA"))
#'validateLocus(loci = c("E","F","G"), source = "gDNA")
validateLocus<-function(loci, source){
      valid <- TRUE

      source <- checkSource(source)

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
#'Apply validateLocus() to Multiple Loci
#'
#'Applies validateLocus() to a vector of locus names, validates them against HLAgazeteer$gen, and returns a vector of validated locus names. 
#'
#'@param loci A character vector of locus names found in the current HLAgazeteer. 
#'
#'@param source A character vector describing the type of alignment 'loci' should be validated against. Allowed options are 'AA' (protein), 'cDNA' (nucleotide) and 'gDNA' (genomic). Only a single value should be provided. The default value is 'gen'.
#'
#'@param verbose A logical value. If 'verbose' = TRUE, messages describing invalid 'loci' or 'source' values are generated. If 'verbose' = FALSE, no messages are generated.
#'
#'@return A character vector of locus names that are present in HLAgazeteer$gen.
#'
#'@export
#'
#'@note The results of this check should only be considered valid for the IPD-IMGT/HLA Database release version of the current HLAgazeteer.
#'
#'@examples
#'multiLocusValidation(loci = c("DRB1","DPB1","DQB8"))
#'multiLocusValidation(loci = c("A","B","C","D","Q"))
#'
multiLocusValidation <- function(loci, source = "gDNA", verbose = TRUE) {
    lociTest <- rep(FALSE,length(loci))
    
    funCaller <- as.character(sys.call(sys.parent())) ## checking to see if alignmentFull() called mLV(), because its source is always "gDNA", but user may specify "DRB" as a locus.

    source <- checkSource(source)
    if(length(source) != 1) {stop(if(verbose){paste("Please provide only a single 'source' value.\n")})}
    
    subSource <- source ### record the submitted source value for subverting 'source=gDNA' for cases when 'loci=DRB' when called from alignmentFull() and 'loci != all'.
    
      for(i in 1:length(loci)) {
            if(!loci[i] == "all") { 

                                    if(funCaller[1] == "alignmentFull" && loci[i] == "DRB") { 
                                            source <- "cDNA"} ### change source to "cDNA" for "DRB" only
              
                if(suppressMessages(validateLocus(loci[i],source))) { 
                   
                      lociTest[i] <- TRUE 
                      } 
                                    if(subSource != source ) { source <- subSource} ### return source to its original submitted version
                }
              }
        if(any(lociTest == FALSE)) {if(verbose) {
            message(paste("The",loci[lociTest == FALSE],"locus is invalid in version", HLAgazeteer$version ,"and has been removed.\n",sep=" "))
              }
         }
    
      loci <- loci[lociTest == TRUE]

    loci
}

