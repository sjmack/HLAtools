## v1.2.0 15 March 2024

####################
##checkgDNAStart
#'Identify gDNA Alignments in Which the First Feature Boundary is not Identified as Position +1.
#'
#'Checks the position of the first feature boundary (the 5' UTR - Exon 1 boundary in expressed genes) in each gDNA alignment, and identifies those alignments in which that position is not 1.
#'
#'This function reviews the gDNA atlases in the HLAatlas object and returns a data frame of the first feature boundary position for each locus. For expressed genes and some pseudogenes, this is the position of the start of Exon 1.
#'
#'@param verbose A logical indicating if loci with first feature boundary positions that are not 1 should be reported in the console (verbose = TRUE).
#'
#'@return A one-row data frame with one column for each locus with a gDNA alignment.
#'
#'@note This function requires that the HLAalignments object has been populated with alignments via the alignmentFull() function.
#'
#'@export
#'
#'@note For internal HLAtools use.
#'
checkgDNAstart <- function(verbose=FALSE){

  gLen <- length(HLAatlas$gen)
  findings <- matrix(rep(1,gLen),ncol=gLen)
  colnames(findings) <- names(HLAatlas$gen)

    tic <- 0
      for (i in 1:gLen){

        if(HLAatlas$gen[[i]][1,1] != 1) {

        if(verbose) {cat(paste("The protein start position for the ",names(HLAatlas$gen)[i]," locus is ",HLAatlas$gen[[i]][1,1],".","\n",sep=""))}
      findings[i] <- as.numeric(HLAatlas$gen[[i]][1,1])
    tic <- tic + 1
    }

        if(tic != 0 && verbose) {paste("All loci start with position 1.")}
        }
      as.data.frame(findings)
    }
