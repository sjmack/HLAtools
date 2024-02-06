## UpdateAll v01.1.1 5 February 2024

################
##updateAll
#'Updates all of the data objects derived from IPD-IMGT/HLA sources.
#'
#'Applies extractGeneTypes(), buildGazeteer(), ffN(), atlasFull() and alignmentsFull() to build a set of alignments for all loci supported in the ANHIG/IMGTHLA GitHub repository.
#'
#'@param version The version of the ANHIG/IMGTHLA Github repository to build alignments from. The default value is to call the getLatestVersion() function to generate alignments for the most recent IPD-IMGT/HLA Database release.
#'
#'@return A list object containing data frames of protein (prot), coding nucleotide (nuc), and genomic nucleotide (gen) alignments for all genes in the specified IPD-IMGT/HLA Database release, along with the pertinent reference database release.
#'
#'
#'@examples
#' updateAll()
#'
#'@note
#'For internal HLAtools use.
#'@export
updateAll <-function(version = getLatestVersion()){

  #build new local objects

  IMGTHLAGeneTypes <- buildIMGTHLAGeneTypes() # this function takes no version argument.
  HLAgazeteer <- buildGazeteer(version)
  alleleListHistory <- updateAlleleListHistory()
  fragmentFeatureNames <- ffN(version)
  HLAatlas <- atlasFull(version)
  HLAalignments <- alignmentFull(version)

}
