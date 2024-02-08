### Alignment Full v1.2.0 1 February 2024

################
##alignmentFull
#'Build a complete set protein, coding nucleotide and genomic nucleotide alignments for all IPD-IMGT/HLA Loci.
#'
#'Applies buildAlignments() to build a set of alignments for all loci supported in the ANHIG/IMGTHLA GitHub repository.
#'
#'@param version The version of the ANHIG/IMGTHLA Github repository to build alignments from. The default value, "Latest", generates alignments for the most recent IPD-IMGT/HLA Database release.
#'
#'@return A list object containing data frames of protein (prot), coding nucleotide (nuc), and genomic nucleotide (gen) alignments for all genes in the specified IPD-IMGT/HLA Database release, along with the pertinent reference database release.
#'
#'
#'@examples
#'\dontrun{
#' HLAalignments <- alignmentFull()
#'}
#'
#'@note
#'For internal HLAtools use.
#'@export

alignmentFull <- function(version = "Latest") {

  if(version != "Latest"){ #
    if(!validateVersion(version)){stop(paste(version," is not a valid IPD-IMGT/HLA Database release version."))}
  }else{ version <- getLatestVersion()}

  #list of loci names
  NL1 <- as.list(HLAgazeteer$nuc) #cDNA
  NL2 <- as.list(HLAgazeteer$gen) #gDNA
  NL3 <- as.list(HLAgazeteer$prot) #AA

  #creating nested lists for cDNA, gDNA and amino acid
  cList <- vector(mode='list', length=length(NL1))
  gList <- vector(mode='list', length=length(NL2))
  protList <- vector(mode='list', length=length(NL3))

  #filling nested lists with alignments
  for(i in 1:length(NL1)){
    cList[i] <- buildAlignments(as.character(NL1[i]), "cDNA", version = version)[[1]][1]
  }
  for(i in 1:length(NL2)){
    gList[i] <- buildAlignments(as.character(NL2[i]), "gDNA", version = version)[[1]][1]
  }
  for(i in 1:length(NL3)){
    protList[i] <- buildAlignments(as.character(NL3[i]), "AA", version = version)[[1]][1]
  }
  #naming inside of nested lists
  names(cList) <- NL1
  names(gList) <- NL2
  names(protList) <- NL3

  #placing nested lists inside of larger list
  AllAlignment<-list(firstList <- protList,
                 secondList <- cList,
                 thirdlist <- gList,
                 versionString <- version)
#naming nested lists

  names(AllAlignment) <- c("prot","nuc","gen","version")

  AllAlignment
}
