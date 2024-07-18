##atlasFull v1.2.0 - 16 April 2024

################
##atlasFull
#'Generate a Complete set of Protein, Coding nucleotide and Genomic Nucleotide Atlases
#'
#'Applies atlasMaker() to build a set of 'atlases' identifying the gene-feature (exon, intron, and UTR) boundaries in the protein, coding nucleotide, and genomic nucleotide alignments of all HLA loci.
#'
#'@param version A character string identifying the version of the ANHIG/IMGTHLA Github repository used to generate atlases. The default value, "Latest", generates an atlas for the most recent IPD-IMGT/HLA Database release.
#'
#'@return A list object containing a list of data frames of protein (prot), coding nucleotide (nuc), and genomic nucleotide (gen) atlases for all genes in the IPD-IMGT/HLA Database release, along with a character vector identifying the pertinent IPD-IMGT/HLA Database release.
#'
#'@note Data informing the atlases were downloaded from the ANHIG/IMGTHLA Github repository.
#'
#'@note
#'For internal HLAtools use.
#'
#'@export
#'
atlasFull <- function(version = "Latest") {

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
    cList[i] <- suppressWarnings(atlasMaker(as.character(NL1[i]), "cDNA", version = version)[[1]][1])
  }
  for(i in 1:length(NL2)){
    gList[i] <- suppressWarnings(atlasMaker(as.character(NL2[i]), "gDNA", version = version)[[1]][1])
  }
  for(i in 1:length(NL3)){
    protList[i] <- suppressWarnings(atlasMaker(as.character(NL3[i]), "AA", version = version)[[1]][1])
  }

  #naming inside of nested lists
  names(cList) <- NL1
  names(gList) <- NL2
  names(protList) <- NL3

  #placing nested lists inside of larger list
  HLAatlas<-list(firstList <- protList,
                      secondList <- cList,
                      thirdlist <- gList,
                      versionString <- version)
                    # versionString <- ifelse(version == "Latest",getLatestVersion(),version))
  #naming nested lists
  names(HLAatlas) <- c("prot","nuc","gen","version")

  HLAatlas

}
