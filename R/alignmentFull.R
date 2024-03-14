### Alignment Full v3.0.0 13 March 2024

################
##alignmentFull
#'Build sets of protein, codon, coding nucleotide and genomic nucleotide alignments for specified  Loci.
#'
#'Applies buildAlignments() to build a set of alignments for loci supported in the ANHIG/IMGTHLA GitHub repository.
#'
#'@param loci A vector of the locus names for which alignments should be built. The default value ("all") generates alignments for all loci. 
#'@param alignType A vector of alignment types. The allowed values are "all", "prot", "codon", "nuc", and "gen", which specify either all available alignments for a given locus or the respective protein, codon, nucleotide and genomic alignments, as determined by the HLAgazeteer.
#'@param version The version of the ANHIG/IMGTHLA Github repository alignments are obtained from. The default value, "Latest", generates alignments for the most recent IPD-IMGT/HLA Database release.
#'
#'@return A list object containing data frames of protein (prot), codons (codon), coding nucleotide (nuc), or genomic nucleotide (gen) alignments, for specified genes in the specified IPD-IMGT/HLA Database release, along with the pertinent reference database release.
#'
#'
#'@examples
#'\dontrun{
#' HLAalignments <- alignmentFull()
#'}
#'
#'@note
#'This function requires internet access. Depending on local download speeds, building all available alignments for all loci can take several minutes.
#'@export

alignmentFull <- function(loci = "all", alignType = "all", version = "Latest") {

  if(version != "Latest"){ #
    if(!validateVersion(version)){stop(paste(version," is not a valid IPD-IMGT/HLA Database release version."))}
  }else{ version <- getLatestVersion()}

  if(length(loci) == 1 && loci == "all" && alignType == "all") {
  
  #list of loci names
  NL4 <- NL1 <- as.list(HLAgazeteer$nuc) #cDNA and codon
  NL2 <- as.list(HLAgazeteer$gen) #gDNA
  NL3 <- as.list(HLAgazeteer$prot) #AA
  
  } else {
  
      if(alignType == "all" || alignType == "nuc")  { NL1 <- as.list(loci[loci %in% as.list(HLAgazeteer$nuc)]) } #nuc 
      if(alignType == "all" || alignType == "gen")  { NL2 <- as.list(loci[loci %in% as.list(HLAgazeteer$gen)]) } #gen 
      if(alignType == "all" || alignType == "prot")  { NL3 <- as.list(loci[loci %in% as.list(HLAgazeteer$prot)]) } # prot
      if(alignType == "all" || alignType == "codon")  { NL4 <- as.list(loci[loci %in% as.list(HLAgazeteer$nuc)]) } # codon
    
      }

  #creating nested lists for cDNA, gDNA and amino acid
  cList <- vector(mode='list', length=length(NL1))
  gList <- vector(mode='list', length=length(NL2))
  protList <- vector(mode='list', length=length(NL3))
  codonList <- vector(mode='list', length=length(NL4))

  #filling nested lists with alignments
  if(length(NL1) > 0) { for(i in 1:length(NL1)){
    cList[i] <- buildAlignments(as.character(NL1[i]), "cDNA", version = version)[[1]][2] } 
    
    names(cList) <- NL1
    
    }
  if(length(NL2) > 0) { for(i in 1:length(NL2)){
    gList[i] <- buildAlignments(as.character(NL2[i]), "gDNA", version = version)[[1]][1] }
    
    names(gList) <- NL2
  }
  if(length(NL3) > 0) { for(i in 1:length(NL3)){
    protList[i] <- buildAlignments(as.character(NL3[i]), "AA", version = version)[[1]][1] }
    
    names(protList) <- NL3
  }
  if(length(NL4) > 0) { for(i in 1:length(NL4)){
    codonList[i] <- buildAlignments(as.character(NL4[i]), "cDNA", version = version)[[1]][1] }
    
    names(codonList) <- NL4
  }
  
  #naming inside of nested lists
  #names(cList) <- NL1
  #names(gList) <- NL2
  #names(protList) <- NL3
  #names(codonList) <- NL4

  #placing nested lists inside of larger list
  AllAlignment<-list(firstList <- protList,
                 fourthlist <- codonList,
                 secondList <- cList,
                 thirdlist <- gList,
                 versionString <- version)
#naming nested lists

  names(AllAlignment) <- c("prot","codon","nuc","gen","version")

  AllAlignment
}
