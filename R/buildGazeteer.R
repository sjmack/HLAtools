#buildGazeteer v01.2.0 7FEB2024

##############
##buildGazeteer
#'Build lists of genes
#'
#'Queries the ANHIG/IMGTHLA GitHub Repository to identify HLA genes that do and do not have amino acid, cDNA, and genomic (gDNA) alignments.
#'
#'@param version A string identifying of the desired IPD-IMGT/HLA Database release version to which the alleles should be updated. The default value is most recent IPD-IMGT/HLA Database release version.
#'
#'@return A list object of vectors identifying those loci in the IPD-IMGT/HLA Database for which amino acid (prot), cDNA (nuc), and gDNA (gen) alignments are available. In addition, subsets of these loci identified as pseudogenes (pseudo) and gene fragments (frag) for which nucleotide and genomic alignments are and are not available (the latter denoted as 'noNuc' or 'noGen') are also identified.
#'
#'@importFrom stringr fixed
#'@importFrom utils read.table
#'@importFrom HLAtools IMGTHLAGeneTypes
#'
#'@note Pseudogenes and gene fragments are defined by the IPD-IMGT/HLA Database at https://hla.alleles.org/genes/index.html, and are described in the HLAtools::IMGTHLAGeneTypes data object.
#'@note This function requires internet access to function.
#'
#'@examples
#'\dontrun{
#'HLAgazeteer <- buildGazeteer("3.34.0")
#'}
#'
#'@export
buildGazeteer <- function(version = getLatestVersion()) {

  URL <- "https://github.com/ANHIG/IMGTHLA/tree/Latest/alignments"
  rough <- suppressWarnings(read.table(URL))
  rough <- strsplit(rough[[1]], split = "[,]|[:]|[\\{]|[\\}]")
  narrow <- gsub("alignments/","",rough[[1]][grep("alignments/",rough[[1]])]) # pulls only loci in the 'alignments' directory
  narrow <- substr(narrow,2,nchar(narrow)-1)

  nucList <- vector("list", 1)
  genList <- vector("list", 1)
  protList <- vector("list", 1)

 for(i in 1:length(narrow)) {
        if (grepl("Class", narrow[i], fixed = TRUE) == FALSE) {
           if (grepl("nuc", narrow[i], fixed = TRUE)) {
          if(!narrow[i] == "DRB_nuc.txt") { #Exclude combined "DRB" alignments
          pure <- strsplit(narrow[i], "_")
          nucList[[1]] <- append(nucList[[1]], fixed(pure[[1]][1]), after = length(nucList[[1]]))
          }
        }
        if (grepl("gen", narrow[i], fixed = TRUE)) {
          pure <- strsplit(narrow[i], "_")
          genList[[1]] <- append(genList[[1]], fixed(pure[[1]][1]), after = length(genList[[1]]))
         }
        if (grepl("prot", narrow[i], fixed = TRUE)) {
          if(!narrow[i] == "DRB_prot.txt") { #Exclude combined "DRB" alignments
          pure <- strsplit(narrow[i], "_")
          protList[[1]] <- append(protList[[1]], fixed(pure[[1]][1]), after = length(protList[[1]]))
          }
        }
      }
  }
  locList <- list(firstList <- unlist(protList),  ## prot
                      secondList <- unlist(nucList), ## nuc
                      thirdlist <- unlist(genList)) ## gen

  names(locList) <- c("prot","nuc","gen")

  pseudo <- IMGTHLAGeneTypes$Names[IMGTHLAGeneTypes$`Pseudogene/Fragment` == "Pseudogene"]
  pseudo <- pseudo[grep("HLA-",pseudo,fixed=TRUE)]
  pseudo <- sub("HLA-","",pseudo,fixed=TRUE)

  frag <- IMGTHLAGeneTypes$Names[IMGTHLAGeneTypes$`Pseudogene/Fragment` == "Fragment"]
  frag <- sub("HLA-","",frag,fixed=TRUE)

  locList$pseudo <- pseudo
  locList$frag <- frag
  locList$pseudo.noNuc <- locList$pseudo[is.element(locList$pseudo,locList$nuc) == FALSE]
  locList$pseudo.noGen <- locList$pseudo[is.element(locList$pseudo,locList$gen) == FALSE]
  locList$frag.noNuc <- locList$frag[is.element(locList$frag,locList$nuc) == FALSE]
  locList$frag.noGen <- locList$frag[is.element(locList$frag,locList$gen) == FALSE]

  locList
}

