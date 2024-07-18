#buildGazeteer v6.0.0 10JUL2024

##############
##buildGazeteer
#'Define Categories of Genes Supported by the IPD-IMGT/HLA Database
#'
#'@description
#'Consumes information in the ANHIG/IMGTHLA GitHub repository and at hla.alleles.org/genes to define specific categories of genes supported by the IPD-IMGT/HLA Database, which are represented as nineteen elements of the HLAgazeteer object.
#'
#'Elements: 
#'*  All genes with alignments ($align)
#'*  Genes that that do and do not have amino acid ($prot/$noprot), nucleotide ($nuc/$nonuc), and genomic ($gen/$nogen) alignments 
#'*  HLA genes ($hla), pseudogenes ($pseudo), gene fragments ($frag)
#'*  Genes that are expressed ($expressed) or not expressed ($notexpressed)
#'*  Genes in the Class I region ($classireg), Class I HLA genes ($classIhla), Genes in the Class II region ($classiireg), and Class II HLA genes ($classiihla)
#'*  Classical HLA genes ($classical) and non-classical exprssed HLA genes ($nonclassical)
#'*  All genes presented in map order ($map)
#'
#'The twentieth element ($version) identifies the IPD-IMGT/HLA Database version used to build the HLAgazeteer.
#'
#'@param version A string identifying of the desired IPD-IMGT/HLA Database release version to which the gazeteer should be updated. The default value is most recent IPD-IMGT/HLA Database release version.
#'
#'@note The *$prot* and *$nuc* vectors include a 'DRB' "gene". While 'DRB' is not a gene name, the DRB_prot.txt file includes combined alignments for the DRB1, DRB3, DRB4, and DRB5 genes, and the DRB_nuc.txt file includes combined alignments for the DRB1, DRB2, DRB3, DRB4, DRB5, DRB6, DRB7, DRB8, and DRB9 genes. 'DRB' is included in these vectors for the purpose of validation when these combined alignments are desired.
#'
#'@return A list object of vectors organizing the genes in the IPD-IMGT/HLA Database into specific categories.
#'
#'@note For interal HLAtools use.
#'@note These elements are constructed using data compiled from hla.alleles.org/genes/, github.com/ANHIG/IMGTHLA/tree/Latest/alignments, DOI:10.2741/a317, DOI:10.1111/tan.15180, and Human Genome Assembly GRCh38.p14 reference assembly NC_000006.12.
#'@note Additional information about these genes can be found in HLAtools::IMGTHLAGeneTypes.
#'@note The $map element does not distinguish the order of the DRB3, DRB4 and DRB5 genes (DRB3/4/5) or the DR6 and DR7 (DR6/7) genes, as the individual elements of these gene sets are found on different DRB haplotypes. 
#'@note This function requires internet access to function.
#'
#'@importFrom stringr fixed
#'@importFrom utils read.table
#'
#'@source Andersson G. Evolution of the human HLA-DR region. Front. Biosci. 1998 Jul 27:3:d739-45.
#'@source Alexandrov et al. HLA-OLI: A new MHC class I pseudogene and HLA-Y are located on a 60 kb indel in the human MHC between HLA-W and HLA-J. HLA 2023 Nov; 102(5):599-606.
#'
#'@export
#'
buildGazeteer <- function(version = getLatestVersion()) {
  
  if(!validateVersion(version)) {stop(paste(version,"is not a valid IPD-IMGT/HLA Database release.",sep=" "))}
  
  URL <- paste("https://github.com/ANHIG/IMGTHLA/tree/",repoVersion(version),"/alignments",sep="")
  
  # Fix for version 3.13.0/3.13.1, where "3130" is used in the URL, but "3.13.1" is used in the AlleleListHistory file and alignments
  if(version == "3.13.1") { URL <- "https://github.com/ANHIG/IMGTHLA/tree/3130/alignments" }
  
  narrow <- getAlignmentNames(URL)
  
  nucList <- vector("list", 1)
  genList <- vector("list", 1)
  protList <- vector("list", 1)

 for(i in 1:length(narrow)) {
        if (grepl("Class", narrow[i], fixed = TRUE) == FALSE) {
           if (grepl("nuc", narrow[i], fixed = TRUE)) {
       #   if(!narrow[i] == "DRB_nuc.txt") { #Exclude combined "DRB" alignments
          pure <- strsplit(narrow[i], "_")
          nucList[[1]] <- append(nucList[[1]], fixed(pure[[1]][1]), after = length(nucList[[1]]))
      #    }
        }
        if (grepl("gen", narrow[i], fixed = TRUE)) {
          pure <- strsplit(narrow[i], "_")
          genList[[1]] <- append(genList[[1]], fixed(pure[[1]][1]), after = length(genList[[1]]))
         }
        if (grepl("prot", narrow[i], fixed = TRUE)) {
      #    if(!narrow[i] == "DRB_prot.txt") { #Exclude combined "DRB" alignments
          pure <- strsplit(narrow[i], "_")
          protList[[1]] <- append(protList[[1]], fixed(pure[[1]][1]), after = length(protList[[1]]))
     #     }
        }
      }
  }
  locList <- list(zerothList <- sort(unique(c(unlist(genList),unlist(nucList),unlist(protList),c("DRB2", "DRB6", "DRB7", "DRB9")))), #all aligned; these last four are in the DRB_nuc.txt file
                      firstList <- unlist(genList),  ## prot
                      secondList <- unlist(nucList), ## nuc
                      thirdlist <- unlist(protList)) ## gen

  names(locList) <- c("align","gen","nuc","prot")
  
  locList
  
  locList$align <- locList$align[!locList$align %in% "DRB"] # remove "DRB" for the list of genes with alignments

  pseudo <- IMGTHLAGeneTypes$GeneTypes$Names[IMGTHLAGeneTypes$GeneTypes$`Pseudogene/Fragment` == "Pseudogene"]
  pseudo <- sub("HLA-","",pseudo,fixed=TRUE)

  frag <- IMGTHLAGeneTypes$GeneTypes$Names[IMGTHLAGeneTypes$GeneTypes$`Pseudogene/Fragment` == "Fragment"]
  frag <- sub("HLA-","",frag,fixed=TRUE)
  
  hla <-IMGTHLAGeneTypes$GeneTypes$Names[grep("HLA-",IMGTHLAGeneTypes$GeneTypes$Names,fixed=TRUE)]
  hla <- sort(sub("HLA-","",hla,fixed=TRUE))

  locList$nogen <- sort(c(pseudo[is.element(pseudo,locList$gen) == FALSE],frag[is.element(frag,locList$gen) == FALSE]))
  locList$nonuc <- sort(c(pseudo[is.element(pseudo,locList$nuc) == FALSE],frag[is.element(frag,locList$nuc) == FALSE]))
  locList$noprot <- unique(sort(c(locList$gen[!locList$gen %in% locList$prot],locList$nonuc)))

  locList$pseudo <- pseudo
  locList$frag <- frag
  
  locList$hla <- hla
  
  expressed <- IMGTHLAGeneTypes$GeneTypes$Names[IMGTHLAGeneTypes$GeneTypes$`Pseudogene/Fragment`==""]
  expressed <- expressed[!expressed %in% "HLA-DQB3"]
  locList$expressed <- sub("HLA-","",expressed,fixed=TRUE)
  locList$notexpressed <- locList$hla[!locList$hla %in% locList$expressed]
  
  classireg <- sort(c(IMGTHLAGeneTypes$GeneTypes$Names[grep("Class I",IMGTHLAGeneTypes$GeneTypes$`Molecular Characteristics`,fixed=TRUE)],IMGTHLAGeneTypes$GeneTypes$Names[grep("class I",IMGTHLAGeneTypes$GeneTypes$`Molecular Characteristics`,fixed=TRUE)]))
  classiireg <- IMGTHLAGeneTypes$GeneTypes$Names[!IMGTHLAGeneTypes$GeneTypes$Names %in% classireg]
  classihla <- sub("HLA-","",classireg[grep("HLA-",classireg,fixed=TRUE)],fixed = TRUE)
  classiihla <- sub("HLA-","",classiireg[grep("HLA-",classiireg,fixed=TRUE)],fixed = TRUE)
  classireg <- sort(sub("HLA-","",classireg[classireg != c("HLA-Z")],fixed=TRUE)) # HLA-Z is a class I gene fragment in the class II region
  classiireg <- sort(c(sub("HLA-","",classiireg,fixed=TRUE),"Z"))
  
  locList$classireg <- classireg
  locList$classihla <- classihla
  locList$classiireg <- classiireg
  locList$classiihla <- classiihla
  
  locList$classical <- sort(c("A","B","C","DRA","DRB1","DRB3","DRB4","DRB5","DQA1","DQB1","DPA1","DPB1"))
  locList$nonclassical <- sort(c("F","G","E","DQA2","DQB2","DOB","DMB","DMA","DOA","DPA2","DPB2"))
  
  locList$map <- c("HFE","F","MICE","V","P","G","H","T","K","U","A","W","MICD","Y","R","J","L","N","MICC","E","C","B","S","MICA","X","MICB","DRA","DRB9","DRB3/4/5","DRB8","DRB6/7","DRB2","DRB1","DQA1","DQB1","DQB3","DQA2","DQB2","DOB","TAP2","PSMB8","TAP1","PSMB9","Z","DMB","DMA","DOA","DPA1","DPB1","DPA2","DPB2","DPA3")

  locList$version <- version
    
  locList
}

##############
##getAlignmentNames
#'Retrieve Alignment Filenames for HLA Genes
#'
#'Retrieves the filenames of the protein, nucleotide and genomic alignments available a specific branch of the IMGTHLA GitHub Repository
#'
#'@param URL A chararacter string containing a Uniform Resource Locator (URL) identifying of the desired IPD-IMGT/HLA Database release version from which the alignment filenames should be retrieved.
#'
#'@return A character vector of all of the filenames. 
#'
#'@importFrom stringr str_detect
#'
#'@note For internal HLAtools use.
#'
#'@export
#'
getAlignmentNames <- function(URL){
  on.exit(closeAllConnections())
  rough <- readLines(URL)
  rough <- strsplit(rough, split = "[,]|[:]|[\\{]|[\\}]")
  
  roughList <- c()
  k <- 0
  
  for(i in 1:length(rough)) {
    if(length(rough[[i]])>0) {
      for(j in 1:length(rough[[i]])) {
        if(nchar(rough[[i]][j])!= 0){
          if(any(!is.na(str_detect(rough[[i]][j],c("alignments"))))) {
            if(str_detect(rough[[i]][j],c("alignments")) && str_detect(rough[[i]][j],c(".txt"))) {
              roughList <- append(roughList,rough[[i]][j],length(roughList))
              k <- k+1
              roughList[k] <- substr(roughList[k],13,nchar(roughList[k])-1)
            }
          }        
        }
      }
    }
  }
  
  roughList <- roughList[!is.na(roughList)]
  roughList <- unique(roughList)
  roughList[!str_detect(roughList[],"Link--primary")]
}
