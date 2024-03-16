##data v1.0.0 15MAR2024

##alleleListHistory
#'Identifies HLA loci that have amino acid, cDNA, or gDNA alignments.
#'
#'A large list object of two elements. The first element contains Allele List version information. The second element contains all releases of allele lists on the IMGT/HLA database.
#'This object is built by the UpdateAlleleListHistory() function in the HLAtools package.
#' @docType data
#' @name alleleListHistory
#' @usage data(alleleListHistory)
#' @format A large list of 2 elements that identify the table of allele names for each IPD-IMGT/HLA reference Database release and the source release
#' \itemize{
#'    \item(Version: the Date and IPD-IMGT HLA Release version under which the alleleListHistory object was made)
#'    \item(AlleleListHistory: a dataframe version of the alleleListHistory.txt file)
#' }
#' @source https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/Allelelist_history.txt
"alleleListHistory"

##fragmentFeatureNames
#'Identifies the gene features of HLA pseduogenes and gene fragments
#'
#'A list object of 15 elements. Each of the 15 elements corresponds to one of the 15 HLA pseudogenes and gene fragments, and contains two items. The first item identifies the gene features for that locus. The second item contains an annotation detailing information regarding non-standard gene structures.
#'This object is built by the ffN() function in the HLAtools package.
#' @docType data
#' @name fragmentFeatureNames
#' @usage data(fragmentFeatureNames)
#' @format A list of 15 elements that identify the gene features for each locus wih a genomic alignment and annotate non-standard features
#' @source https://hla.alleles.org/genes/index.html
#' @source https://github.com/ANHIG/IMGTHLA/tree/Latest/alignments
"fragmentFeatureNames"

##GLSC.ex
#'Example data frame containing Genotype List String Code data.
#'
#'A two-column data frame (modified from pould::hla.hap.demo) including GL String data in Gl String Code format.
#'Column one identifies each subject's status, while column two identifies each subject's HLA genotype in GL String Code format.
#'This data is provided for use in examples and demonstrations.
#' @docType data
#' @name GLSC.ex
#' @usage data(GLSC.ex)
#' @format A dataframe with 419 rows and 2 columns.
#' @source pould::hla.hap.demo
#' @references Mack et al. HLA 2023 Oct;102(4):501-507 https://doi.org/10.1111/tan.15145
"GLSC.ex"

##GLstring.ex
#'Example data frame containing Genotype List String data.
#'
#'A two-column example data frame (from pould::hla.hap.demo) including GL String Data.
#'Column one identifies each subject's status, while column two identifies each subject's HLA genotype in GL String format.
#'This data is provided for use in examples and demonstrations.
#' @docType data
#' @name GLstring.ex
#' @usage data(GLstring.ex)
#' @format A dataframe with 419 rows and 2 columns.
#' @source pould::hla.hap.demo
#' @references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
"GLstring.ex"

##HLAgazeteer
#'Identifies HLA loci that have amino acid, cDNA, or gDNA alignments.
#'
#'A list object of vectors that identify loci with amino acid (HLAgazeteer$prot), cDNA (HLAgazeteer$nuc), gDNA alignments (HLAgazeteer$gen), pseudogenes (HLAgazeteer$pseudo), gene fragments (HLAgazeteer$frag), as well as pseudogenes and gene fragments for which no cDNA ($noNuc) and gDNA ($noGen) alignments are available).
#'This object is built by the buildGazeteer() function in the HLAtools package.
#'
#' @note Pseudogenes and gene fragments derived from the IMGTHLAGeneTypes object.
#' @docType data
#' @name HLAgazeteer
#' @usage data(HLAgazeteer)
#' @format A list object containing 9 character vectors
#' @source https://github.com/ANHIG/IMGTHLA/tree/Latest/alignments
#' @source HLAtools::IMGTHLAGeneTypes
"HLAgazeteer"

##IMGTHLAGeneTypes
#'Describes the molecular characteristics of the genes curated by the IPD-IMGT/HLA Database, and identifies gene fragments and pseudogenes.
#'
#'A data frame of three columns identifying each gene supported by the IPD-IMGT/HLA Database, its molecular characteristics, and its status as gene fragment or pseudogene.
#'This object is built by the BuildIMGTHLAGeneTypes() function.
#'
#' @docType data
#' @name IMGTHLAGeneTypes
#' @usage data(IMGTHLAGeneTypes)
#' @format A data frame of three columns
#' @source https://github.com/ANHIG/IMGTHLA/tree/Latest/alignments
"IMGTHLAGeneTypes"

##UNIFORMAT.example
#'Example data frame containing UNIFORMAT data.
#'
#'An two-column example example data frame including UNIFORMAT data.
#'Column one contains sample identifiers, while column two identifies each subject's HLA genotype in UNIFORMAT format.
#'To be used for example runs and demonstrations.
#' @docType data
#' @name UNIFORMAT.example
#' @usage data(UNIFORMAT.example)
#' @format A dataframe with 20 rows and 2 columns.
#' @source https://hla-net.eu/wp/wp-content/uploads/example-three-loci.unif_.txt
#' @references Nunes Tissue Antigens 2007;69(s1):203-205 https://doi.org/10.1111/j.1399-0039.2006.00808.x
"UNIFORMAT.example"

##HLAatlas
#'identifies the boundary positions of exons, introns and UTRs in the amino acid, cDNA and gDNA alignments in the HLAalignments data object. 
#'
#'A list object of sub-lists of R dataframes (atlases) for each locus with a protein (prot), cDNA (nuc), and gDNA (gen) alignment. Each atlas identifies the position of the exon, intron or UTR boundaries in an alignment.
#'This object is built by the atlasFull() function in the HLAtools package.
#' @docType data
#' @name HLAatlas
#' @usage data(HLAatlas)
#' @format A list of 4 elements that include the gene features boundaries in each ANHIG/IMGTHLA sequence alingment. 
#' \itemize{
#'    \item(prot: peptide-alignment atlases)
#'    \item(nuc: cDNA-alignment atlases)
#'    \item(gen: gDNA-alignment atlases)
#'    \item(version: The IPD-IMGT/HLA Database release version under which these data were generated)
#'    }
#' @source https://github.com/ANHIG/IMGTHLA/tree/Latest/alignments
"HLAatlas"
