##data v2.0.0 10JUL2024

##alleleListHistory
#'Allele Names Across All Release Versions
#'
#'A large list object of two elements. The first element contains version information. The second element contains all releases of allele name lists in the IMGT/HLA database.
#'This object is built by the UpdateAlleleListHistory() function.
#'
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
#'Gene Features of HLA Pseduogenes and Gene Fragments
#'
#'A list object of 15 elements. Each of the 15 elements corresponds to one of the 15 HLA pseudogenes and gene fragments, and contains two items. The first item identifies the gene features for that locus. The second item contains an annotation detailing information regarding non-standard gene structures.
#'This object is built by the ffN() function.
#'
#' @docType data
#' @name fragmentFeatureNames
#' @usage data(fragmentFeatureNames)
#' @format A list of 15 elements that identify the gene features for each locus wih a genomic alignment and annotate non-standard features
#' @source https://hla.alleles.org/genes/index.html
#' @source https://github.com/ANHIG/IMGTHLA/tree/Latest/alignments
"fragmentFeatureNames"

##GLSC.ex
#'Example Data Frame of Genotype List String Code Data
#'
#'A two-column data frame (modified from pould::hla.hap.demo) including GL String data in Gl String Code format.
#'Column one identifies each subject's status, while column two identifies each subject's HLA genotype in GL String Code format.
#'This data is provided for use in examples and demonstrations.
#'
#' @docType data
#' @name GLSC.ex
#' @usage data(GLSC.ex)
#' @format A dataframe with 419 rows and 2 columns.
#' @source pould::hla.hap.demo
#' @references Mack et al. HLA 2023 Oct;102(4):501-507 https://doi.org/10.1111/tan.15145
"GLSC.ex"

##GLstring.ex
#'Example Data Frame of Genotype List String Data.
#'
#'A two-column example data frame (from pould::hla.hap.demo) including GL String Data.
#'Column one identifies each subject's status, while column two identifies each subject's HLA genotype in GL String format.
#'This data is provided for use in examples and demonstrations.
#'
#' @docType data
#' @name GLstring.ex
#' @usage data(GLstring.ex)
#' @format A dataframe with 419 rows and 2 columns.
#' @source pould::hla.hap.demo
#' @references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
"GLstring.ex"

##HLAgazeteer
#'Functional and Organizational Categories of Genes Supported by the IPD-IMGT/HLA Database
#'
#'A list object of nineteen vectors that identify genes in the HLA region that share specific commonalities. 
#'This object is built by the buildGazeteer() function.
#'
#' @docType data
#' @name HLAgazeteer
#' @usage data(HLAgazeteer)
#' @format A large list of 19 vectors that define specific categories of genes supported by the IPD-IMGT/HLA Database
#' \itemize{
#'.   \item(align: all genes with alignments in the IPD/IMGT-HLA GitHub Repository)
#'    \item(gen: genes with genomic alignments in the IPD/IMGT-HLA GitHub Repository)
#'    \item(nuc: genes with nucleotide alignments in the IPD/IMGT-HLA GitHub Repository)
#'    \item(prot: genes with protein alignments in the IPD/IMGT-HLA GitHub Repository)
#'    \item(nogen: genes with no genomic alignments)
#'    \item(nonuc: genes with no nucleotide alignments)
#'    \item(noprot: genes with no protein alignments)
#'    \item(pseudo: pseudogenes)
#'    \item(frag: gene fragments)
#'    \item(hla: HLA genes)
#'    \item(expressed: genes that are expressed)
#'    \item(notexpressed: genes that are not expressed)
#'    \item(classireg: genes found in the HLA class I region)
#'    \item(classihla: class I HLA genes)
#'    \item(classiireg: genes found in the HLA class II region)
#'    \item(classiihla: class II HLA genes)
#'    \item(classical: classical HLA genes)
#'    \item(nonclassical: non-classical HLA genes)
#'    \item(map: all genes organized by 5' to 3' map order on the genomic reference + strand)
#'    \item(version: IPD-IMGT/HLA Database version used to build the HLAgazeteer)
#' }
#' @source https://hla.alleles.org/genes
#' @source https://github.com/ANHIG/IMGTHLA/tree/Latest/alignments
#' @source Andersson Front. Biosci. 1998, 3(4), 739–745. https://doi.org/10.2741/a317
#' @source Alexandrov et al. HLA 2023, Vol.102(5), p.599-606. https://doi.org/10.1111/tan.15180
#' @source https://www.ncbi.nlm.nih.gov/nucleotide/NC_000009.12
"HLAgazeteer"

##IMGTHLAGeneTypes
#'Molecular characteristics of the Genes Curated by the IPD-IMGT/HLA Database
#'
#'A data frame of three columns identifying each gene supported by the IPD-IMGT/HLA Database, its molecular characteristics, and its status as a gene fragment or pseudogene.
#'This object is built by the BuildIMGTHLAGeneTypes() function.
#'
#' @docType data
#' @name IMGTHLAGeneTypes
#' @usage data(IMGTHLAGeneTypes)
#' @format A list object of two elements:
#' \itemize{
#'    \item(GeneTypes: a data frame of three columns)
#'    \item(version: a character string identifying the date that the source file was written)
#' }
#' @source https://hla.alleles.org/genes/index.html
"IMGTHLAGeneTypes"

##UNIFORMAT.example
#'Example Data Frame of UNIFORMAT Data.
#'
#'A two-column example data frame including UNIFORMAT data.
#'Column one contains sample identifiers, while column two identifies each subject's HLA genotype in UNIFORMAT format.
#'To be used for example runs and demonstrations.
#'
#' @docType data
#' @name UNIFORMAT.example
#' @usage data(UNIFORMAT.example)
#' @format A dataframe with 20 rows and 2 columns.
#' @source https://hla-net.eu/wp/wp-content/uploads/example-three-loci.unif_.txt
#' @references Nunes Tissue Antigens 2007;69(s1):203-205 https://doi.org/10.1111/j.1399-0039.2006.00808.x
"UNIFORMAT.example"

##HLAatlas
#'Boundary Positions of Exons, Introns and UTRs in Amino Acid, cDNA and gDNA Alignments 
#'
#'A list object of sub-lists of R dataframes (atlases) for each locus with a protein (prot), cDNA (nuc), and gDNA (gen) alignment. Each atlas identifies the position of the exon (E), intron (I) or UTR (U) gene-feature boundaries in an alignment. Boundaries for non-standard hybrid (H), join (J), novel (N) and segment (S) features may be included in gene fragment and pseudogene atlases.
#'This object is built by the atlasFull() function.
#'
#' @docType data
#' @name HLAatlas
#' @usage data(HLAatlas)
#' @format A list of 4 elements that include the gene features boundaries in each ANHIG/IMGTHLA sequence alignment. 
#' \itemize{
#'    \item(prot: peptide-alignment atlases)
#'    \item(nuc: cDNA-alignment atlases)
#'    \item(gen: gDNA-alignment atlases)
#'    \item(version: The IPD-IMGT/HLA Database release version under which these data were generated)
#'    }
#' @source https://github.com/ANHIG/IMGTHLA/tree/Latest/alignments
"HLAatlas"

##sHLAdata
#'Synthetic HLA Data for use with Package Examples
#'
#'A BIGDAWG-formatted data frame of 18 columns and 47 rows containing synthetic HLA-A, -C, -B, -DRB1, -DQA1, -DQB1, -DPA1, and -DPB1 genotype data for 24 control subjects and 23 case subjects. Allele name data were recorded under IPD-IMGT/HLA Database release version 3.56.0. These are synthetic data generated for the imaginary UchiTelle population, and do not represent biologically-derived HLA genotypes. 
#'
#' @docType data
#' @name sHLAdata
#' @usage data(sHLAdata)
#' @format A data frame of 18 columns and 47 rows.
#' @source These synthetic data were generated as part of the 13th International HLA Workshop to demonstrate PyPop functions. See <http://pypop.org/docs/guide-chapter-usage.html#data-minimal-noheader> for additional details.
"sHLAdata"