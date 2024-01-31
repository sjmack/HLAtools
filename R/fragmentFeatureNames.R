##Fragment Feature Names v3.5 January 25, 2024

################
##FragmentFeatureNames (ffN)
#'Identifies and annotates the gene features present in genomic alignments of HLA pseudogenes and gene fragments.
#'
#'HLA pseudogenes and gene fragments many not share all of the gene features of their expressed homologs. The ffN() function generates a list object identifying the features present in a given gene fragment or pseudogene with an alignment file. This includes annotations of the differences between the pseudogenes and gene fragments, and their expressed homologs.
#'Pseuogenes and gene fragments are identified in the IMGTHLAGeneTypes data object.
#'
#'@param version A character string identifying the pertinent IPD-IMGT/HLA Database release version (e.g., "3.55.0") under which the returned object is generated. This parameter does not impact the generation of the feature names or annotations, and is only included to provide IPD-IMGT/HLA Database release version context.
#'
#'@return A list object where each element is the name of a pseudogene or gene fragment. Each of these elements is a list of 'features' and 'annotation'. "Features" identifies the gene features for that gene, ordered from 5' to 3'. Because these genes are not-expressed, the three standard gene features (U, Untranslated Region (UTR); E, Exon; and I, Intron) may not always be present. In these cases, H, J, N and S features (defined below) are returned. Each feature identifier is followed by an identifying number (e.g. U.5 is the 5' UTR, and E.3 is Exon 3). "Annotation" provides the composition of the non-standard H, J, N and S features for each gene.
#'
#'@note Standard Features
#'@note E - Exon, a peptide-encoding sequence
#'@note I - Intron, an intervening sequence found between Exons
#'@note U - UTR, an untranslated region of sequence preceding the first Exon or following the last Exon
#'
#'@note Additional Features used in these annotations, based on boundaries indicated in the sequence alignment, are:
#'@note H - Hybrid, a sequence that includes at least one known feature sequence and one nucleotide sequence that does not correspond to a known feature
#'@note J - Join, a sequence that includes two or more features that are separated by a feature boundary in the reference
#'@note N - Novel, a novel sequence that does correspond to a known feature sequence
#'@note S - Segment, a subset of a longer feature sequence
#
#'@note Reference Genes
#'@note HLA-C genomic sequence is used as the reference for class I pseudogenes and gene fragments.
#'@note HLA-DPA1 and -DPB1 genomic sequences are used as the references for class II pseudogenes -DPA2 and -DPB2, respectively.
#'
#'@export
#'
#'@examples fragmentFeatureNames <- ffN("3.35.0")
#'
#'@note Features and their annotations have been identified manually, and this function is for internal HLAtools use.

ffN <- function(version){ #Fragment Feature Names

  if(paste("X",gsub("\\.","",as.character(version)),sep="") %in% colnames(alleleListHistory$AlleleListHistory)) {

  featureNames <- rep(list(list()), length(c(HLAgazeteer$pseudo[!is.element(HLAgazeteer$pseudo,HLAgazeteer$pseudo.noGen)],HLAgazeteer$frag[!is.element(HLAgazeteer$frag,HLAgazeteer$frag.noGen)])))
  names(featureNames) <-sort(c(HLAgazeteer$pseudo[!is.element(HLAgazeteer$pseudo,HLAgazeteer$pseudo.noGen)],HLAgazeteer$frag[!is.element(HLAgazeteer$frag,HLAgazeteer$frag.noGen)]))


  featureNames$DPA2$features <- c("U.5","E.1","I.1","E.2","I.2","E.3","I.3","E.4","U.3") #  8 boundaries -- **Unconfirmed
  featureNames$DPB2$features <- c("U.5","E.1","I.1","E.2","I.2","E.3","I.3","E.4","I.4","E.5","U.3") # 10 boundaries
  featureNames$H$features <- c("U.5","E.1","I.1","E.2","I.2","E.3","I.3","E.4","I.4","E.5","I.5","E.6","I.6","E.7","I.7","E.8","U.3") # 16 boundaries full length
  featureNames$J$features <- c("U.5","E.1","I.1","E.2","I.2","E.3","I.3","E.4","I.4","E.5","I.5","E.6","I.6","E.7","I.7","E.8","U.3") # 16 boundaries full length
  featureNames$K$features <- c("U.5","E.1","I.1","E.2","I.2","E.3","I.3","E.4","I.4","E.5","I.5","E.6","I.6","E.7","I.7","E.8","U.3") # 16 boundaries full length
  featureNames$L$features <- c("U.5","E.1","I.1","E.2","I.2","E.3","I.3","E.4","I.4","E.5","I.5","E.6","I.6","E.7","I.7","E.8","U.3") # 16 boundaries full length
  featureNames$N$features <- c("N.1","H.1","N.2") #  2 boundaries
  featureNames$P$features <- c("J.1","E.3","I.3","E.4","I.4","E.5","I.5","E.6","I.6","E.7","I.7","U.3") #11 boundaries
  featureNames$S$features <- c("H.1","S.1","J.1","E.7","I.7","E.9","S.2") #  6 boundaries
  featureNames$T$features <- c("H.1","E.4","I.4","E.5","I.5","E.6","I.6","E.7","U.3") #  8 boundaries
  featureNames$U$features <- c("N.1","J.1","S.1") #  2 boundaries
  featureNames$V$features <- c("U.5","E.1","I.1","E.2","I.2","S.1","N.1") #6 boundaries
  featureNames$W$features <- c("J.1","E.3","I.3","E.4","I4","E.5","I.5","E.6","I.6","E.7","I.7","N.1","S.1") #12 boundaries
  featureNames$Y$features <- c("U.5","E.1","I.1","E.2","I.2","E.3","I.3","E.4","I.4","E.5","I.5","E.6","I.6","E.7","I.7","E.8","U.3") #16 boundaries full length

  featureNames$DPA2$annotation <- "All of the reference gene features are present. E.1 starts 28 nucleotides before the reference, and ends 88 nucleotides before the reference."
  featureNames$DPB2$annotation <- "All of the reference gene features are present. I.2 contains a very large insertion relative to the reference."
  featureNames$H$annotation <- "All of the reference gene features are present."
  featureNames$J$annotation <- "All of the reference gene features are present."
  featureNames$K$annotation <- "All of the reference gene features are present."
  featureNames$L$annotation <- "All of the reference gene features are present."
  featureNames$N$annotation <- "N.1 is 51 nucleotides of novel sequence that does not align to Exon 3. (the last 40 nucleotides align poorly). The first 34 nucleotides of H.1 do not align to the Exon 4 reference; the last 135 nucleotides align to Exon 4. N.2 does not align to Exon 4 or Intron 4."
  featureNames$P$annotation <- "J.1 is ~350 nucleotides of novel sequence followed by ~120 nucleotides from the 5' end of Intron 2. Exon 8 is absent."
  featureNames$S$annotation <- "H.1 is 37 nucleotides of novel sequence, followed by the last 185 nucleotides of Intron 5. S.1 is the first 27 nucleotides of Exon 6. J.1 is the last 4 nucleotides of Exon 6 (2 nucleotides in the reference have been deleted), followed by the last 100 nucleotides of Intron 6 (6 nucleotides in the refernece have been deleted). E.9 in a 191 nucleotide long Exon in what is the 5' end of the 3' UTR in the reference. S.2 is the 3' end of the 3' UTR."
  featureNames$T$annotation <- "H.1 is 31 nucleotides of novel sequence, followed the last 93 nucleotides of Exon 3, and all of Intron 3."
  featureNames$U$annotation <- "N.1 is 54 nucleotides of novel sequence. J.1 includes 196 nucleotdes that align to Exon 3, and the first 10 nucleotides of Intron 3. S.1 is 301 nucleotides of Intron 3. HLA-U*01:01:01:02, *01:02, 01:03 and 01:04 include a 41 nucleotide novel sequence insertion. U*01:04 includes an additional 145 nucleotides of novel sequence at the 3' end of S.1."
  featureNames$V$annotation <- "S.1 is the first 195 nucleotides of Exon 3. N.1 is 524 nucleotides of novel sequence that does not align to Exon or Intron 3."
  featureNames$W$annotation <- "J.1 is ~365 nucleotides of novel sequence followed by ~120 nucleotides from the 5' end of Intron 2. N1 is a novel exon that follows the deletion of Exon 8. S1 is the remainder of the 3' UTR"
  featureNames$Y$annotation <- "U.5 is a single nucleotide, and the U.3 is 171 nucleotides long in HLA-Y*02:01. There is no 5' or 3' UTR in HLA-Y*01:01 or *03:01."

  featureNames$version <- version

  featureNames  } else {message(paste(version, "is not a recognized release version.",sep=" "))}

}

