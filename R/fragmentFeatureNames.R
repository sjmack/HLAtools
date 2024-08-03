##Fragment Feature Names v4.0 July 12, 2024

################
##FragmentFeatureNames (ffN)
#'Identify and Annotate Gene Features in Pseudogenes and Gene Fragments.
#'
#'@description
#'HLA pseudogenes and gene fragments many not share all of the gene features of their expressed homologs. The ffN() function generates a list object identifying the features present in a given gene fragment or pseudogene with a genomic alignment file. This includes annotations of the differences between the pseudogenes and gene fragments, and their expressed homologs.
#'Pseudogenes and gene fragments are identified in the IMGTHLAGeneTypes data object.
#'
#'Standard Features
#'* E - Exon, a peptide-encoding sequence
#'* I - Intron, an intervening sequence found between Exons
#'* U - UTR, an untranslated region of sequence preceding the first Exon or following the last Exon
#'
#'Additional Features used in these annotations, based on boundaries indicated in the sequence alignment, are:
#'* H - Hybrid, a sequence that includes at least one known feature sequence and one nucleotide sequence that does not correspond to a known feature
#'* J - Join, a sequence that includes two or more features that are separated by a feature boundary in the reference
#'* N - Novel, a novel sequence that does correspond to a known feature sequence
#'* S - Segment, a subset of a longer feature sequence
#'
#'Reference Genes
#'* HLA-C genomic sequence is used as the reference for class I pseudogenes and gene fragments.
#'* HLA-DPA1 and -DPB1 genomic sequences are used as the references for class II pseudogenes -DPA2 and -DPB2, respectively.
#'
#'@param version A character string identifying the pertinent IPD-IMGT/HLA Database release version (e.g., "3.55.0") under which the returned object is generated. This parameter does not impact the generation of the feature names or annotations, and is only included to provide IPD-IMGT/HLA Database release version context.
#'
#'@return A list object where each element is the name of a pseudogene or gene fragment. Each of these elements is a list of 'features' and 'annotation'. "Features" identifies the gene features for that gene, ordered from 5' to 3'. Because these genes are not-expressed, the three standard gene features (U, Untranslated Region (UTR); E, Exon; and I, Intron) may not always be present. In these cases, H, J, N and S features (defined below) are returned. Each feature identifier is followed by an identifying number (e.g. U.5 is the 5' UTR, and E.3 is Exon 3). "Annotation" provides the composition of the non-standard H, J, N and S features for each gene.
#'#'
#'@export
#'
#'@examples 
#'fragmentFeatureNames <- ffN("3.35.0")
#'
#'@note Features and their annotations have been identified manually. Feature annotations will not change unless a new pseudogene or gene fragment is added in a future release, in which case new annotations will be generated.
#'@note The H, J, N and S features are described for class I pseudogenes and gene fragment sequences, all of which are described relative to the HLA-C reference sequence. Feature length differences for DPA2 and DPB2, relative to the DPA1 and DPB1 references, are noted in annotations of standard feature abbreviations (E, I and U).
#'@note No annotations are included for the DRB2, DRB6, DRB7, and DRB9 genes, as genomic alignments for these genes are not available.
#'@note For internal HLAtools use.
#'
ffN <- function(version){ #Fragment Feature Names
  
  if(version != "Latest"){ #
    if(!validateVersion(version)){stop(paste(version," is not a valid IPD-IMGT/HLA Database release version."))}
  }else{ version <- getLatestVersion()}

  if(paste("X",gsub("\\.","",as.character(version)),sep="") %in% colnames(alleleListHistory$AlleleListHistory)) {

  featureNames <- rep(list(list()), length(HLAgazeteer$align[HLAgazeteer$align %in% unique(c(HLAgazeteer$pseudo,HLAgazeteer$frag))][!HLAgazeteer$align[HLAgazeteer$align %in% unique(c(HLAgazeteer$pseudo,HLAgazeteer$frag))] %in% HLAgazeteer$nogen]))
  names(featureNames) <-sort(HLAgazeteer$align[HLAgazeteer$align %in% unique(c(HLAgazeteer$pseudo,HLAgazeteer$frag))][!HLAgazeteer$align[HLAgazeteer$align %in% unique(c(HLAgazeteer$pseudo,HLAgazeteer$frag))] %in% HLAgazeteer$nogen])


  featureNames$DPA2$features <- c("U.5","E.1","I.1","E.2","I.2","E.3","I.3","E.4","U.3") #  8 boundaries -- **Unconfirmed
  featureNames$DPB2$features <- c("U.5","E.1","I.1","E.2","I.2","E.3","I.3","E.4","I.4","E.5","U.3") # 10 boundaries
  featureNames$H$features <- c("U.5","E.1","I.1","E.2","I.2","E.3","I.3","E.4","I.4","E.5","I.5","E.6","I.6","E.7","I.7","E.8","U.3") # 16 boundaries full length
  featureNames$J$features <- c("U.5","E.1","I.1","E.2","I.2","E.3","I.3","E.4","I.4","E.5","I.5","E.6","I.6","E.7","I.7","E.8","U.3") # 16 boundaries full length
  featureNames$K$features <- c("U.5","E.1","I.1","E.2","I.2","E.3","I.3","E.4","I.4","E.5","I.5","E.6","I.6","E.7","I.7","E.8","U.3") # 16 boundaries full length
  featureNames$L$features <- c("U.5","E.1","I.1","E.2","I.2","E.3","I.3","E.4","I.4","E.5","I.5","E.6","I.6","E.7","I.7","E.8","U.3") # 16 boundaries full length
  featureNames$N$features <- c("N.1","H.1","N.2") #  2 boundaries
  featureNames$P$features <- c("H.1","E.3","I.3","E.4","I.4","E.5","I.5","E.6","I.6","E.7","I.7","U.3") # 11 boundaries
  featureNames$R$features <- c("H.1","E.3","I.3","E.4","I.4","E.5","I.5","E.6","I.6","J.1","J.2") # 10 boundaries
  featureNames$S$features <- c("H.1","S.1","J.1","E.7","I.7","E.9","S.2") #  6 boundaries
  featureNames$T$features <- c("H.1","E.4","I.4","E.5","I.5","E.6","I.6","E.7","U.3") #  8 boundaries
  featureNames$U$features <- c("N.1","J.1","S.1") #  2 boundaries
  featureNames$V$features <- c("U.5","E.1","I.1","E.2","I.2","S.1","N.1") #6 boundaries
  featureNames$W$features <- c("H.1","E.3","I.3","E.4","I4","E.5","I.5","E.6","I.6","E.7","I.7","N.1","S.1") #12 boundaries
  featureNames$Y$features <- c("U.5","E.1","I.1","E.2","I.2","E.3","I.3","E.4","I.4","E.5","I.5","E.6","I.6","E.7","I.7","E.8","U.3") #16 boundaries full length

  featureNames$DPA2$annotation <- "All of the reference gene features are present. E.1 starts 28 nucleotides before the reference, and ends 88 nucleotides before the reference."
  featureNames$DPB2$annotation <- "All of the reference gene features are present. I.2 contains a very large insertion relative to the reference."
  featureNames$H$annotation <- "All of the reference gene features are present."
  featureNames$J$annotation <- "All of the reference gene features are present."
  featureNames$K$annotation <- "All of the reference gene features are present."
  featureNames$L$annotation <- "All of the reference gene features are present."
  featureNames$N$annotation <- "N.1 is 51 nucleotides of novel sequence that does not align to Exon 3. The last 40 nucleotides align poorly. The first 34 nucleotides of H.1 do not align to the Exon 4 reference; the last 135 nucleotides align to Exon 4. N.2 does not align to Exon 4 or Intron 4."
  featureNames$P$annotation <- "H.1 is 287 nucleotides of novel sequence followed by 168 nucleotides of the 5' end of Intron 2, which include a 15 nucleotide section of novel sequence. Exon 8 is absent."
  featureNames$R$annotation <- "H.1 is 287 nucleotides of novel sequence followed by 156 nucleotides of the 5'end of Intron 2, which includes a 5 nucleotide section of novel sequence. J.1 is all of Exon 7 and the first 12 nucleotides of Intron 7. J.2 is the remainder of Intron 7, and the 3' UTR. Exon 8 is absent."
  featureNames$S$annotation <- "H.1 is 37 nucleotides of novel sequence, followed by the last 185 nucleotides of Intron 5. S.1 is the first 27 nucleotides of Exon 6. J.1 is the last 4 nucleotides of Exon 6 (2 nucleotides in the reference have been deleted), followed by the last 100 nucleotides of Intron 6 (6 nucleotides in the reference have been deleted). E.9 in a 191 nucleotide long Exon in what is the 5' end of the 3' UTR in the reference. S.2 is the 3' end of the 3' UTR."
  featureNames$T$annotation <- "H.1 is 31 nucleotides of novel sequence, followed the last 93 nucleotides of Exon 3, and all of Intron 3."
  featureNames$U$annotation <- "N.1 is 54 nucleotides of novel sequence. J.1 includes 196 nucleoitdes that align to Exon 3, and the first 10 nucleotides of Intron 3. S.1 is 301 nucleotides of Intron 3. HLA-U*01:01:01:02, *01:02, 01:03 and 01:04 include a 41 nucleotide novel sequence insertion. U*01:04 includes an additional 145 nucleotides of novel sequence at the 3' end of S.1."
  featureNames$V$annotation <- "S.1 is the first 195 nucleotides of Exon 3. N.1 is 524 nucleotides of novel sequence that does not align to Exon or Intron 3."
  featureNames$W$annotation <- "H.1 is 327 nucleotides of novel sequence followed by 148 nucleotides of the 5' end of Intron 2, which include a 15 nucleotide section of novel sequence. N1 is a novel exon that follows the deletion of Exon 8. S1 is the remainder of the 3' UTR"
  featureNames$Y$annotation <- "U.5 is a single nucleotide, and U.3 is 171 nucleotides long in HLA-Y*02:01. There is no 5' or 3' UTR in HLA-Y*01:01 or *03:01."

  featureNames$version <- version

  featureNames  } else {message(paste(version, "is not a recognized release version.",sep=" "))}

}

