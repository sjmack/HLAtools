##AminoAcidSequenceSearch v0.5 21AUG23

##########
##AAalign
#'Generate an amino acid alignment.
#'
#'Generates an amino acid alignment table for user-specified HLA alleles at user-specified amino acid positions.
#'
#'@param alleles A vector of un-prefixed HLA locus names.
#'@param positions either a vector of amino acid positions, against which all loci will be aligned, or a list of vectors of AA positions, exactly one vector for each allele, against which each corresponding allele will be aligned.
#'
#'@return A data frame of allele names and the corresponding peptide sequence for each amino acid position. an error message is returned if input loci is not available in the ANHIG/IMGTHLA Github Repository.
#'
#'@export
#'
#'@examples
#'AAalign(c("DRB1*01:01","DQB1*02:01","DPB1*01:01"),c(1,2,3,7,8,9,13,14,15))
#'AAalign(c("DQA1*01:01","DQB1*05:01","DPB1*01:01","DPA1*01:03"),list(32:58,33:59,31:57,29:55))
#'
AAalign <- function(alleles,positions){
  #print(alleles)
  #makes sure no duplicate alleles will be present in table
  alleles <- unique(alleles)
  #print(alleles)
  if(is.list(positions)) {
    #calls function below
    multipleAAalign(alleles,positions)
  }
  else {
    #calls function below
    singleAAalign(alleles,positions)
  }
}

############
##singleAAalign
#'Generate an amino acid alignment for specific alleles at specific positions.
#'
#'Generates an amino acid alignment at a single set of amino acid positions for HLA alleles at one or more loci.
#'
#'@param alleles A vector of un-prefixed HLA locus names.
#'@param positions A vector of amino acid positions, against which all loci will be aligned.
#'
#'@return A data frame of allele names and the corresponding peptide sequence for each amino acid position. a peptide sequence is not returned for a specific allele if input allele is not available in the ANHIG/IMGTHLA Github Repository. position will be left empty if specific allele does not have a peptide at the input position.
#'
#'@note For internal HLAtools use only.
#'
#'@export
#'
#'@examples
#'singleAAalign(c("DRB1*01:01","DQB1*02:01","DPB1*01:01"),c(32:58))
#'
singleAAalign <- function(alleles,positions){

  #making an array with correct dimensions
  #positions on the top, alleles on the side
  #as many rows as alleles, as many columns as positions
  align <- data.frame(matrix(" ", nrow = length(alleles), ncol = length(positions)+1),stringsAsFactors = FALSE)
  #naming first column
  #spaces under this column will be filled with actual Allele names
  colnames(align)[1] <- "Allele"
  #naming remaining columns, running through positions array and naming each column
  #spaces under this column will be filled with amino acids???????????????????????????????????????????
  colnames(align)[2:(length(positions)+1)] <- positions

  #running through array of alleles
  for(i in 1:length(alleles)){
    #setting the names for first column
    align[i,1] <- alleles[i]
    seqStr <- suppressMessages(AAsearch(alleles[i],positions,prefix=FALSE))
    #filling in the table (ex. Y N Q)
    align[i,2:(length(positions)+1)] <- if(!is.null(seqStr)) {strsplit(seqStr,split="~",fixed = TRUE)[[1]]} else {rep(" ",length(positions))}

  }
  #returning this array
  align
}

#############
##multipleAAalign
#'Generate an amino acid alignment for specific HLA alleles at different sets of positions.
#'
#'Generates an amino acid alignment for HLA alleles allowing each allele to be aligned to a different set of positions.
#'
#'@param alleles A vector of un-prefixed HLA locus names.
#'@param positions A list of vectors of AA positions, exactly one vector for each allele, against which each corresponding allele will be aligned.
#'
#'@return A data frame of allele names and the corresponding peptide sequence for each amino acid position for the specified positions designated for an allele. a peptide sequence is not returned for a specific allele if input allele is not available in the ANHIG/IMGTHLA Github Repository. position will be left empty if specific allele does not have a peptide at the input position.
#'
#'@note For internal HLAtools use only.
#'
#'@export
#'
#'@examples
#'multipleAAalign(c("DQA1*01:01","DQB1*05:01","DPB1*01:01"),list(32:58,33:59,31:57))
#'
multipleAAalign <- function(alleles,positions){
  #initial parameter check
  if(length(alleles) != length(positions)) {return(warning("The number of sets of positions must exactly match the number of alleles, please adjust input accordingly."))}
  #creating an array
  #as many rows as 2*alleles-1, as many columns as positions
    #rows are 2x because there is new frame for each, so there are multiple header columns
    #ex. DQA1 has its own frame (32,58)
    #ex. DQB1 has its own frame (33, 59)

  align <- data.frame(matrix(" ", nrow = (length(alleles)*2)-1, ncol = max(lengths(positions))+1),stringsAsFactors = FALSE)
  colnames(align)[1] <- strsplit(alleles[1],split="*",fixed=TRUE)[[1]][1]
  colnames(align)[2:(length(positions[[1]])+1)] <- positions[[1]]

  align[1,1] <- alleles[1]
  seqStr <- suppressMessages(AAsearch(alleles[1],positions[[1]],prefix=FALSE))
  align[1,2:(length(positions[[1]])+1)] <- if(!is.null(seqStr)) {strsplit(seqStr,split="~",fixed=TRUE)[[1]]} else {rep(" ",length(positions[[1]]))}

  for(n in 2:length(alleles)) {
    align[(n-1)*2,1] <- strsplit(alleles[n],split="*",fixed=TRUE)[[1]][1]
    align[(n-1)*2,2:(length(positions[[n]])+1)] <- positions[[n]]

    align[((n-1)*2)+1,1] <- alleles[n]
    seqStr <- suppressMessages(AAsearch(alleles[n],positions[[n]],prefix=FALSE))
    align[((n-1)*2)+1,2:(length(positions[[n]])+1)] <- if(!is.null(seqStr)) {strsplit(seqStr,split="~",fixed=TRUE)[[1]]} else {rep(" ",length(positions[[n]]))}

  }
  align
}

##########
##AAsearch
#'Search HLA alignments for specific alleles and amino acid positions.
#'
#'Searches for alleles in the ANHIG/IMGT HLA database, searches for amino acids for submitted allele name(s) and position(s).
#'
#'@param allelename A full or 2-field HLA allele name, excluding the "HLA-" prefix.
#'@param positions A vector of amino-acid positions (e.g., c(-17,1,22,130)); indel positions are named using decimals. So the first indel between positions 26 and 27 is named 26.1, and the second indel between 26 and 27 is 26.2, etc.
#'@param prefix A logical that indicates if the position number(s) should be included in the result.
#'@param sep The value that should be used to separate the sequences for each position. No value can be specified (sep=""). If prefix=FALSE and sep="", a sequence block is returned.
#'
#'@return A character string containing the corresponding amino acid sequence for each position in 'positions' for 'allelename'. An amino acid sequence is not returned if 'allelename' is not not found in the pertinent alignment. A position will be empty if 'allelename' does not have an amino acid at the specified position.
#'
#'@export
#'
#'@examples
#'AAsearch("DRB1*15:07:01",11:22)
#'AAsearch("DRB1*15:07",11:22)
#'
AAsearch <- function(allelename,positions,prefix=TRUE,sep="~"){

  trimmed <- FALSE

  # check required parameters
  if((missing(allelename) || missing(positions))) {stop("Please provide an allele name and at least one amino-acid position.\n")}

  # parse the allele name
  split <- unlist(gregexpr(pattern="\\*",allelename))
  locus <- substr(allelename,1,split-1)
  allele <- substr(allelename,split+1,nchar(allelename))

  # assess allelename validity
  if(!locus %in% names(HLAtools.data::HLAalignments$protAlignments) || length(split) != 1) {stop(allelename," is not a valid HLA allele name.")}

  # exclude aa positions that do not exist for a locus
  checkpos <- positions[positions %in% colnames(HLAtools.data::HLAalignments$protAlignments[[locus]]) == TRUE]
  posdiff <- length(positions)-length(checkpos)
  if(posdiff != 0) {message("There ", ifelse(posdiff==1,"is ","are "), xfun::numbers_to_words(posdiff), " listed position",ifelse(posdiff==1,"","s")," (",paste(positions[!positions %in% checkpos],collapse=","), ") ", ifelse(posdiff==1,"that does","that do")," not exist at the ", locus, " locus.",sep="")}
  if(length(checkpos) == 0) {return(message("There are no valid positions."))}
  positions <- unique(checkpos)

  # assess existence of allele names; check two-field column for 2-field truncates
  if(!allele %in% HLAtools.data::HLAalignments$protAlignments[[locus]]$allele) {
    message(allelename," is not a known full-length allele name.")
    if(sum(charToRaw(allelename) == charToRaw(":")) == 1) {message("Checking 2-field allele names.")
      if(allelename %in% HLAtools.data::HLAalignments$protAlignments[[locus]]$trimmed_allele) {
        trimmed <- TRUE
      } else {message(allelename," is not a known 2-field allele name.")
        return(cat(""))}
    } else {return(cat(""))}
  }

  # numeric sort of aa positions
  if(any(posSort(positions,locus) != positions)) {
    positions <- posSort(positions,locus)
    message("Sorting positions to reflect sequence order.",if(prefix==FALSE){" Setting 'prefix=TRUE' for clarity."})
    prefix=TRUE
    positions <- as.numeric(positions)
  }

  return(multipleAAsearch(locus=locus,allele=allele,positions=positions,prefix=prefix,sep=sep,trimmed=trimmed))
}

################
##numFields
#'Identifies the number of fields in a colon-delimited HLA allele name.
#'
#'Returns the number of fields in a colon-delimited HLA allele name. A value of 1 is returned for digit-delimited HLA allele names.
#'
#'@param allele A colon-delimited HLA allele name.
#'
#'@return The number of fields in the allele name.
#'
#'@export
#'
#'@examples
#'numFields("HLA-A*01:01")
#'numFields("DRB1*04:03:01")
#'numFields("13:02:01:01")
numFields <- function(allele) {
  length(strsplit(allele,":",fixed=TRUE)[[1]])
}

###########
##singleAAsearch
#'Search HLA amino acid sequences at single specified position for a specified HLA allele
#'
#'Generates a character string of the amino acid (AA)  sequence at the specified position for the specified HLA allele name.
#'
#'@param locus Specific locus of allele.
#'@param allele Name of allele to search for.
#'@param position Amino-acid position.
#'@param prefix A logical that indicates if the position number(s) should be included in the result.
#'@param trimmed A logical identifying whether 'allele' is a two-field, "trimmed" allele name (trimmed=TRUE), or a full-length name (trimmed=FALSE). The default value is FALSE.
#'
#'@return  Amino-acid (possibly with position number) at specified position for specified allele as available in the ANHIG/IMGTHLA Github Repository. position will be left empty if specific allele does not have a peptide at the input position.
#'
#'@note For internal HLAtools use only.
#'
#'@export
#'
#'@examples
#'singleAAsearch("DRB1", "15:07:01", 11)
#'
singleAAsearch <- function(locus, allele, position, prefix=TRUE,trimmed=FALSE){
  ## Initial parameter validations -- SJM
  if(length(position) > 1) {position <- position[1]; warning("More than one position was specified. Results will be returned for position ",position," only.")}
  if(trimmed == TRUE && numFields(allele) > 2) {return(warning("Two-field allele names are required when 'searchfield = trimmed_allele'."))}
  ## End -- SJM
  if(trimmed == TRUE) {
    aa_seq <- HLAtools.data::HLAalignments$protAlignments[[locus]][HLAtools.data::HLAalignments$protAlignments[[locus]]$trimmed_allele %in% paste(locus,allele,sep="*"),colnames(HLAtools.data::HLAalignments$protAlignments[[locus]]) %in% as.character(position)][1]
  }
  else {
    aa_seq <- HLAtools.data::HLAalignments$protAlignments[[locus]][HLAtools.data::HLAalignments$protAlignments[[locus]]$allele_name %in% paste(locus,allele,sep="*"),colnames(HLAtools.data::HLAalignments$protAlignments[[locus]]) %in% as.character(position)]
    ## Catching 'character(0)' results -- SJM
    if(identical(aa_seq, character(0))) {return(warning(paste("No value was found for ", allele,". This may not be a full allele name. Please redo your search using trimmed=TRUE.",sep="")))}
    ## End -- SJM
  }
  if(prefix){paste0(position,aa_seq)}else{aa_seq}
}

#############
##multipleAAsearch
#'Search HLA amino acid sequences at multiple positions for a specified HLA allele
#'
#'Generates a character string of amino acid (AA) sequences at the specified positions for the specified HLA allele name.
#'
#'@param locus Specific locus of allele.
#'@param allele Name of allele to search for.
#'@param positions Amino-acid positions to search.
#'@param prefix A logical that indicates if the position number(s) should be included in the result.
#'@param sep Defines the separator between position sequences
#'@param trimmed A logical identifying whether 'allele' is a two-field, "trimmed" allele name (trimmed=TRUE), or a full-length name (trimmed=FALSE). The default value is FALSE.
#'
#'@return  Amino-acids (possibly with position number) at specified positions for specified allele as available in the ANHIG/IMGTHLA Github Repository. position will be left empty if specific allele does not have a peptide at the input position.
#'
#'@note For internal HLAtools use only.
#'
#'@export
#'
#'@examples
#'multipleAAsearch("DRB1", "15:07:01", 11:22)
#'
multipleAAsearch <- function(locus, allele, positions, prefix=TRUE,sep="~", trimmed=FALSE){

  ## Initial parameter validations -- SJM
  if(length(allele) > 1) {allele <- allele[1]; warning("More than one allele was specified. Results will be returned for the first allele only.")}
  if(length(locus) > 1) {locus <- locus[1]; warning("More than one locus was specified. Results will be returned for the ", locus, " locus only.")}
  if(trimmed == TRUE && numFields(allele) > 2) {return(warning("Two-field allele names are required when 'trimmed=TRUE'."))}
  ## End -- SJM

  motif <- ""

  if(trimmed==TRUE) {
    fullname <- HLAtools.data::HLAalignments$protAlignments[[locus]][HLAtools.data::HLAalignments$protAlignments[[locus]]$trimmed_allele %in% paste(locus,allele,sep="*"),4][1]
    message("Returning sequence for ",fullname,", the first of ", nrow(HLAtools.data::HLAalignments$protAlignments[[locus]][HLAtools.data::HLAalignments$protAlignments[[locus]]$trimmed_allele %in% paste(locus,allele,sep="*"),])," ",paste(locus,allele,sep="*")," alleles.",sep="")
  }

  for(i in 1:length(positions)) {

    #motif <- paste(motif,paste0(as.character(positions[i]),singleAAsearch(locus,allele,positions[i],prefix)),sep="~")
    motif <- paste(motif,singleAAsearch(locus,allele,positions[i],prefix,trimmed=trimmed),sep=sep)
  }

  substr(motif,nchar(sep)+1,nchar(motif))
}

#########
##posSort
#'Correctly sorts a sequence alignment positions that contain indels.
#'
#'Correctly sorts a sequence alignment positions that contain indels identified as X.N, eg., Positions L, K, L, N, N.1, N.2 ... N.N, M, M, O, where N.1 - N.N are insertions between positions N and M.
#'
#'@param posVec a vector of amino acid positions.
#'@param locus a locus in HLAalignments$protalignments.
#'
#'@return A correctly sorted sequence.
#'
#'@export
#'
#'@examples
#'posSort(c(2,4,3,1,5), "DRB1")
#'
posSort <- function(posVec,locus){
    tab <- as.data.frame(cbind(match(posVec,names(HLAtools.data::HLAalignments$protAlignments[[locus]])),posVec))
    tab <- tab[order(tab$V1),]
    tab$posVec
}
