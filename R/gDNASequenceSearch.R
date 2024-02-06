##gDNASequenceSearch v1.0.0 5FEB24

################
##gDNAsearch
#'Search genomic HLA alignments for specific alleles and nucleotide positions.
#'
#'Generates a character string of nucleotide variants at specified positions in the genomic (gDNA) sequence of an HLA allele, using IPD-IMGT/HLA genomic alignments in the HLAtools.data package.
#'
#'@param allelename A full-length or 2-field HLA allele name, excluding the "HLA-" prefix, of an HLA allele with a genomic alignment.
#'@param positions A vector of nucleotide positions (e.g., 'c(-17,1,22,130)' or '25:35') in the pertinent alignment; indel positions are identified using decimals. E.g., the first indel between positions 26 and 27 is named 26.1, and the second indel between 26 and 27 is 26.2, etc.
#'@param prefix A logical that indicates if the position number(s) should be included in the returned string.
#'@param sep The value that will be used to separate the sequences for each position. No value can be specified (sep=""). If prefix=FALSE and sep="", a sequence block is returned. The default value of 'sep' is "~". If a trimmed (two-field version of a 3- or 4- field name) is provided, the sequence of the lowest-numbered full-length name is returned.
#'
#'@return A character string containing the corresponding nucleotide sequence for each position in 'positions' for 'allelename'. A nucleotide sequence is not returned if 'allelename' is not not found in the pertinent alignment. A position will be empty if 'allelename' does not have a nucleotide at the specified position.
#'
#'@export
#'
#'@examples
#'gDNAsearch("DRB1*15:07:01",11:22)
#'gDNAsearch("DRB1*15:07:01",11:22,prefix = FALSE,sep="")
#'gDNAsearch("DRB1*15:07",11:22)
#'gDNAsearch("DRB1*15:07:01",c(321,321.1,322),prefix = FALSE,sep="")
gDNAsearch <- function(allelename,positions,prefix=TRUE,sep="~"){

  trimmed <- FALSE

  #make sure all necessary parameters are provided(allelename and positions)
  if((missing(allelename) || missing(positions))) {
    #message if one or both of the parameters are not provided
    stop("Please provide an allele name and at least one amino-acid position.\n")
  }
  #splits allele at this point, two sections are subsequently assigned to variables for later use
  split <- unlist(gregexpr(pattern="\\*",allelename))
  #name of allele without numbers(ex DRB1*17:07:01 is just DRB1)
  locus <- substr(allelename,1,split-1)
  #opposite of previous (ex DRB1*17:07:01 is just 17:07:01)
  allele <- substr(allelename,split+1,nchar(allelename))

  # assess allelename validity(makes sure allele is named in data source)
  if(!locus %in% names(HLAtools.data::HLAalignments$gDNAAlignments) || length(split) != 1) {
    #if it is not, return message
    stop(allelename," is not a valid HLA allele name.")
  }

  # exclude aa positions that do not exist for a locus
  checkpos <- positions[positions %in% colnames(HLAtools.data::HLAalignments$gDNAAlignments[[locus]]) == TRUE]
  posdiff <- length(positions)-length(checkpos)
  #if there is a difference between input positions and checked positions, do this
  if(posdiff != 0) {
    message("There ", ifelse(posdiff==1,"is ","are "), xfun::numbers_to_words(posdiff), " listed position",ifelse(posdiff==1,"","s")," (",paste(positions[!positions %in% checkpos],collapse=","), ") ", ifelse(posdiff==1,"that does","that do")," not exist at the ", locus, " locus.",sep="")
  }
  #if no positions are found in check, do this
  if(length(checkpos) == 0) {
    return(message("There are no valid positions."))
  }
  positions <- unique(checkpos)

  # assess existence of allele names; check two-field column for 2-field truncates
  if(!allele %in% HLAtools.data::HLAalignments$gDNAAlignments[[locus]]$allele) {
    message(allelename," is not a known full-length allele name.")
    if(sum(charToRaw(allelename) == charToRaw(":")) == 1) {message("Checking 2-field allele names.")
      if(allelename %in% HLAtools.data::HLAalignments$gDNAAlignments[[locus]]$trimmed_allele) {
        trimmed <- TRUE
      }
      else {message(allelename," is not a known 2-field allele name.")
        return(cat(""))
      }
    }
    else {
      return(cat(""))
    }
  }

  # numeric sort of nucleotide positions
  if(any(posSort(positions,locus) != positions)) {
    positions <- posSort(positions,locus)
    message("Sorting positions to reflect sequence order.",if(prefix==FALSE){" Setting 'prefix=TRUE' for clarity."})
    prefix=TRUE
    positions <- as.numeric(positions)
  }
  return(multigDNAsearch(locus=locus,allele=allele,positions=positions,prefix=prefix,sep=sep,trimmed=trimmed))
}

################
##unigDNAsearch
#'Search genomic HLA nucleotide sequences at single specified position for a specified HLA allele
#'
#'Generates a character string of the genomic (gDNA) nucleotide sequence at the specified position for the specified HLA allele name.
#'
#'@param locus An HLA locus in the HLAalignments$gDNAalignments object in the HLAtools.data package.
#'@param allele The name of an allele at 'locus'.
#'@param position A numeric value identifying a single nucleotide position. If more than one position is specified, results will be generated for the first position.
#'@param prefix A logical that indicates if the position number should be included in the result.
#'@param trimmed A logical identifying whether 'allele' is a two-field, "trimmed" allele name (trimmed=TRUE), or a full-length name (trimmed=FALSE). The default value is FALSE.
#'
#'@return The nucleotide sequence at the specified position for the specified allele. If sequence is not available for a position, no value will be returned.
#'
#'@note For internal HLAtools use only.
#'
#'@export
#'
#'@examples
#'unigDNAsearch("DRB1", "15:07", 57, trimmed=TRUE)
#'unigDNAsearch("DRB1", "15:07:01", 57)

unigDNAsearch <- function(locus, allele, position, prefix=TRUE,trimmed=FALSE){
  sec <- ""
  ## Initial parameter validations -- SJM
  if(length(allele) > 1 ) {allele <- allele[1]; warning("More than one allele was specified. Results will be returned for the ",paste(locus,allele,sep="*")," allele only.")}
  if(length(position) > 1) {position <- position[1]; warning("More than one position was specified. Results will be returned for position ",position," only.")}
  # if(!(searchfield %in% c("allele_name","trimmed_allele"))) {return(warning("The value of 'searchfield' must be 'allele_name' or 'trimmed_allele'."))}
  if(trimmed == TRUE && numFields(allele) > 2) {return(warning("Two-field allele names are required when 'trimmed = TRUE'."))}
           ## End -- SJM
  #if input says search trimmed allele...
  #finding 1 nucleotide at specified position
  if(trimmed == TRUE) {
    gDNA_seq <- HLAtools.data::HLAalignments$gDNAAlignments[[locus]][HLAtools.data::HLAalignments$gDNAAlignments[[locus]]$trimmed_allele %in% paste(locus,allele,sep="*"),colnames(HLAtools.data::HLAalignments$gDNAAlignments[[locus]]) %in% as.character(position)]
  ## trimmed_allele result validations -- SJM
    gDNA_seq <- unique(gDNA_seq)
    if(length(gDNA_seq) > 1) {return(warning(paste("There are ",length(gDNA_seq)," nucleotide variants for position ", position, " for the two-field allele name ", paste(locus,allele,sep="*"),". Please redo this search using 'trimmed=FALSE'.", sep="")))}
           ## End -- SJM
    sec <- paste0(sec,gDNA_seq)
    #includes position number in front of nucleotide
    numsec <- paste0(position, sec) #if input says allele_name...  #finding 1 nucleotide at specified position
    }
    else {
    gDNA_seq <- HLAtools.data::HLAalignments$gDNAAlignments[[locus]][HLAtools.data::HLAalignments$gDNAAlignments[[locus]]$allele_name %in% paste(locus,allele,sep="*"),colnames(HLAtools.data::HLAalignments$gDNAAlignments[[locus]]) %in% as.character(position)]
   ## Catching 'character(0)' results -- SJM
    if(identical(gDNA_seq, character(0))) {return(warning(paste("No value was found for ", allele,". This may not be a full allele name. Please redo your search using 'trimmed=TRUE'.",sep="")))}
          ## End -- SJM
    sec <- paste0(sec,gDNA_seq)
    numsec <- paste0(position, sec)
  }
  #prefix decides how much information to return
  if(prefix){numsec}else{sec}
}

################
##multigDNAsearch
#'Search genomic HLA nucleotide sequences at multiple positions for a specified HLA allele
#'
#'Generates a character string of genomic (gDNA) nucleotide sequences at the specified positions for the specified HLA allele name.
#'
#'@param locus An HLA locus in the HLAalignments$gDNAalignments object in the HLAtools.data package.
#'@param allele The name of an allele at 'locus'.
#'@param positions A set of consecutive (e.g., 1:10) or non consecutive (e.g., c(1,3,34,344)) nucleotide positions in the 'locus' alignment.
#'@param prefix A logical that indicates if the position number(s) should be included in the result.
#'@param sep Defines the separator between position sequences. The default separator is '~'.
#'@param trimmed A logical identifying whether 'allele' is a two-field, 'trimmed' allele name (trimmed=TRUE), or a full-length name (trimmed=FALSE). The default value is FALSE.
#'
#'@return nucleotides (possibly with position number) at 'positions' for 'allele'. If 'allele' does not have a nucleotide at 'position', '' will be returned.
#'
#'@note For internal HLAtools use only.
#'
#'@export
#'
#'@examples
#'multigDNAsearch("DRB1", "15:07",c(23:26), trimmed=TRUE)
#'multigDNAsearch("DRB1", "15:07:01",c(23:26))

multigDNAsearch <- function(locus, allele, positions, prefix=TRUE,sep="~",trimmed=FALSE){

  motif <- ""
  #go through array
  for(i in 1:length(positions)) {
    #using single search to find each set of nucleotides at a position
    motif <- paste(motif,unigDNAsearch(locus,allele,positions[i],prefix,trimmed=trimmed),sep=sep)
  }
  if(trimmed == TRUE) {
    fullname <- HLAtools.data::HLAalignments$gDNAAlignments[[locus]][HLAtools.data::HLAalignments$gDNAAlignments[[locus]]$trimmed_allele %in% paste(locus,allele,sep="*"),4][1]
   # print(fullname) # SJM I Turned This Off
    message("Returning sequence for ",fullname,", the first of ", nrow(HLAtools.data::HLAalignments$gDNAAlignments[[locus]][HLAtools.data::HLAalignments$gDNAAlignments[[locus]]$trimmed_allele %in% paste(locus,allele,sep="*"),])," ",paste(locus,allele,sep="*")," alleles.",sep="")}
  #separates letters w "~"
  substr(motif,nchar(sep)+1,nchar(motif))

}


################
##gDNAalign
#'Generate a genomic DNA sequence alignment.
#'
#'Generates an gDNA alignment table for user-specified HLA alleles at user-specified genomic nucleotide positions.
#'
#'@param alleles A vector of un-prefixed HLA locus names
#'@param positions Either a vector of nucleotide positions, against which all loci will be aligned, or a list of vectors of nucleotide positions, exactly one vector for each allele, against which each corresponding allele will be aligned.
#'
#'@return A data frame of allele names and the corresponding nucleotide sequences for each desired nucleotide position. an error message is returned if input loci is not available in the ANHIG/IMGTHLA Github Repository
#'
#'@export
#'
#'@examples
#'gDNAalign(c("DRB1*01:01","DQB1*02:01","DPB1*01:01"),c(1,2,3,7,8,9,13,14,15))
#'gDNAalign(c("DQA1*01:01:01:01","DQB1*05:01:01:01","DPB1*01:01:01:01"),list(32:58,33:59,31:57))

gDNAalign <- function(alleles,positions){
  #makes sure no duplicate alleles will be present in table
  alleles <- unique(alleles)
  if(is.list(positions)) {
    #calls multigDNAalign if multiple sets of positions are given
    multigDNAalign(alleles,positions)
  }
  else {
    #calls unigDNAalign if only one set of positions are given
    unigDNAalign(alleles,positions)
  }
}

################
##unigDNAalign
#'Generate a genomic alignment for specific HLA alleles at a set of positions
#'
#'Generates a genomic (gDNA) nucleotide alignment at a single set of nucleotide positions for HLA alleles at one or more loci.
#'
#'@param alleles A vector of un-prefixed HLA locus names
#'@param positions a vector of nucleotide positions, against which all loci will be aligned.
#'
#'@return A data frame of allele names and the corresponding nucleotide sequence for specified position. a nucleotide sequence is not returned for a specific allele if input allele is not available in the ANHIG/IMGTHLA Github Repository. position will be left empty if specific allele does not have a nucleotide at the input position.
#'
#'@note for internal use only
#'
#'@export
#'
#'@examples
#'unigDNAalign("DQA1*01:01:01:01",c(32:58))

unigDNAalign <- function(alleles,positions){

  #making an array with correct dimensions
  align <- data.frame(matrix(" ", nrow = length(alleles), ncol = length(positions)+1),stringsAsFactors = FALSE)
  #naming first column
  colnames(align)[1] <- "Allele"
  #naming remaining columns, running through positions array and naming each column (top row)
  colnames(align)[2:(length(positions)+1)] <- positions

  for(i in 1:length(alleles)){
    #setting the names for first column
    align[i,1] <- alleles[i]
    #getting nucleotides using gDNAsearch
    seqStr <- suppressMessages(gDNAsearch(alleles[i],positions,prefix=FALSE))
    #filling in the table
    #these positions will be filled if seqStr is not empty.
    align[i,2:(length(positions)+1)] <- if(!is.null(seqStr)) {
      strsplit(seqStr,split="~",fixed = TRUE)[[1]]
    }
    #if value is null, space fills spots
    else {
      rep(" ",length(positions))
    }

  }
  align
}

################
##multigDNAalign
#'Generate a genomic alignment for specific HLA alleles at different sets of positions.
#'
#'Generates a genomic (gDNA) nucleotide alignment for HLA alleles allowing each allele to be aligned to a different set of positions.
#'
#'@param alleles A vector of un-prefixed HLA allele names.
#'@param positions A list of vectors of nucleotide positions, exactly one vector for each allele, against which each corresponding allele will be aligned.
#'
#'@return A data frame of 'allele' and the corresponding nucleotide sequence for 'positions'. a nucleotide sequence is not returned for a specific allele if input allele is not available in the ANHIG/IMGTHLA Github Repository. position will be left empty if specific allele does not have a nucleotide at the input position.
#'
#'@note for internal use only
#'
#'@export
#'
#'@examples
#'multigDNAalign(c("DQA1*01:01:01:01", "DRB1*01:01:01:01"),list(c(25,46,50,78,88), c(26,47,51,79,89)))
#'multigDNAalign(c("DQA1*01:01:01:01", "DRB1*01:01:01:01"),list(32:58, 33:59))
#'multigDNAalign(c("DQA1*01:01:01:01", "DQA1*05:01:01:01", "DRB1*01:01:01:01", "DRB1*11:01:02:02"),list(32:58, 32:58, 33:59, 33:59))

multigDNAalign <- function(alleles,positions){
  #creating an array of correct dimensions
  #as many rows as 2*alleles-1, as many columns as positions
  align <- data.frame(matrix(" ", nrow = (length(alleles)*2)-1, ncol = max(lengths(positions))+1),stringsAsFactors = FALSE)
  #first column name is first locus (ex. DQA1*01:01 is column name as DQA1(it cuts off back part of name))
  colnames(align)[1] <- strsplit(alleles[1],split="*",fixed=TRUE)[[1]][1]
  #remaining columns are given positions for first allele(first set in position array)
  colnames(align)[2:(length(positions[[1]])+1)] <- positions[[1]]

  #first spot in first column is named the full first allele
  align[1,1] <- alleles[1]
  #using gDNAsearch to search for nucleotides for the first allele at the first set of positions
  seqStr <- suppressMessages(gDNAsearch(alleles[1],positions[[1]],prefix=FALSE))
  #filling in the remainder of row 1
  align[1,2:(length(positions[[1]])+1)] <- if(!is.null(seqStr)) {
    #if return is not empty, fill in row with the split return from gDNAsearch method
    strsplit(seqStr,split="~",fixed=TRUE)[[1]]
  }
  else {
    #if array is empty, fill with spaces
    rep(" ",length(positions[[1]]))
  }
  #filling in remaining rows(rows 2 and on)
  for(n in 2:length(alleles)) {
    #effectively the column name, same as above, will be in same row as positions
    align[(n-1)*2,1] <- strsplit(alleles[n],split="*",fixed=TRUE)[[1]][1]
    #filling in row with positions
    align[(n-1)*2,2:(length(positions[[n]])+1)] <- positions[[n]]
    #full allele name, will be in same row as nucleotides
    align[((n-1)*2)+1,1] <- alleles[n]
    #using gDNAsearch to find nucleotides for specified positions
    seqStr <- suppressMessages(gDNAsearch(alleles[n],positions[[n]],prefix=FALSE))
    #filling in remainder of row
    align[((n-1)*2)+1,2:(length(positions[[n]])+1)] <- if(!is.null(seqStr)) {
      #if function returns something, splits it and fills in row
      strsplit(seqStr,split="~",fixed=TRUE)[[1]]
    }
    #if function returns nothing, fill with space
    else {
      rep(" ",length(positions[[n]]))
    }

  }
  align
}
