#### unified alignment search and construction functions v4.3.0 20 March 2024 Ryan Nickens & Steven Mack 

################
#AlignmentSearch
#'Search Alignments for Specific Positions in a Specific Allele
#'
#'Searches ANHIG/IMGT-HLA alignments and returns protein, codons or nucleotide sequences for a submitted allele name and position(s).
#'
#'@param alignType The type of alignment being searched. Allowed values are "codon","gen", nuc" and "prot". Only one 'alignType' value is allowed.
#'@param allelename A full or 2-field HLA allele name, excluding the "HLA-" prefix.
#'@param positions A vector of sequence positions (e.g., c(-17,1,22,130)); in nucleotide and genomic alignments, indel positions are named using decimals. So the first indel between positions 26 and 27 is named 26.1, and the second indel between 26 and 27 is 26.2, etc.
#'@param prefix A logical that indicates if the position number(s) should be included in the result.
#'@param sep The value that will be used to separate the sequences for each position. No value can be specified (sep=""). If prefix=FALSE and sep="", a sequence block is returned.
#'
#'@return A character string containing the corresponding peptide, codon or nucleotide sequence for each position in 'positions' for 'allelename'. A sequence is not returned if 'allelename' is not not found in the pertinent alignment. A position will be empty if 'allelename' does not have a value at the specified position.
#'
#'@importFrom xfun numbers_to_words
#'@export
#'
#'@note This function requires that the HLAalignments object has been populated with alignments via the alignmentFull() function.
#'@note Indel positions must be text-formatted (e.g. "607.12"). C.f., alignmentSearch("nuc","DRB1*15:07:01",11:22) vs alignmentSearch("nuc","DRB1*11:250N",c(605,"607.1","607.12",608)).
#'
alignmentSearch <- function(alignType,allelename,positions,prefix=TRUE,sep="~"){

  alignType <- checkAlignType(alignType)
  if(length(alignType)!=1) {stop(paste("Please specify only one 'alignType'."))}
  
  trimmed <- FALSE
  
  #make sure all necessary parameters are provided(allelename and positions)
      if((missing(allelename) || missing(positions))) {
        #message if one or both of the parameters are not provided
        stop("Please provide an allele name and at least one codon position.\n")
      }
    #splits allele at this point, two sections are subsequently assigned to variables for later use
    split <- unlist(gregexpr(pattern="\\*",allelename))
    #name of allele without numbers(ex DRB1*17:07:01 is just DRB1)
    locus <- substr(allelename,1,split-1)
    #opposite of previous (ex DRB1*17:07:01 is just 17:07:01)
    allele <- substr(allelename,split+1,nchar(allelename))
  
    # assess allelename validity(makes sure allele is named in data source)
     if(!locus %in% names(HLAalignments[[alignType]]) || length(split) != 1) {
    #if it is not, return message
    stop(allelename," is not a valid HLA allele name.")
  }
  
  # exclude aa positions that do not exist for a locus
  checkpos <- positions[positions %in% colnames(HLAalignments[[alignType]][[locus]]) == TRUE]
  posdiff <- length(positions)-length(checkpos)
  #if there is a difference between input positions and checked positions, do this
    if(posdiff != 0) {
      message("There ", ifelse(posdiff==1,"is ","are "), xfun::numbers_to_words(posdiff), " listed position",ifelse(posdiff==1,"","s")," (",paste(positions[!positions %in% checkpos],collapse=","), ") ", ifelse(posdiff==1,"that does","that do")," not exist at the ", locus, " locus for the ", alignType," alignment.",sep="")
     }
  #if no positions are found in check, do this
    if(length(checkpos) == 0) {
    return(message("There are no valid positions."))
      }
      positions <- unique(checkpos)
  
  # assess existence of allele names; check two-field column for 2-field truncates
    if(!allele %in% HLAalignments[[alignType]][[locus]]$allele) {
    message(allelename," is not a known full-length allele name.")
    if(sum(charToRaw(allelename) == charToRaw(":")) == 1) {message("Checking 2-field allele names.")
        if(allelename %in% HLAalignments[[alignType]][[locus]]$trimmed_allele) {
          
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
      if(any(posSort(positions,alignType,locus) != positions)) {
        positions <- posSort(positions,alignType,locus)
        message("Sorting positions to reflect sequence order.",if(prefix==FALSE){" Setting 'prefix=TRUE' for clarity."})
        prefix=TRUE
        #positions <- as.numeric(positions)
    }
  return(multiSearch(alignType=alignType,locus=locus,allele=allele,positions=positions,prefix=prefix,sep=sep,trimmed=trimmed))
}

################
##multiSearch
#'Search Alignment Sequences at Multiple Positions for a Specified Allele
#'
#'Generates a character string of multiple protein, codon or nucleotide or genomic sequences at the specified positions for the specified allele name.
#'
#'@param alignType The type of alignment being searched. Allowed values are "prot", codon", "nuc" and "gen". Only one 'alignType' value is allowed.
#'@param locus A specific locus.
#'@param allele An allele name.
#'@param positions Nucleotide positions to search.
#'@param prefix A logical that indicates if the position number(s) should be included in the result.
#'@param sep Defines the separator between position sequences.
#'@param trimmed A logical identifying whether 'allele' is a two-field, 'trimmed' allele name (trimmed=TRUE), or a full-length name (trimmed=FALSE). The default value is FALSE.
#'
#'@return codons (possibly with position number) at 'positions' for 'allele'. If 'allele' does not have a codon at 'position', '' will be returned.
#'
#'@note For internal HLAtools use only.
#'
#'@export
#'
multiSearch <- function(alignType, locus, allele, positions, prefix=TRUE,sep="~",trimmed=FALSE){

    alignType <- checkAlignType(alignType)
    if(length(alignType)!=1) {stop(paste("Please specify only one 'alignType'."))}
  
    if(validateAllele(paste(locus,allele,sep="*"))) {
  if(trimmed) {
      if( paste(locus,allele,sep="*") %in% HLAalignments[[alignType]][[locus]]$allele_name) {
                  message(paste(allele,"is a full length alllele name.",sep=" "))
                  trimmed <- FALSE }
      } else {
                if( paste(locus,allele,sep="*") %in% HLAalignments[[alignType]][[locus]]$trimmed_allele && !paste(locus,allele,sep="*") %in% HLAalignments[[alignType]][[locus]]$allele_name) {
                  trimmed <- TRUE
                }
      }
  } else {stop(paste(paste(locus,allele,sep="*"),"is not a full-length or trimmed allele name.",sep=" "))}
  
  motif <- ""
  if(trimmed == TRUE) {
    fullname <- HLAalignments[[alignType]][[locus]][HLAalignments[[alignType]][[locus]]$trimmed_allele %in% paste(locus,allele,sep="*"),4][1]
 
    message("Returning sequence for ",fullname,", the first of ", nrow(HLAalignments[[alignType]][[locus]][HLAalignments[[alignType]][[locus]]$trimmed_allele %in% paste(locus,allele,sep="*"),])," ",paste(locus,allele,sep="*")," alleles.",sep="")

    allele <- strsplit(fullname,"*",fixed = TRUE)[[1]][2]
    trimmed <- FALSE
    }
    #go through array
      for(i in 1:length(positions)) {
        #using single search to find each set of nucleotides at a position
        motif <- paste(motif,uniSearch(alignType,locus,allele,positions[i],prefix,trimmed=trimmed),sep=sep)
    #   motif <- paste(motif,unicDNAsearch(locus,allele,positions[i],prefix,trimmed=trimmed),sep=sep)
     }
    #separates letters w "~"
    substr(motif,nchar(sep)+1,nchar(motif))
}

################
##uniSearch
#'Search Sequences at a Single Position for an Allele
#'
#'Generates a character string of the peptide, codon or nucleotide sequence at the specified position for the specified HLA allele.
#'
#'@param alignType The type of alignment being searched. Allowed values are "prot", "codon", "nuc" and "gen". Only one 'alignType' value is allowed.
#'@param locus A specific HLA locus.
#'@param allele The name of the allele being searched.
#'@param position The specified position.
#'@param prefix A logical that indicates if the position number(s) should be included in the result.
#'@param trimmed A logical identifying whether 'allele' is a two-field, 'trimmed' allele name (trimmed=TRUE), or a full-length name (trimmed=FALSE). The default value is FALSE.
#'
#'@return The peptide residue, codon or nucleotide sequence at specified position for specified allele as available in the ANHIG/IMGTHLA Github Repository.
#'
#'@note For internal HLAtools use.
#'
#'@export
#'
uniSearch <- function(alignType, locus, allele, position, prefix=TRUE, trimmed=FALSE){

    alignType <- checkAlignType(alignType)
    if(length(alignType)!=1) {stop(paste("Please specify only one 'alignType'."))}
  
    sec <- ""
  ## Initial parameter validations -- SJM
  if(length(position) > 1) {position <- position[1]; warning("More than one position was specified. Results will be returned for position ",position," only.")}
  if(trimmed == TRUE && numFields(allele) > 2) {return(warning("Two-field allele names are required when 'trimmed = TRUE'."))}
  ## End -- SJM
  #if input says trimmed=true...
  #finding 3 nucleotides at specified position
  if(trimmed == TRUE) {
    if(alignType == "codon") {
      cDNA_seq <- HLAalignments[[alignType]][[locus]][HLAalignments[[alignType]][[locus]]$trimmed_allele %in% paste(locus,allele,sep="*"),colnames(HLAalignments[[alignType]][[locus]]) %in% as.character(position)][1,]
    } else {
      cDNA_seq <- HLAalignments[[alignType]][[locus]][HLAalignments[[alignType]][[locus]]$trimmed_allele %in% paste(locus,allele,sep="*"),colnames(HLAalignments[[alignType]][[locus]]) %in% as.character(position)][1]
    }
    for(i in 1:length(cDNA_seq)) {
      sec <- paste0(sec,cDNA_seq[[i]])
    }
    #only use the first allele with shortened name
    sec <- sec[1]
    #includes position number in front of nucleotides
    numsec <- paste0(position, sec) } else {
    #if input says allele_name...
    #finding 3 nucleotides at specified position.  ##### SM 2/16/2024 ADDED AS.DATA.FRAME BELOW FOR NUC ALIGNMENTS
    cDNA_seq <- as.data.frame(HLAalignments[[alignType]][[locus]][HLAalignments[[alignType]][[locus]]$allele_name %in% paste(locus,allele,sep="*"),colnames(HLAalignments[[alignType]][[locus]]) %in% as.character(position)])
    #if no codon is found at the position
        if(nrow(cDNA_seq) == 0) {
          return(warning(paste("No value was found for ", allele,". This may not be a full allele name. Please redo your search using trimmed=TRUE.",sep="")))
        }
      for(i in 1:length(cDNA_seq)) {
        sec <- paste0(sec,cDNA_seq[[i]])
      }
    numsec <- paste0(position, sec)
    }
  #prefix decides how much information to return
  if(prefix){numsec}else{sec}
}

#### Alignment Generation Functions

################
##customAlign
#'Generate a Customized Peptide, Codon or Nucleotide Sequence Alignment.
#'
#'Generates a peptide, codon, coding or genomic nucleotide alignment table for user-specified HLA alleles at user-specified positions.
#'
#'@param alignType The type of alignment being searched. Allowed values are "prot", codon", "nuc" and "gen". Only one 'alignType' value is allowed.
#'@param alleles A vector of un-prefixed HLA allele names.
#'@param positions Either a vector of variant positions, against which all loci will be aligned, or a list of vectors of nucleotide positions, exactly one vector for each allele, against which each corresponding allele will be aligned.
#'
#'@return A data frame of allele names and the corresponding nucleotide sequences for each desired nucleotide position. an error message is returned if input locus is not available in the ANHIG/IMGTHLA Github Repository.
#'
#'@note This function requires that the HLAalignments object has been populated with alignments via the alignmentFull() function. C.f., customAlign("codon",c("DRB1*01:01","DQB1*02:01","DPB1*01:01"),c(1,2,3,7,8,9,13,14,15)) and customAlign("codon",c("DPB1*01:01:01:01","DQA1*01:01:01:01","DQB1*05:01:01:01"),list(19:35,1:4,6:9)).
#'
#'@export
#'
customAlign <- function(alignType,alleles,positions){

  alignType <- checkAlignType(alignType)
  if(length(alignType)!=1) {stop(paste("Please specify only one 'alignType'."))}
  
  #makes sure no duplicate alleles will be present in table
    alleles <- unique(alleles)
      if(is.list(positions)) {
      #calls multicDNAalign if multiple sets of positions are given
      multiAlign(alignType,alleles,positions)
      }
      else {
      #calls unicDNAalign if only one set of positions are given
        uniAlign(alignType,alleles,positions)
      }
  }

################
##uniAlign
#'Generate an Alignment for Specific Alleles at Specific Positions
#'
#'Generates a peptide, codon, coding or genomic nucleotide alignment at a single set of positions for HLA alleles at one or more loci.
#'
#'@param alignType The type of alignment being searched. Allowed values are "prot", codon", "nuc" and "gen". Only one 'alignType' value is allowed.
#'@param alleles A vector of un-prefixed HLA locus names.
#'@param positions A vector of codon positions, against which all loci will be aligned.
#'
#'@return A data frame of allele names and the corresponding codon sequence for specified position. a codon sequence is not returned for a specific allele if input allele is not available in the ANHIG/IMGTHLA Github Repository. position will be left empty if specific allele does not have a codon at the input position.
#'
#'@note For internal HLAtools use.
#'
#'@export
#'
uniAlign <- function(alignType, alleles,positions){

  alignType <- checkAlignType(alignType)
  if(length(alignType)!=1) {stop(paste("Please specify only one 'alignType'."))} 

  #making an array with correct dimensions
  align <- data.frame(matrix(" ", nrow = length(alleles), ncol = length(positions)+1),stringsAsFactors = FALSE)
  #naming first column
  colnames(align)[1] <- "Allele"
  #naming remaining columns, running through positions array and naming each column (top row)
  colnames(align)[2:(length(positions)+1)] <- positions
  
  for(i in 1:length(alleles)){
    #setting the names for first column
    align[i,1] <- alleles[i]
    #getting nucleotides using cDNAsearch
    seqStr <- suppressMessages(alignmentSearch(alignType, alleles[i],positions,prefix=FALSE))
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
##multiAlign
#'Generate an Alignment for Specific Alleles at Different Positions
#'
#'Generates a peptide, codon, coding nucleotide or genomic alignment for HLA alleles allowing each allele to be aligned to a different set of positions.
#'
#'@param alignType The type of alignment being searched. Allowed values are "prot", codon", "nuc" and "gen".  Only one 'alignType' value is allowed.
#'@param alleles A vector of un-prefixed HLA locus names.
#'@param positions A list of vectors of nucleotide positions, exactly one vector for each allele, against which each corresponding allele will be aligned.
#'
#'@return A data frame of 'allele' and the corresponding nucleotide sequence for specified positions designated for an allele. a nucleotide sequence is not returned for a specific allele if input allele is not available in the ANHIG/IMGTHLA Github Repository. position will be left empty if specific allele does not have a nucleotide at the input position.
#'
#'@note For internal HLAtools use.
#'
#'@export
#'
multiAlign <- function(alignType,alleles,positions){

  alignType <- checkAlignType(alignType)
  if(length(alignType)!=1) {stop(paste("Please specify only one 'alignType'."))}
  
    #initial parameter check
      if(length(alleles) != length(positions)) {return(warning("The number of sets of positions must exactly match the number of alleles, please adjust input accordingly."))}
      #creating an array of correct dimensions
      #as many rows as 2*alleles-1, as many columns as positions
      align <- data.frame(matrix(" ", nrow = (length(alleles)*2)-1, ncol = max(lengths(positions))+1),stringsAsFactors = FALSE)
      #first column name is first locus (ex. DQA1*01:01 is column name as DQA1(it cuts off back part of name))
      colnames(align)[1] <- strsplit(alleles[1],split="*",fixed=TRUE)[[1]][1]
      #remaining columns are given positions for first allele(first set in position array)
      colnames(align)[2:(length(positions[[1]])+1)] <- positions[[1]]
  
  #first spot in first column is named the full first allele
      align[1,1] <- alleles[1]
      #using cDNAsearch to search for nucleotides for the first allele at the first set of positions
      seqStr <- suppressMessages(alignmentSearch(alignType,alleles[1],positions[[1]],prefix=FALSE))
      #filling in the remainder of row 1
      align[1,2:(length(positions[[1]])+1)] <- if(!is.null(seqStr)) {
      #if return is not empty, fill in row with the split return from cDNAsearch method
      strsplit(seqStr,split="~",fixed=TRUE)[[1]]
        } else {
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
              #using cDNAsearch to find nucleotides for specified positions
              seqStr <- suppressMessages(alignmentSearch(alignType,alleles[n],positions[[n]],prefix=FALSE))
              #filling in remainder of row
              align[((n-1)*2)+1,2:(length(positions[[n]])+1)] <- if(!is.null(seqStr)) {
              #if function returns something, splits it and fills in row
              strsplit(seqStr,split="~",fixed=TRUE)[[1]]
          } else { #if function returns nothing, fill with space
      rep(" ",length(positions[[n]]))
          }
    
        }
      align
}

