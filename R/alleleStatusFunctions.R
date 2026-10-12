### Allele Status Query & Lookup Functions v1.0 29SEP2026

  ####################
  #queryAlleleStatus
  #'Query the allele status data for a particular allele in the alleleStatus$AlleleStatus data object.
  #'
  #'Returns the allele status information for a specific allele. 
  #'
  #'@param allele A 2- to 4-field HLA allele name, excluding the "HLA-" prefix. If a truncated allele name is provided, values for the lowest-numbered full-length version of that allele will be returned.
  #'@param column A vector identifying the alleleStatus columns for which data should be returned. Valid column names are "Cells", "Groups", "Confirmed", "Start", "End", "Partial" and "Type". If "all" or "" is specified, data for all columns will be returned. The default value is 'all'.
  #'@param vector A logical value that identifies if a vector should be returned (TRUE), or if a one-line data.frame should be returned (FALSE).
  #'
  #'@return A data.frame or vector containing the corresponding values of 'column' for 'allele'.
  #'
  #'@export
  #'
  #'@examples
  #'queryAlleleStatus("DQB1*05:07",c("Cells","Type"))
  #'
  queryAlleleStatus <- function(allele, column="all",vector = FALSE){
    
          if(!allele %in% alleleStatus$AlleleStatus$Allele) {
               allele <- validationByTranslation(allele)
            }
              
             if(length(column) == 1 && column %in% c("all","")){
                    column <- c("Cells","Groups","Confirmed","Start","End","Partial","Type")
                    }
    
            if(!all( column %in% (colnames(alleleStatus$AlleleStatus[2:ncol(alleleStatus$AlleleStatus)])) )){
                message("Please select columns from among ", paste(colnames(alleleStatus$AlleleStatus)[2:ncol(alleleStatus$AlleleStatus)],sep=" ",collapse = ", "),".")
                  return(NA)
                  }
    
                 status <- alleleStatus$AlleleStatus[alleleStatus$AlleleStatus$Allele == allele,colnames(alleleStatus$AlleleStatus) %in% column]
                 rownames(status) <- allele
                 
                 if(vector) { return(append(rownames(status),unlist(status[1,],use.names = FALSE)))
                   }
                 
                 status
                   
  }
  
  ####################
  #multiQueryAlleleStatus
  #'Query the the allele status data for multiple alleles.
  #'
  #'Returns a table of allele status information for a set of specified alleles. 
  #'
  #'@param alleles A vector of 2- to 4-field HLA allele names, excluding the "HLA-" prefixes. If truncated allele names are provided, values for the lowest-numbered full-length version of each allele will be returned.
  #'@param column A vector identifying the 'alleleStatus' columns for which data should be returned. Valid column names are "Cells", "Groups", "Confirmed", "Start", "End", "Partial" and "Type". If "all" or "" is specified, data for all columns will be returned. The default value is 'all'.
  #'
  #'@return A table containing the corresponding values of 'column' for 'alleles'.
  #'
  #'@export
  #'
  #'@examples
  #'multiQueryAlleleStatus(c("DQB1*05:07","B*57:01"),c("Cells","Type"))
  #'multiQueryAlleleStatus(c("DRB1*08:07","DQA1*01:66:01:01"),c("Cells","Partial"))
  #'
  multiQueryAlleleStatus <- function(alleles,column) {
    
        alleles <- unique(alleles)
    
       if(length(column) == 1 && column %in% c("all","")) {column <- colnames(alleleStatus$AlleleStatus[2:ncol(alleleStatus$AlleleStatus)]) }
    
         outframe <- data.frame(matrix(NA, nrow = length(alleles), ncol = length(column)+1))
         
         colnames(outframe) <- c("Allele",column)
    
         for(i in 1:length(alleles)) {
          
          outframe[i,] <- queryAlleleStatus(alleles[i],column,TRUE)
 
        }  
    
        outframe
  }
  
  ####################
  #AlleleStatusLookup
  #'#'Lookup alleles that share alleleStatus$AlleleStatus table values 
  #'
  #'Identify the alleles that share values in a column of the alleleStatus$AlleleStatus table
  #'
  #'@param loci A character string identifying a specific loci to look up values of term in. If loci is "" or "all", values 'term' for all loci will be looked up. 
  #'@param term A character string identifying the antigen, exception or specificity to be looked up in 'column'.
  #'@param column A character string identifying the desired column ("Cells", "Groups", "Confirmed", "Start", "End", "Partial" or "Type") in the AlleleStatus table in which to look up alleles that share 'term'.
  #'@param vector A boolean identifying the desired format of the returned object -- either a vector of alleles (vector = TRUE) or a slash-delimited string of alleles (vector = FALSE). The default value is TRUE.
  #'
  #'@return A vector of individual alleles that have 'term' in 'column' for 'loci', or a slash-delimited string of all the alleles that have 'term' in 'column' for 'loci'. 
  #'
  #'@examples
  #'alleleStatusLookup(c("A","F"),20,"Cells")
  #'alleleStatusLookup("all",120,"Cells")
  #'alleleStatusLookup("all",30,"Cells",TRUE)
  #'alleleStatusLookup("all",20,"Cells",FALSE)
  #'@export
  #'
  alleleStatusLookup <- function(loci,term,column,vector=TRUE){
    
            loci <- as.character(loci)
            term <- as.character(term)
            column <- as.character(column)
            
            locus <- "" ## holder for the table of all loci in AlleleStatus
            
            AStable <- alleleStatus$AlleleStatus
            
            for(i in 1:nrow(AStable)){   ## add locus column
              locus[i] <- strsplit(AStable[i,1],"*",fixed=T)[[1]][1]
            }
            
            AStable <- cbind(locus,AStable) ## added locus column on the left for narrowing searches
            
            if("" %in% loci || "all" %in% loci) {
              loci <- unique(AStable$locus)
            }
            
            if(!all(loci %in% unique(AStable$locus))){ # Stop for non-HLA loci
              message("The ", loci[!loci %in% unique(AStable$locus)], ifelse(length(loci[!loci %in% unique(AStable$locus)] == 1)," locus is","loci are"), " not found in the alleleStatus table.")
                return(NA)
                 }
            
            AStable <- AStable[AStable$locus %in% loci,] ## only the loci in 'loci'

                  repColumn <- column
                  searchCol <- switch(column,
                      "Cells" = 3,
                      "Groups" = 4,
                      "Confirmed" = 5,
                      "Start" = 6,
                      "End" = 7,
                      "Partial" = 8,
                      "Type" = 9)
  
                  if(is.null(searchCol)) {
                    message("'",column, "' is not a valid value for 'column'.")
                    return(NA)
                  }
                  
                 intReturn <-  AStable[AStable[,searchCol] == term,2] ## generate a vector for the alleles that match 'term' in 'column'
             
                 if(vector) {
                   
                   return(intReturn)
                   
                        } else {
                           paste(intReturn,collapse = "/")
                 }
                 
             }   
               
  ####################
  #alleleStatusSummary
  #'Identify the unique values in the alleleStatus table
  #'
  #'Returns a sorted vector of unique values for each column of the alleleStatus$AlleleStatus data object.
  #'
  #'@param column A character string identifying the column ("Cells", "Groups", "Confirmed", "Start", "End", "Partial" or "Type") to query.
  #'
  #'@return A vector of the unique values in 'column'.
  #'
  #'@export
  #'
  #'@examples
  #'alleleStatusSummary("Start")
  #'
  alleleStatusSummary <- function(column) {
    
         searchCol <- switch(column,
                        "Cells" = 2,
                        "Groups" = 3,
                        "Confirmed" = 4,
                        "Start" = 5,
                        "End" = 6,
                        "Partial" = 7,
                        "Type" = 8)
    
                    uniCol <- sort(unique(alleleStatus$AlleleStatus[,searchCol]))
                    
                    if(searchCol %in% c(2,3,5,6)) {
                      return(sort(as.numeric(uniCol)))
                    }
    
            sort(uniCol)
  }           
                 
  ####################
  #validationByTranslation
  #'Validate an HLA allele name by attempting to translate it to its most recent version.
  #'
  #'Applies translateAllele() to convert an allele name to the most recent IPD-IMGT/HLA release version. When 'allele' is not valid, NA is returned.
  #'
  #'@param allele An HLA nomenclature version 3 HLA allele name of any length.
  #'@param expand A logical value that determines if the full length version of 'allele' (expand = TRUE) or the provided value of 'allele' should be returned when 'allele' is valid (expand = FALSE). The default value is TRUE.
  #'
  #'@return The original or full length version of 'allele', when 'allele' is valid, or NA, 'when' allele is not valid, with an accompanying message.
  #'
  #'@examples
  #'validationByTranslation("DRB1*08:03")
  #'validationByTranslation("A*02")
  #'validationByTranslation("A*02",FALSE)
  #'validationByTranslation("C*0")
  #'validationByTranslation("B*")
  #'
  #'@note The most recent IPD-IMGT/HLA release version is determined by the loaded version of the alleleListHistory data object. 
  #'
  #'@export
  #'
  validationByTranslation <- function(allele, expand = TRUE){
    
              if(length(strsplit(allele,split="*",fixed = T)[[1]]) == 1) {
                message(allele," is not a valid allele name.")
                return(NA)
              }  
              if(nchar(strsplit(allele,split="*",fixed = T)[[1]][2]) == 1) {
                message(allele," is not a valid allele name.")
                return(NA)
              }
    
             alleleOrig <- allele    
                  # validation by translation
                      allele <- suppressMessages(translateAllele(allele,alleleListHistory$Version$Version,alleleListHistory$Version$Version,expand = TRUE))
    
             if(is.na(allele)) {message("Allele ", alleleOrig," is not a valid allele name in release version ",alleleListHistory$Version$Version,".",sep="")
                         return(NA)}
                      
             ifelse(expand,strsplit(allele,"/",fixed = TRUE)[[1]][1],alleleOrig)
             
  }
  