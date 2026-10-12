### Query and lookup functions for rel_dna_ser and rel_ser_ser v2.5.0 September 22, 2026 SJM

  ################
  ##queryRelDNASer
  #'Query the WMDASerologialRelationships$rel_dna_ser table values for a specific allele
  #'
  #'Identify the antigen, expert exception or HATS classification of a specific HLA allele 
  #'  
  #'@param allele A colon-delimited HLA allele name. When a truncated allele name is provided as the value for 'allele', values are returned for the lowest numbered full-length version of that allele present in the loaded *alleleListHistory* object.
  #'@param columns A character vector identifying the desired column(s), from among "USA", "PSA", "ASA", "Expert" and "HATS", in the rel_dna_ser table, which will be queried for the value of 'allele'. The default value, "all", returns values for all columns. See Notes for additional details. 
  #'@param collapse A boolean determining if multiple values for an allele in a column should be slash-delimited (TRUE), so that only one line is returned, or if lines for each individual value should be returned (FALSE). The default value is FALSE. 
  #'
  #'@return A one-row data.frame identifying the values of "USA", "PSA", "ASA", "Expert" and/or "HATS" for 'allele'.
  #'
  #'@note "USA" stands for 'Unambiguous Serological Antigen', "PSA" stands for 'Possible Serological Antigen', "ASA" stands for 'Assumed Serological Antigen', "Expert" stands for 'Expert assigned exceptions', "HATS" stands for 'HLA to Serotype'.
  #'@note Details of the rel_dna_ser table can be found at https://hla.alleles.org/pages/wmda/wmda/, https://hla.alleles.org/pages/antigens/hla_antigens/ and  https://github.com/ANHIG/IMGTHLA/tree/Latest/wmda#rel_dna_sertxt.
  #'@note Details of HATS can be found at https://github.com/kosoegawa/HATS. The "HATS" column is only available in IPD-IMGT/HLA Database release versions 3.64.0 and later.
  #'
  #'@export
  #'
  #'@examples
  #'\dontrun{queryRelDNASer("B*07:155")
  #'queryRelDNASer("B*07:155","HATS")
  #'queryRelDNASer("A*24:02")
  #'queryRelDNASer("A*68:259")
  #'queryRelDNASer("DRB1*11:16")}
  queryRelDNASer <- function(allele,columns = "all",collapse = F) {
      
      if(columns == "") {message("Please provide a value for 'columns'.")
                          return(NA)}
    
      noHATS <- FALSE
      if(as.numeric(squashVersion(WMDASerologicalRelationships$version)) < 3640) {noHATS <- TRUE}
      
      if(noHATS) { 
        allowedColumns <- c("USA","PSA","ASA","Expert","all") } else {
          allowedColumns <- c("USA","PSA","ASA","Expert","HATS","all")
        }
      
      if(missing(allele)) {message("Please provide an value for 'allele'.")
          return(NA)}
      
      if(!length(allele) == 1) {message("The value of 'allele' must be a single allele.")
          return(NA)}
      
      if(!columns %in% allowedColumns) {
        ifelse(length(columns) == 1,
        message("The column value '",paste(columns[which(!(columns %in% allowedColumns))],collapse = ", "),"' is not a valid rel_dna_ser column.",sep=""),
        message("The column values '",paste(columns[which(!(columns %in% allowedColumns))],collapse = ", "),"' are not valid rel_dna_ser columns.",sep=""))
        return(NA)
      }
      
      alleleOrig <- allele    
      allCol <- FALSE
                # validation by translation
                allele <- translateAllele(allele,alleleListHistory$Version$Version,alleleListHistory$Version$Version,expand = TRUE)
      
                if(is.na(allele)) {message("Allele ", alleleOrig," is not a valid allele name in release version ",alleleListHistory$Version$Version,".",sep="")
                                    return(NA)}
                
                if(grepl("/",allele,fixed=T)) {
                    allele <- strsplit(allele,"/",fixed=TRUE)[[1]][1]
                    }
      if(length(columns) == 1) {
                 if(columns == "all") {          
                 dataType <- switch(columns,
                         "USA" = "Unrelated Serological Antigen data ",
                         "PSA" = "Possible Serological Antigen data ",
                         "ASA" = "Assumed Serological Antigen data ",
                         "Expert" = "Expert assigned exceptions in search determinants of some registries ",
                         "HATS" = "HATS calculated assigned specificity data ")
                             
                    ## dataType is only pertinent when a 'columns' contains a single value  
                 
                 columns <- allowedColumns
                 
                 allCol <- T
                 
            } 
       }
      
           combTab <- expandAntigenTerms("rel_dna_ser")
    
           if(allele %in% combTab$fullAllele) {
    
                preResult <- combTab[combTab$fullAllele == allele,] ## this is the return for columns = "all"
                
                rownames(preResult) <- 1:nrow(preResult)
                
                hold.allele <- preResult$fullAllele
                
                preResult <- preResult[,4:ncol(preResult)] 
              
                    if(allCol == F) { preResult <- preResult[,colnames(preResult) %in% columns,drop = FALSE] }

                   if(length(preResult) == 1 && preResult[1,1] == "") {return(paste("No rel_dna_ser ",dataType,"are available for ",allele,".",sep=""))}

                    preResult <- cbind(as.data.frame(hold.allele),preResult)
                    colnames(preResult)[1] <- "Allele"
            
                ## collapse ambiguous antigens into a slash-delimited string 
                    preResult <- flattenTable(preResult,collapse)
                    
                  ## remove any duplicated rows  
                return(preResult[!duplicated(preResult),])      
            
          } else {
            
                message(paste(allele,"is not present in the rel_dna_ser table.")) ## should not ever see this, as validateXtranslate should catch it
                return(NA)  
            
            }
      }

  ################
  ##relDNASerLookup
  #'Lookup alleles that share WMDASerologialRelationships$rel_dna_ser table values 
  #'
  #'Identify all the alleles that share a specific antigen, expert exception or HATS classification in the WMDASerologialRelationships$rel_dna_ser table
  #'
  #'@param column A character string identifying the desired column ("USA", "PSA", "ASA", "Expert" or "HATS") in the rel_dna_ser table in which to lookup alleles that share 'term'.
  #'@param term A character string identifying the antigen, exception or specificity to be looked up in 'column'.
  #'@param locus A character string identifying a specific locus to return values for. If 'locus' is not specified or is set to 'all', all alleles with the specified term, across all loci, will be returned. The default value is "".
  #'@param string A boolean identifying the desired format of the returned object -- either a vector of alleles (string = FALSE) or a slash-delimited string of alleles (string = TRUE). The default value is FALSE.
  #'
  #'@return A vector or a slash delimited string of all of the alleles that share 'term' in 'column' for 'locus'.
  #'
  #'@export
  #'
  #'@note "USA" stands for 'Unambiguous Serological Antigen', "PSA" stands for 'Possible Serological Antigen', "ASA" stands for 'Assumed Serological Antigen', "Expert" stands for 'Expert assigned exceptions', "HATS" stands for 'HLA to Serotype'.
  #'@note Details of the rel_dna_ser table can be found at https://hla.alleles.org/pages/wmda/wmda/, https://hla.alleles.org/pages/antigens/hla_antigens/ and https://github.com/ANHIG/IMGTHLA/tree/Latest/wmda#rel_dna_sertxt.
  #'@note Details of HATS can be found at https://github.com/kosoegawa/HATS. The "HATS" column is only available in IPD-IMGT/HLA Database release versions 3.64.0 and later.
  #'
  #'@examples
  #'\dontrun{relDNASerLookup("HATS",3801)
  #'relDNASerLookup("USA",0)
  #'relDNASerLookup("USA",0,"A")
  #'}
  #'
  relDNASerLookup <- function(column,term,locus="",string=FALSE){
    
    noHATS <- FALSE
    if(as.numeric(squashVersion(WMDASerologicalRelationships$version)) < 3640) {noHATS <- TRUE}
    
    allowedColumns <- if(noHATS) { allowedColumns <- c("USA","PSA","ASA","Expert") 
                            } else { 
                              allowedColumns <- c("USA","PSA","ASA","Expert","HATS") }
    
    
            if(missing(column)){message("Please provide a value for 'column'.")
                              return(NA)}
            if(missing(term)){message("Please provide a valiue for 'term'.")
                              return(NA)}
    
            if(!locus %in% c("","all")) {
            if(!validateLocus(locus,"gDNA")) {
                                message("Locus ",locus," is not valid.")
                                return(NA)}
                 }

            if(!column %in% allowedColumns) {
                    message(paste(column,"is not a rel_dna_ser column."))
                    return(NA)}
    
            if(!string %in% c(TRUE,FALSE)) {
                    message("The value of the 'string' parameter must be either TRUE or FALSE.")
                    return(NA)
              }
    
            term <- as.character(term)
              combTab <- expandAntigenTerms("rel_dna_ser")    # searchable version of the table      
           
                repColumn <- column
                column <- switch(column,
                      "USA" = 4,
                      "PSA" = 5,
                      "ASA" = 6,
                      "Expert" = 7,
                      "HATS" = 8)
                
              
                   if(!locus %in% c("","all")) {  ## want all alleles for a specific locus
                     star <- ifelse(combTab[1,2] == "A*","*","") ## do we need to append an asterisk to the 'naked' locus name to search? 
                     combTab <- combTab[combTab$Locus == paste(locus,star,sep=""),]
                   } 
               
                      outSet <- combTab$fullAllele[combTab[,column] == term] ###
                      
                      if(length(outSet) == 0) {return(paste("No alleles have the '",term,"' term in the '",repColumn,"' column.",sep=""))}
                      
                        if(string){
                            outSet <- paste(outSet,collapse = "/")}
                      
                            return(outSet)
      }
       
  ## Pick a term at a gene and get the other two terms --parameters = locus, and column = antigen, split or associated
  
  #############
  ##reSerSerLookup
  #'Lookup a specific antigen, split antigen or associated antigen, and identify the values for the other two
  #'
  #'@param locus locus A character string identifying a specific locus to return values for. Only one locus can be specified.
  #'@param column A character string identifying the desired column ("Antigen", "Split", or "Associated") in the rel_ser_ser table in which to lookup 'term'.
  #'@param term A number or character string identifying the antigen, exception or specificity to be looked up in 'column'.
  #'@param collapse A boolean determining if multiple values in a column should be slash-delimited (TRUE), so that only one line is returned, or if lines for each individual value should be returned (FALSE). The default value is FALSE. 
  #'
  #'@export
  #'
  #'@note Details of the rel_ser_ser table can be found at https://hla.alleles.org/pages/wmda/wmda/, https://hla.alleles.org/pages/antigens/hla_antigens/ and and https://github.com/ANHIG/IMGTHLA/tree/Latest/wmda#rel_ser_sertxt.
  #'
  #'@examples
  #'relSerSerLookup("Cw","Antigen","7")
  #'relSerSerLookup("B","Antigen","7")
  #'relSerSerLookup("B","Associated",4005)
  #'
  relSerSerLookup <- function(locus,column,term,collapse=FALSE){
    
    RSSTab <- expandAntigenTerms("rel_ser_ser")
    
    if(!locus %in% RSSTab$Locus) {
      message(locus," is not a rel_ser_ser locus.")
      return(NA)}
    
    if(missing(column)){message("Please provide a value for 'column'.")
      return(NA)}
    
    if(missing(term)){message("Please provide a term to lookup.")
      return(NA)}
    
    if(!column %in% c("Antigen","Split","Associated")) {message("Please provide, ''Antigen'', ''Split'' or ''Associated' for the value of 'column'.")
      return(NA)}
    
    if(missing(term)){message("Please provide a valiue for 'term'.")
      return(NA)}
    
    term <- as.character(term)

          column <- switch(column,
                           "Antigen" = 2,
                           "Split" = 3,
                           "Associated" = 4)
          
          RSSTab <- RSSTab[RSSTab$Locus == locus,]
          RSSTab <- RSSTab[RSSTab[,column] == term,]
          
          if(nrow(RSSTab) == 0) {
                  column <- switch(column - 1, "Antigen","Split","Associated") ## switch back for error reporting
                  message(paste("No matches for term",term,"in the",column,"column for the",locus,"locus were found.",sep=" "))
                  return(NA)}
          
          rownames(RSSTab) <- c(1:nrow(RSSTab))
          
          ### insert collapse here
          ## collapse ambiguous antigens into a slash-delimited string 
          RSSTab <- flattenTable(RSSTab,collapse)
          
          return(RSSTab)
  }
  
  ################
  ##appendWMDA
  #'Increase Searchability of WMDASerologicalRelationships's rel_dna_ser table
  #'
  #'Append a column of full HLA allele names to a table in the WMDASerologicalRelationships$rel_dna_ser table format.
  #'
  #'@param rel_dna_ser A version of the 'rel_dna_ser' table.
  #'
  #'@return An eight-column data.frame in which the first column includes full HLA allele names constructed from the "Locus" and "Allele" columns of the WMDASerologicalRelationships$rel_dna_ser table.
  #'
  #'@note For internal HLAtools use.
  #'
  #'@export
  #'
  appendWMDA <- function(rel_dna_ser) {
    
    star <- ifelse(rel_dna_ser[1,1] == "A*","","*")
    
    fullAllele <- paste(rel_dna_ser$Locus,rel_dna_ser$Allele,sep=star)
    
    cbind(fullAllele,rel_dna_ser)
  }
  
  ################
  ##expandAntigenTerms
  #'Reformat the rel_ser_ser and rel_dna_ser tables to enable searches within slash-delimited strings of antigen terms.
  #'
  #'Expands the rows of the rel_ser_ser and rel_dna_ser tables that contain slash-delimited antigen terms to contain unique values for searching.
  #'
  #'@param rel.table A string identifying the WMDASerologicalRelationships table -- 'rel_dna_ser' or 'rel_ser_ser' to be expanded.
  #'
  #'@return A data.frame  in which all of the slash-delimited terms have been split into individual rows. The returned 'rel_dna_ser' table will include a 'locus*allele' column.
  #'   
  #'@note For internal HLAtools use.
  #'@note The first time that the 'rel_dna_ser' table is expanded for a given IPD-IMGT/HLA Database release, the process can take several seconds.
  #'   
  #'@export
  #'   
  expandAntigenTerms <- function(rel.table) { # rel.table <- "rel_ser_ser" or "rel_dna_ser" ## only parameter
      
      if(!rel.table %in% c("rel_ser_ser","rel_dna_ser")) {message("The value of 'rel.table' must be either 'rel_ser_ser' or 'rel_dna_ser'.")
                return(NA)}
    
      if(file.exists(file.path(tempdir(),paste(rel.table,"rda",sep=".")))) {
                
          load(file.path(tempdir(),paste(rel.table,"rda",sep="."))) ## Retrieve a pregenerated version from tempdir()
        
          if(exists("rel_ser_ser") && rel.table == "rel_ser_ser") {
                  if(rel_ser_ser$version == WMDASerologicalRelationships$version) {
                    return(rel_ser_ser$rel_ser_ser)
                            }
                         }
          if(exists("rel_dna_ser") && rel.table == "rel_dna_ser") {
                 if(rel_dna_ser$version == WMDASerologicalRelationships$version) {
                   return(rel_dna_ser$rel_dna_ser)
                      }              
                    }
                  }
      
      if(rel.table == "rel_dna_ser") {
            ser.tab <- WMDASerologicalRelationships$rel_dna_ser
            n.rep <- 5 # of columns that could be replaced 
            span <- 3:7 # column numbers of the replaceable columns
            subtract <- 2 # value to subtract from the total number of columns to ID the target column with the slashes
            n.col <- ncol(WMDASerologicalRelationships$rel_dna_ser) # total number of columns  
               }
      
        if(rel.table == "rel_ser_ser") {
             ser.tab <- WMDASerologicalRelationships$rel_ser_ser
             n.rep <- 3
             span <- 2:4
             subtract <- 1
             n.col <- ncol(WMDASerologicalRelationships$rel_ser_ser)
        }
                to.expand <- as.list(rep("",n.rep))  #  ID all of the cells with /
                names(to.expand) <- span
      
                    for(i in span){
                        to.expand[[i-subtract]] <- which(grepl("/", ser.tab[,i]) == T)
                    }
      
                for(i in 1:length(to.expand)){  ## these are the columns with sero data  
              
                    if(length(to.expand[[i]]) > 0){ # skip columns with no slashes
                  
                        cur.col <- as.numeric(names(to.expand)[[i]])  ## to.expand column X with slashes in it 
                      
                              for(j in 1:length(to.expand[[i]])) {
                        
                                  orig.row <- unname(unlist(ser.tab[to.expand[[i]][j],]))  ## original version of the jth row 
                                  cur.split <- strsplit(orig.row[cur.col],split = "/",fixed = T)[[1]] ## the N individual slash-delimited antigens
                                  
                                  for(k in 1:length(cur.split)) { ## put each new row at the end of the df
                                    
                                      orig.row[cur.col] <- cur.split[k]
                                      ser.tab <- rbind(ser.tab,orig.row)
                                    
                                      } # end k loop
                                  } # end j loop
                             } # set of columns with any slashes 
                         } # end to.expand i loop
                
                ## combine all of the columns in to.expand
                # and all of the rows that had slashes
                
                     ser.tab <- ser.tab[-unlist(to.expand),]
                     ser.tab <- ser.tab[order(ser.tab$Locus),]
               
                     if(rel.table == "rel_dna_ser") {ser.tab <- appendWMDA(ser.tab)} ## add a locus*allele names column
               
                      ## save the generated table and its version in the tempdir for reloading in the future 
                        version.n <- WMDASerologicalRelationships$version                     
                     
                        if(rel.table == "rel_dna_ser") {
                          outfile <- file.path(tempdir(), paste("rel_dna_ser","rda",sep=".")) ## naming the object
                          rel_dna_ser <- list(ser.tab,version.n)
                          names(rel_dna_ser) <- c(rel.table,"version")
                          save(rel_dna_ser,file=outfile)
                        }
                        
                        if(rel.table == "rel_ser_ser") {
                          outfile <- file.path(tempdir(), paste("rel_ser_ser","rda",sep=".")) ## naming the object
                          rel_ser_ser <- list(ser.tab,version.n)
                          names(rel_ser_ser) <- c(rel.table,"version")
                          save(rel_ser_ser,file=outfile)
                        }
                     
               return(ser.tab) 
  }
  
  ################
  ##flattenTable
  #'Reduce a mulirow data.frame of antigen data to a single row.  
  #'
  #'Combines all of the antigens in a single column into a single slash-delimited string. 
  #'
  #'@param antTab A multi-row data.frame of antigen data generated in the relSerSerLookup() or queryRelDNASer() functions.
  #'
  #'@param coll A boolean value determining if the data.frame should be flattened. There is no default value. 
  #'
  #'@return When coll == TRUE, a one-row data.frame containing slash-delimited antigen values is returned. When coll == FALSE, the antTab is returned unchanged.
  #'
  #'@note For internal HLAtools use. 
  #'
  #'@export
  #'
  flattenTable<- function(antTab,coll) {
  
          if(nrow(antTab) > 1 && coll == T) { 
                    boolCol <- !sapply(antTab, function(x) length(unique(x)) == 1)
  
                    if(any(boolCol == T)) {
    
                    antTab[1,which(boolCol == T)] <- paste(antTab[,boolCol == T],collapse="/")
    
                    antTab <- antTab[1,]
                    }
              }
          return(antTab)
          }
  