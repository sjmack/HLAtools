#GLupdater v2.0.0 27FEB2025 R. Nickens & S.J. Mack

###############
##multiUpdateGL
#'Update a Vector of GL String Code Data to a Desired IPD-IMGT/HLA Database Version
#'
#'Updates a column from a data frame in GL String Code format to a desired reference database version.
#'
#'@param GLSCs A vector of GL String codes. 
#'@param to A character string identifying the IPD-IMGT/HLA Database release version to translate 'GLSCs' to. Values can range from version 1.05.0 to the loaded version of the alleleListHistory object.
#'@param expand A logical that indicates if the lowest-numbered truncated allele names that match truncated allele names in 'GLSCs' should be returned (expand = FALSE), or if a slash-delimited string of all matching full-length allele names should be returned (expand = TRUE). The default value is FALSE.
#'@param verbose A logical that indicates if messages regarding the update process should be sent to the console (TRUE) or not (FALSE). The default value is FALSE.
#'
#'@return A version GLSCs (in a data frame) updated to the desired version. NA values are returned in place of alleles that are not present in 'to'.
#'
#'@export
#'
#'@examples
#'multiUpdateGLs <- multiUpdateGL(GLSC.ex$GL.String.Code[1:50],"2.25.0")
#'multiUpdateGLs <- multiUpdateGL(GLSC.ex$GL.String.Code[1:50],"3.58.0")
multiUpdateGL <- function(GLSCs,to,expand=FALSE,verbose=FALSE){
  
  if(!checkVersion(to)) {stop(paste("Release version",to,"is not loaded in the HLAtools package. Please provide a loaded release version for 'to', or update the alleleListHistory object."))}
  if(!is.logical(expand)) {stop(paste(expand,"is not a valid value for 'expand'. Please provide a logical value."))}
  if(!is.logical(verbose)) {stop(paste(verbose,"is not a valid value for 'verbose'. Please provide a logical value."))}
  
  outVec <-rep(NA,length(GLSCs))
  
  for(i in 1:length(GLSCs)){
    
    outVec[i] <- updateGL(GLSCs[i],to,expand,verbose)
    
  }
  
  outVec
}        

###############
##multiTranslateGLstring
#'Translate a Vector of GL Strings to a Desired IPD-IMGT/HLA Database Version
#'
#'Translates a column from a data frame in GL String Code format to a desired reference database version.
#'
#'@param GLstrings A vector of GL Strings encoded in the same IPD-IMGT/HLA Database release version. 
#'@param from A character string identifying the IPD-IMGT/HLA Database release version of 'GLstrings'.
#'@param to A character string identifying the IPD-IMGT/HLA Database release version to translate 'GLstrings' to. Values can range from version 1.05.0 to the loaded version of the alleleListHistory object.
#'@param expand A logical that indicates if the lowest-numbered truncated allele names that match truncated allele names in 'GLstrings' should be returned (expand = FALSE), or if a slash-delimited string of all matching full-length allele names should be returned (expand = TRUE). The default value is FALSE.
#'@param verbose A logical that indicates if messages regarding the trnslation process should be sent to the console (TRUE) or not (FALSE). The default value is FALSE.
#'
#'@return A version GLSCs (in a data frame) updated to the desired version. NA values are returned in place of alleles that are not present in 'to'.
#'
#'@export
#'
#'@examples
#'multiGLs <- multiTranslateGLstring(GLstring.ex$Gl.String[1:50],"3.01.0","2.15.0")
#'multiGLs <- multiTranslateGLstring(GLstring.ex$Gl.String[1:50],"3.30.0","3.58.0")
multiTranslateGLstring <- function(GLstrings,from,to,expand=FALSE,verbose=FALSE){
  
  if(!checkVersion(to)) {stop(paste("Release version",to,"is not loaded in the HLAtools package. Please provide a loaded release version for 'to', or update the alleleListHistory object."))}
  if(!checkVersion(from)) {stop(paste("Release version",from,"is not loaded in the HLAtools package. Please provide a loaded release version for 'from', or update the alleleListHistory object."))}
  if(!is.logical(expand)) {stop(paste(expand,"is not a valid value for 'expand'. Please provide a logical value.\n"))}
  if(!is.logical(verbose)) {stop(paste(verbose,"is not a valid value for 'verbose'. Please provide a logical value.\n"))}
  
  outVec <-rep(NA,length(GLstrings))
  
  for(i in 1:length(GLstrings)){
    
    outVec[i] <- translateGLstring(GLstrings[i],from,to,expand,verbose)
    
  }
  
  outVec
}        

###############
##updateGL
#'Update a GL String Code to a Specified IPD-IMGT/HLA Database Version
#'
#'Updates the elements of a GLString Code across IPD-IMGT/HLA Database release versions. Truncated allele names can be expanded to the list of all allele names that contain the truncated name.
#'
#'@param GLSC A GL String Code containing HLA allele names. 
#'@param to A character string identifying the IPD-IMGT/HLA Database release version to translate 'GLSC' to. Values can range from version 1.05.0 to the loaded version of the alleleListHistory object.
#'@param expand A logical that indicates if the lowest-numbered truncated allele names that match truncated allele names in 'GLSC' should be returned (expand = FALSE), or if a slash-delimited string of all matching full-length allele names should be returned (expand = TRUE). The default value is FALSE.
#'@param verbose A logical that indicates if messages regarding the update process should be sent to the console (TRUE) or not (FALSE). The default value is FALSE.
#'
#'@return A GL String Code updated the 'to' release version.
#'
#'@export
#'
#'@examples
#'updateGL(GLSC.ex$GL.String.Code[1],"2.15.0",FALSE,FALSE)
#'updateGL(GLSC.ex$GL.String.Code[1],"2.15.0",TRUE,FALSE)
#'
updateGL <- function(GLSC,to,expand = FALSE,verbose = FALSE) {
  
  if(!checkVersion(to)) {stop(paste("Release version",to,"is not loaded in the HLAtools package. Please provide a loaded release version for 'to', or update the alleleListHistory object."))}
  if(!is.logical(expand)) {stop(paste(expand,"is not a valid value for 'expand'. Please provide a logical value."))}
  if(!is.logical(verbose)) {stop(paste(verbose,"is not a valid value for 'verbose'. Please provide a logical value."))}
  
  parts <- strsplit(GLSC,"#",TRUE)[[1]]
  fromVers <- parts[2] # -- version to translate from
  if(!checkVersion(fromVers)) {stop(paste("Release version",fromVers,"is not loaded in the HLAtools package. Please update the alleleListHistory object."))}
  
  inString <- parts[3] # -- GL string to translate to 
  outString <- rep(NA,length(inString)) ## -- the translated string 
  
  outString <- translateGLstring(inString,fromVers,to,expand,verbose)
  
  paste("hla",to,outString,sep="#")
  
}

###############
##GLupdate
#'Update a GL String Code to a Specified IPD-IMGT/HLA Database Version
#'
#'Updates the elements of a GLString Code across IPD-IMGT/HLA Database release versions. Truncated allele names can be expanded to the list of all allele names that contain the truncated name.
#'
#'@param GLSC A GL String Code containing HLA allele names. 
#'@param to A character string identifying the IPD-IMGT/HLA Database release version to translate 'GLSC' to. Values can range from version 1.05.0 to the loaded version of the alleleListHistory object.
#'@param expand A logical that indicates if the lowest-numbered truncated allele names that match truncated allele names in 'GLSC' should be returned (expand = FALSE), or if a slash-delimited string of all matching full-length allele names should be returned (expand = TRUE). The default value is FALSE.
#'@param verbose A logical that indicates if messages regarding the update process should be sent to the console (TRUE) or not (FALSE). The default value is FALSE.
#'
#'@export
#'
#'@return A GL String Code updated the 'to' release version.
#'
#'@note This function is maintained for compatibility with previous releases, and will be removed in a future release. HLAtools users are advised to use updateGL().
#'
#'@examples
#'GLupdate(GLSC.ex$GL.String.Code[1],"2.15.0",FALSE,FALSE)
#'GLupdate(GLSC.ex$GL.String.Code[1],"2.15.0",TRUE,FALSE)
#'
GLupdate <- function(GLSC,to,expand = FALSE,verbose = FALSE) {
  
  updateGL(GLSC,to,expand,verbose)
  
}
###############
##translateGLstring
#'Translate a GL String across IPD-IMGT/HLA Database Release Versions
#'
#'Translates a single Genotype List (GL) string across IPD-IMGT/HLA Database release versions. Truncated allele names can be expanded to the list of all allele names that contain the truncated name.  
#'
#'@param GLstring A GL String of HLA allele names. 
#'@param from A character string identifying the IPD-IMGT/HLA Database release version of 'GLstring'.
#'@param to A character string identifying the IPD-IMGT/HLA Database release version to translate 'GLstring' to. Values can range from version 1.05.0 to the loaded version of the alleleListHistory object.
#'@param expand A logical that indicates if the lowest-numbered truncated allele names that match truncated allele names in the GL String should be returned (expand = FALSE), or if a slash-delimited string of all matching full-length allele names should be returned (expand = TRUE). The default value is FALSE.
#'@param verbose A logical that indicates if messages regarding allele translations should be sent to the console (TRUE) or not (FALSE). The default value is FALSE.
#'
#'@export
#'
#'@return A GL String of allele names in the 'to' release version. If the GL String is determined to be invalid, FALSE is returned. 
#'
#'@examples
#'translateGLstring(GLstring.ex$Gl.String[1],"3.01.0","2.15.0")
#'translateGLstring(GLstring.ex$Gl.String[1],"3.01.0","2.15.0",TRUE)
#'translateGLstring(GLstring.ex$Gl.String[1],"3.01.0","2.15.0",FALSE,TRUE)
#'translateGLstring(GLstring.ex$Gl.String[1],"3.01.0","2.15.0",TRUE,TRUE)
#'
translateGLstring <- function(GLstring,from,to,expand=FALSE,verbose=FALSE) {      
  
  if(!checkVersion(to)) {stop(paste("Release version",to,"is not loaded in the HLAtools package. Please provide a loaded release version for 'to', or update the alleleListHistory object."))}
  if(!checkVersion(from)) {stop(paste("Release version",from,"is not loaded in the HLAtools package. Please provide a loaded release version for 'from', or update the alleleListHistory object."))}
  if(!is.logical(expand)) {stop(paste(expand,"is not a valid value for 'expand'. Please provide a logical value.\n"))}
  if(!is.logical(verbose)) {stop(paste(verbose,"is not a valid value for 'verbose'. Please provide a logical value.\n"))}
  
      if(verbose) {
        if(!validateGLstring(GLstring,"1.1")) {return(FALSE)}
        } else {
          if(suppressMessages(validateGLstring(GLstring,"1.1")) == FALSE) {return(FALSE)}
        }
  
        #set up all of the parts
          bits <- strsplit(GLstring,"[+,/,|,?,^,~]")[[1]] ## the alleles to be translated, in order of appearance 
          
            if(length(bits) == 1) {
              return(translateAllele(bits,from,to,expand,verbose))
                  }
          
          pieces <- rep(NA,length(bits)-1) ## the operators that delimit the alleles, in order.
          nuts <- ifelse(grep("HLA-",GLstring),TRUE,FALSE) ## determine if the returned vector elements should be prefixed with "HLA-".
          bolts <- 0 ## spacers for determining the nature and order of delimiters
          screws <- rep(NA,length(bits)) ## the new, translated allele/strings, corresponding to each bit
  
          for(i in 1:length(pieces)) { # determine the nature and order of delimiters
              pieces[i] <- substr(GLstring,bolts+nchar(bits[i])+1,bolts+nchar(bits[i])+1)
              bolts <- bolts+nchar(bits[i])+1
            }
           
            for(i in 1:length(bits)) { ## looping through each allele in the string
              screws[i] <- translateAllele(bits[i],from,to,expand,verbose)
                } 
          
            paste(paste(c(rbind(screws[1:length(screws)-1],pieces)),collapse=""),screws[length(screws)],sep="")
                
}

################
##translateAllele
#'Translate HLA Allele Names Across IPD-IMGT/HLA Database Release Versions
#'
#'Translates a single HLA allele name across IPD-IMGT/HLA Database release versions.  Truncated allele names in one version can be expanded to a list of all allele names that match the provided allele name in a chosen version.
#'
#'@param allele A full or truncated HLA allele name. 
#'@param from A character string identifying the IPD-IMGT/HLA Database release version of 'allele'.
#'@param to A character string identifying the IPD-IMGT/HLA Database release version to translate 'allele' to. Values can range from version 1.05.0 to the loaded version of the the alleleListHistory object. 
#'@param expand A logical that indicates if the lowest-numbered truncated allele name that matches a truncated allele name should be returned (expand = FALSE), or if a slash-delimited string of all matching full-length allele names should be returned (expand = TRUE). The default value is FALSE.
#'@param verbose A logical that indicates if messages regarding allele translations should be sent to the console (TRUE) or not (FALSE). The default value is FALSE.
#'
#'@return A character string of either a single matching allele or all matching alleles in the 'to' release version. If an allele does not have a cognate match in 'to', NA is returned for that allele. 
#'
#'@export
#'         
#'@examples
#'translateAllele("A*01:01","3.01.0","2.20.0",FALSE)
#'translateAllele("A*0101","2.20.0","3.01.0",TRUE)
#'translateAllele("B*57:01","3.12.0","3.12.0",TRUE)
#'
translateAllele <- function(allele,from,to,expand=FALSE,verbose = FALSE) {
  
  if(!checkVersion(to)) {stop(paste("Release version",to,"is not loaded in the HLAtools package. Please provide a loaded release version for 'to', or update the alleleListHistory object."))}
  if(!checkVersion(from)) {stop(paste("Release version",from,"is not loaded in the HLAtools package. Please provide a loaded release version for 'from', or update the alleleListHistory object."))}
  if(!is.logical(expand)) {stop(paste(expand,"is not a valid value for 'expand'. Please provide a logical value.\n"))}
  if(!is.logical(verbose)) {stop(paste(verbose,"is not a valid value for 'verbose'. Please provide a logical value.\n"))}
  
  nuts <- ifelse(substr(allele,1,4) == "HLA-", TRUE,FALSE) ## should the returned vector elements be prefixed with "HLA-".
  
  origField <- origDigit <- NA ## determine if 'allele' is field or digit delimited, and how many fields
  if(length(grep(":",allele)) == 0) {
    origDigit <- nchar(strsplit(allele,"*",fixed = TRUE)[[1]][2])  #vs 1 & 2 allele names
  } else { origField <- length(strsplit(allele,":",fixed=TRUE)[[1]]) # v3 allele names
  }
  
  #---- single allele translation
  loc <- strsplit(stringr::str_sub(allele,-(nchar(allele)-ifelse(nuts,4,0))),"*",fixed = TRUE)[[1]][1]
  check <- queryRelease(from,stringr::str_sub(allele,-(nchar(allele)-ifelse(nuts,4,0))),TRUE) # the set of all allele names matching the current allele in the current release
  check <- sort(check[substr(check,1,nchar(loc)) == loc]) ## exclude alleles from other loci, sort in numerical order
  ## exclude partial-match-loci. e.g, "DMA*01:01" when the allele is "A*01:01"
  
  if(length(check) == 0) {
    if(verbose) {message(paste("Allele",stringr::str_sub(allele,-(nchar(allele)-ifelse(nuts,4,0))),"is not found in release",paste(from,".",sep=""),sep=" "))}
    currFrom <- NA
  } else {
    if(length(check) > 1) {
      if(verbose) {message(paste(length(check), "alleles matching",paste(ifelse(nuts,"HLA-",""),stringr::str_sub(allele,-(nchar(allele)-ifelse(nuts,4,0))),sep=""),"are found in release",paste(from,".",sep=""),"Translation will be performed using the", paste(ifelse(nuts,"HLA-",""),check[1],sep=""),"allele.",sep=" "))}
      currFrom <- check[1]
    } else { ## there is only one allele
      currFrom <- check 
    }
  }
  
  if(is.na(currFrom)) {  # no currFrom -- return NA
    return(NA)
  } else { # there is a currFrom
    
    ## determine which row the current allele is in
    currRow <- match(currFrom,alleleListHistory$AlleleListHistory[,colnames(alleleListHistory$AlleleListHistory) == GLV2(from)])
    
    updAllele <- alleleListHistory$AlleleListHistory[currRow,GLV2(to)] ## only need the exact match 
    
    if(is.na(updAllele)) {
      if(verbose){paste("The",ifelse(nuts,paste("HLA",currFrom,sep="-"),currFrom),"allele has no cognate in release",paste(to,".",sep=""),sep=" ")}
      return(NA)
    } ## allele is not found in 'to'.
    
    updDigit <- updField <- NA ## determine if the updated allele is field or digit delimited
    
    if(length(grep(":",updAllele)) == 0) {
      updDigit <- nchar(strsplit(updAllele,"*",fixed = TRUE)[[1]][2])  #vs 1 & 2 allele names
    } else { 
      updField <- length(strsplit(updAllele,":",fixed=TRUE)[[1]]) # v3 allele names
    }
    
    if(!is.na(origField)) { ## original allele was colon-field delimited
      if(!is.na(updField)) { ## updated allele was colon-field delimited
        updAllele <- alleleTrim(allele = updAllele,resolution = origField,version = as.numeric(strsplit(to,split = ".",fixed = TRUE)[[1]][1]))
      } else { ## updated allele was digit-delimited
        updAllele <- alleleTrim(allele = updAllele,resolution = origField,version = as.numeric(strsplit(to,split = ".",fixed = TRUE)[[1]][1]))
      }
    } else { ## original allele was digit delimited
      if(!is.na(updField)) { ## updated allele was field-delimited
        updAllele <- alleleTrim(allele = updAllele,resolution = ceiling(origDigit/2),version = as.numeric(strsplit(to,split = ".",fixed = TRUE)[[1]][1]))
      } else { ## updated field was digit-delimited
        updAllele <- alleleTrim(allele = updAllele,resolution = ceiling(origDigit/2),version = as.numeric(strsplit(to,split = ".",fixed = TRUE)[[1]][1]))
      }
    }
    
    ## this single allele is what will be returned if expand = FALSE
    if(!expand) {return(ifelse(nuts,paste("HLA-",updAllele,sep=""),updAllele))}  
    
    trans <- updAllele  
    
    if(expand) { ## multiple alleles
      
      check <- queryRelease(to,trans,TRUE) # the set of all allele names matching the current allele in the current release
      loc <- strsplit(stringr::str_sub(check[1],-(nchar(check[1])-ifelse(nuts,4,0))),"*",fixed = TRUE)[[1]][1] ## need to do it again for 'C'vs'Cw'
      updAlleles <- sort(check[substr(check,1,nchar(loc)) == loc]) ## exclude alleles from other loci, sort in numerical order
      if(verbose) {message(paste(length(updAlleles), "version", to, "alleles match",paste(allele,".",sep=""),sep=" "))}
      
      if(nuts){ updAlleles <- paste("HLA-",updAlleles,sep="") }
      
      trans <- paste(updAlleles,collapse="/")
      
    } ## end of multiple allele expansion
  } ## end of 'is there a currFrom?'
  
  trans # return the translated allele
}
######################################################################################

################
##GLV
#'Retrieve Version from Input GL String
#'
#'Extracts the version data from an input GL String Code, or provides appropriate options if version input is not present in the IPD-IMGT/HLA Database.
#'
#'@param GLString A character string of HLA allele names and operators in GL String Code format, signifying their relation with one another and the pertinent IPD-IMGT/HLA Database release version.
#'
#'@return Returns A character string of the version of the release as it appears in the pertinent alleleListHistory$AlleleListHistory-column-header.
#'
#'@export
#'
#'@note For internal use only.
#'
#'@examples
#'GLV("hla#3.25.0#HLA-B15:35")
#'
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
#'@references Mack et al. HLA 05 July 2023 https://doi.org/10.1111/tan.15145
#'
GLV <- function(GLString) {
  ve1 <-""
  ve2 <- ""
  version <- strsplit(strsplit(GLString,"#")[[1]][2], "[.]")[[1]]
  long <- strsplit(gsub("[.]", "", strsplit(GLString,"#")[[1]][2]), "")
  for(i in 1:length(version)) {
    ve1 <- paste0(ve1,version[[i]])
  }
  ve1 <- paste0("X",ve1)
  if (length(long[[1]]) == 3) {
    ve2 <- strsplit(ve1, "")[[1]]
    ve2 <- append(ve2, "0", 2)
    ve2 <- paste0(ve2, collapse = "")
    ve1 <- ve2
  }
  if (!ve1 %in% names(alleleListHistory$AlleleListHistory)) {
    ve1 <- redec(ve1)
    GLVhelper(ve1)
  }
  else {
    ve1
  }
}


################
##GLV2
#'Format GL String Code Version Number
#'
#'Returns a compressed IPD-IMGT/HLA Database release version or a vector of dot-delimited release version options if the specified version is not present in the alleleListHistory object.
#'
#'@param Version A character string of the desired IPD-IMGT/HLA Database version, going back to version 1.05.0.
#'
#'@return Returns either a single, compressed character string-formatted release version, or a character vector of potential release versions.
#'
#'@export
#'
#'@note For internal use only.
#'
#'@examples
#'GLV2("3.34.0")
#'GLV2("3.0.0")
#'
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
#'@references Mack et al. HLA 05 July 2023 https://doi.org/10.1111/tan.15145
#
GLV2 <- function(Version) {
  ve1 <- ve2 <- ""
  version <- strsplit(Version, "[.]")[[1]]
  long <- strsplit(gsub("[.]", "", Version), "")
  for(i in 1:length(version)) {
    ve1 <- paste0(ve1,version[[i]])
  }
  ve1 <- paste0("X",ve1)
  if (length(long[[1]]) == 3) {
    ve2 <- strsplit(ve1, "")[[1]]
    ve2 <- append(ve2, "0", 2)
    ve2 <- paste0(ve2, collapse = "")
    ve1 <- ve2
  }
  
  if (!ve1 %in% names(alleleListHistory$AlleleListHistory)) {
    GLVhelper(Version)
  } else {
    ve1
  }
}

################
##GLVhelper
#'Locate Matches for an Incomplete IPD-IMGT/HLA Database Version
#'
#'Uses the provided version information to locate and return possible matches for an incompletely defined IPD-IMGT/HLA Database release version.
#'
#'@param Version A character string of the desired version to which the alleles should be updated, going back to version 1.05.0.
#'
#'@return A list of character strings of possible matches to a given incomplete input.
#'
#'@export
#'
#'@note For internal use only.
#'
#'@examples
#'GLVhelper("2.25")
#'GLVhelper("3.9.0")
#'
#'
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
#'@references Mack et al. HLA 05 July 2023 https://doi.org/10.1111/tan.15145
#'
GLVhelper <- function(Version) {
  ve5 <-ve4 <-ve3 <-ve2 <- ""
  vte1 <- vte3 <- vector("list", 1)
  fin <- FALSE
  #making list of only the numbers
  version <- strsplit(gsub("[.]","",Version), "")[[1]]
  #check all names second letter, if its same as input, format with decimals and return options to choose from
  if (length(version) == 1) {
    for (i in 2:(length(names(alleleListHistory$AlleleListHistory)))) {
      if(strsplit(names(alleleListHistory$AlleleListHistory)[i], "")[[1]][2] == version) {
        vte1[[1]] <- append(vte1[[1]], names(alleleListHistory$AlleleListHistory)[i])
      }
    }
    if (!is.null(vte1[[1]])) {
      message("Returning all version ", version, " releases. Please re-enter one of these options.")
      for (i in 1:length(vte1[[1]])) {
        vte1[[1]][i] <- redec(vte1[[1]][i])
      }
      return(vte1[[1]])
      fin <- TRUE
    }  else {
      message(version, " is not a valid major release version.")
      fin <- TRUE
    }
  }
  #add 1 zero to the end, test for multiple other options
  if (length(version) == 2) {
    message("Version is truncated. Searching for matches by appending '0'.")
 
    ve2 <- redec(paste0("X",paste(append(version, "0"),collapse="")))

     return(GLVhelper(ve2))
     fin <- TRUE
  }
  if (length(version) == 3) {
    tte3 <- vector("list", 1)

    vv3 <- paste(version,collapse="")
    for(i in 1:length(version)) {
      xt <- append(version, "0", i)
      for (n in 1:length(xt)) {
        ve3 <- paste0(ve3,xt[n])
      }
      vte3[[1]][i] <- paste0("X",ve3)
      ve3 <- ""
    }
    for (i in 1:length(vte3[[1]])) {
      if(vte3[[1]][i] %in% names(alleleListHistory$AlleleListHistory)) {
        #message about what it is returning (its returning earliest version of input)
        tte3[[1]] <- append(tte3[[1]], vte3[[1]][i])
      }
    } 
    if(!is.null(tte3[[1]])) { # added 7/7
    for (i in 2:(length(names(alleleListHistory$AlleleListHistory)))) {
      yy <-names(alleleListHistory$AlleleListHistory)[i]
      yy <- substr(yy, start = 2, stop = 4)
      if(yy == vv3) {
        tte3[[1]] <- append(tte3[[1]], names(alleleListHistory$AlleleListHistory)[i])
          }
        }
    } else { ## when TTE is null, then none of the options exist
      
      message("No direct matches were found for release version '", Version, "', please try one of these options:")
      return(as.vector(sapply(names(alleleListHistory$AlleleListHistory)[substr(names(alleleListHistory$AlleleListHistory)[1:length(alleleListHistory$AlleleListHistory)],2,2) == version[1]],redec)))
         }
    
    if (!is.null(tte3[[1]])) {

      tte3[[1]] <- unique(tte3[[1]])
      for (i in 1:length(tte3[[1]])) {
        tte3[[1]][i] <- redec(tte3[[1]][i])
      }
      message("No direct matches were found for release version '", vv3, "', please try one of these options:")
      return(tte3[[1]])

      fin <- TRUE
    } else {    
      return(GLVhelper(version[1]))
      fin <- TRUE
    }
  }
  
  if(!fin) {
  #remove last element from list, and see if that works
  if (length(version) > 3) {
    
    xnames <- as.vector(sapply(names(alleleListHistory$AlleleListHistory)[-1],redec,simplify = TRUE))
    
    if(Version %in% xnames) {
      
      return(Version)
      
    } else {
      message("No direct matches were found for release version '", Version, "', please try one of these options:")

         return(GLVhelper(version[1]))
              }
          }
      }
  }

################
##redec
#'Reintroduce Version Decimals
#'
#'Correctly put decimals back into an AlleleListHistory column-formatted version name.
#'
#'@param Cname A character string describing a version name in AlleleListHistory column format.
#'
#'@return A character string describing a version name with decimals in appropriate places according to the IPD-IMGT/HLA Database.
#'
#'@export
#'
#'@note For internal use only.
#'
#'@examples
#'redec("X3090")
#'
redec <- function(Cname) {
  Rname <- ""
  Cname <- gsub('^.', '', Cname)
  Dname <- strsplit(Cname, "")[[1]]
  if(length(Dname)==4) {
    Dname <- append(append(Dname, ".", 1),".", 4)
      }
  if(length(Dname)==3) {
    Dname <- append(append(Dname, ".", 1),".", 3)
      }
  if(length(Dname)==2) {
    Dname <- append(Dname, ".", 1)
      }
  for(i in 1:length(Dname)) {
    Rname <- paste0(Rname,Dname[i])
      }
  Rname
}

################
##GLvalidate
#'Validate a GL String Code
#'
#'A lightweight validator that inspects a GL String Code for correct structure. If the namespace field contains a value other than value of 'namespace', the namespace is changed to 'hla'. The version and GL String fields are not evaluated.
#'
#'@param GLString A character string describing a string formatted to GL String Code specifications.
#'@param namespace A character vector of allowed namespace strings. The default value is 'c("hla","kir")'.
#'
#'@return A character string describing wither a well-formatted GL String Code or the value FALSE.
#'
#'@note For internal use only.
#'
#'@export
#'
#'@examples
#'GLvalidate("ha#3.25.0#hla-B15:35") 
#'
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
#'@references Mack et al. HLA 05 July 2023 https://doi.org/10.1111/tan.15145
#'
GLvalidate <- function(GLString,namespace = c("hla","kir")) {
  GLcopy <- GLString
  GLret <- ""
  if (length(strsplit(GLString, "#")[[1]]) > 3) {
    message("The input GL String code is formatted incorrectly, there are too many '#' delimited sections.")
    return(FALSE)
      }
      if (length(strsplit(GLString, "#")[[1]]) < 3) {
        message("The input GL String code is formatted incorrectly, there are too few '#' delimited sections.")
        return(FALSE)
      }
  if (length(strsplit(GLcopy, "#")[[1]]) == 3) {
    whole <- strsplit(GLcopy, '#')[[1]]
    first <- strsplit(GLcopy, '#')[[1]][1]
    se <- GLV(GLString)
    third <- strsplit(GLcopy, "#")[[1]][3]
        if (!first %in% namespace) {
          whole[1] <- "hla"
        }
    third <- gsub("hla-", "HLA-", third)
        for (i in 1:2) {
          GLret <- paste0(GLret, paste0(whole[i], "#"))
        }
    GLret <- paste0(GLret, third)

      return(GLret)
      }
  }
