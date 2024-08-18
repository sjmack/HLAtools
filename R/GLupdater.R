#GLupdater v1.0.0 17AUG24 R. Nickens & S.J. Mack

################
##updateGL
#'Update a GL String Code to a Specified IPD-IMGT/HLA Database Version.
#'
#'A quality control wrapper for GLupdate, which updates a GL String Code to a desired reference database version.
#'
#'@param GLStringCode A character string of HLA allele names and operators in GL String Code format, signifying their relation with one another and the associated IMGT/HLA Database release version.
#'@param Version A character string identifying of the desired IPD-IMGT/HLA Database release version to which the alleles should be updated, going back to version 1.05.0
#'@param expand A logical that indicates whether user would like to return all allele names that contain the input allele name (TRUE), or if only the direct HLA ID match should be returned (FALSE).
#'
#'@return A version of the input GL String code (in the form of a character string) updated to the desired version.
#'
#'@export
#'
#'@examples
#'updateGL("hla#1.05.0#HLA-DPA1*0106", "3.52.0")
#'updateGL("hla#3.36.0#HLA-B*15:35", "3.52.0")
#'updateGL("hla#3.45.0#HLA-DPA1*02:01:01:04", "3.52.0")
#'\donttest{
#'updateGL("hla#3.45.0#HLA-A*02:08", "3.52.0", expand = TRUE)
#'updateGL("hla#1.05.0#HLA-DPA1*0106", "3.52.0", expand = TRUE)
#'updateGL("hla#1.05.0#HLA-DPA1*0106", "2.27.0", expand = TRUE)
#'}
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
#'@references Mack et al. HLA 05 July 2023 https://doi.org/10.1111/tan.15145
updateGL <- function(GLStringCode, Version, expand = FALSE) {
  #making sure desired version is possible
  outpos <- GLV2(Version)
  if(substr(outpos[1],1,1) == "X") {
  #checking GL String code input, making changes if necessary
  GLString <- GLvalidate(GLStringCode)
  #if code is not stopped by previous two lines of code, call GLupdate
  if(!GLString == FALSE){
    return(GLupdate(GLString = GLString, Version = Version, expand = expand))
            } 
  } else {
    return(outpos) }
}

################
##multiUpdateGL
#'Update a column of GL String Code data to a desired IPD-IMGT/HLA Database version.
#'
#'Updates a column from a data frame in GL String Code format to a desired reference database version.
#'
#'@param GLstringArray An array of HLA allele names and operators in GL String code format identifying their relation with one another and the pertinent IPD-IMGT/HLA Database release version.
#'@param Version A character string identifying the desired version to which the alleles should be updated, going back to version 1.05.0.
#'@param expand A logical to determine whether to include only the direct HLA ID match, or all possible allele matches.
#'
#'@return A version of the input array of GL String Codes (in a data frame) updated to the desired version. NA values are returned in place of alleles that are not present in Version.
#'
#'@export
#'
#'@examples
#'\donttest{
#' # Example update of two GL Strings Codes containing truncated alleles from version 3.1.0 to 3.53.0.
#'GLSC.ex[[2]][1:2] 
#'multiUpdateGL(GLSC.ex[[2]][1:2], Version = "3.53.0")
#'}
#'
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
#'@references Mack et al. HLA 05 July 2023 https://doi.org/10.1111/tan.15145
multiUpdateGL <- function(GLstringArray, Version, expand = FALSE) {
  #making sure desired version is possible
      outpos <- GLV2(Version)
      names <- colnames(GLstringArray)
      retArray <- data.frame(y = 1:nrow(as.data.frame(GLstringArray)))
      colnames(retArray) <- colnames(as.data.frame(GLstringArray))
          
      if (is.null(names)) {
            GLstringArray <-as.data.frame(GLstringArray)
        }
          names <- colnames(GLstringArray)
  
          codeSplit <- as.data.frame(do.call(rbind,strsplit(GLstringArray[,1],"#",fixed = TRUE)))
          allVersions <- unique(codeSplit$V2)
          
          stringBlocks <- as.list(rep("",length(allVersions))) ## 1 list element for each release version 1-N
          names(stringBlocks) <- allVersions ## new; each list element is named for the version 
          uniqueVec <- vector()
  
          for(i in 1:length(allVersions)){ ## a vector of the release version names; each i is a relesae vesion
            
              stringBlocks[[i]] <- as.numeric(rownames(codeSplit[codeSplit$V2 == allVersions[i],]))
    
                for(j in 1:length(stringBlocks[[i]])) {
                                                                                                               ## vvvv 1 instead of j    
                  uniqueVec <- unique(c(uniqueVec,strsplit(codeSplit$V3[stringBlocks[[i]][j]],split="[+?~+/^\\|]")[[1]])) 
                }

              uniqueVec <- sort(uniqueVec)
              
              vers <- allVersions[i] ## fixing improperly structured vers values e.g., 3.1.0 instead of 3.01.0
              if(!validateVersion(vers)) {
                vers <- expandVersion(substr(GLV2(vers),2,stop = nchar(GLV2(vers))))
              }
              
              if(expand){ ## creating the lookup table for returned strings
                    lookUp <- as.data.frame(cbind(uniqueVec,paste("HLA-",suppressMessages(GIANT(uniqueVec,vers,Version)),sep=""))) ## silently call GIANT when expand = TRUE
                    } else {
                    lookUp <- as.data.frame(cbind(uniqueVec,paste("HLA-",GIANT(uniqueVec,vers,Version),sep=""))) 
                }
              colnames(lookUp) <- c("from","to")
              
              ### When expand = TRUE; apply queryRelease() to find all of the names that match in release X
              if(expand){
                    lookUp <- cbind(lookUp,as.data.frame(rep(NA,nrow(lookUp))))
                    colnames(lookUp)[3] <- "multi"
                    
                    for(j in 1:nrow(lookUp)) {
                      
                            only <- queryRelease(Version,gsub("HLA-","",lookUp$from[j],fixed = TRUE),TRUE) ## all of the string-matches in that release
                            locPre <- strsplit(gsub("HLA-","",lookUp$from[j],fixed = TRUE),split = "*",fixed = TRUE)[[1]][1] # target locus, for length and match
                            lookUp$multi[j] <- paste(only[substr(only,1,nchar(locPre)+1) == paste(locPre,"*",sep="")],collapse="/") # exclude non-locus matches
                            message("Returning all matches for ", strsplit(lookUp$from[j],"-")[[1]][2],".",sep="")
                            lookUp$multi[j] <- str_replace_all(lookUp$multi[j],"/",paste("/HLA-",sep=""))
                            lookUp$multi[j] <- paste("HLA-",lookUp$multi[j],sep="")
                    }
                }
              
             for(j in 1:length(stringBlocks[[i]])) { ## translating for each set of stringBlocks
                 for(k in 1:nrow(lookUp)) {
                                                                                                          # switch between multi or single lookup 
                        codeSplit$V3[stringBlocks[[i]][j]] <- gsub(pattern = lookUp$from[k],replacement = ifelse(expand, lookUp$multi[k],lookUp$to[k]), x = codeSplit$V3[stringBlocks[[i]][j]],fixed = TRUE)
                                }
                        codeSplit$V2[stringBlocks[[i]][j]] <- Version
  
                        retArray$GLstringArray[stringBlocks[[i]][j]] <- paste(unlist(codeSplit[stringBlocks[[i]][j],]),collapse="#") 
                      }
                 }
          
 return(retArray$GLstringArray)
  
}


################
##GLupdate
#'Update a GL String Code.
#'
#'Updates the allele names in a Genotype List (GL) String Code to the desired reference database version using the IPD-IMGT/HLA Database's Allele List History table.
#'
#'@param GLString A character string of HLA allele names and operators in GL String Code format signifying their relation with one another and the pertinent HLA Allele List version.
#'@param Version A character string of the desired version to which the alleles to be updated, going back to version 1.05.0.
#'@param expand A logical that determines whether to return only the direct HLA ID allele match or all possible HLA allele matches.
#'
#'@return An updated version the GL String Code input (in the form of a character string) updated to the input desired version.
#'
#'@note For internal use only.
#'
#'@importFrom stringr fixed str_replace_all
#'
#'@export
#'
#'@examples
#'GLupdate("hla#3.36.0#HLA-B*15:35", "3.52.0")
#'GLupdate("hla#3.45.0#HLA-A*02:08", "3.52.0")
#'GLupdate("hla#3.45.0#HLA-A*02:08", "3.52.0", expand = TRUE)
#'\donttest{
#'GLupdate("hla#3.1.0#HLA-A*02:01+HLA-A*01:01:01:01", Version = "3.53.0")
#'GLupdate("hla#3.44.0#HLA-DPA1*02:01:01:05", "3.45.0")
#'GLupdate("hla#3.45.0#HLA-DPA1*02:01:01:05", "3.46.0")
#'GLupdate("hla#3.37.0#HLA-A*01:02", "3.52.0", expand = TRUE)
#'}
#'
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
#'@references Mack et al. HLA 05 July 2023 https://doi.org/10.1111/tan.15145
#'@source https://github.com/ANHIG/IMGTHLA/blob/Latest/Allelelist_history.txt
#'
GLupdate <- function(GLString, Version, expand = FALSE) {
  #getting input and output versions formatted correctly
  inpos <- GLV(GLString)
  outpos <- GLV2(Version)
  #this is to determine if a select group of genes need to be treated differently when updating from very old versions to recent versions
  outnum <- gsub("X", "",outpos)
  #creating string to return back
  gBack <- GLString

  #creating list of all alleles
  gBack <- gsub("HLA-","",gBack)
  origV <- strsplit(GLString, "#")[[1]][2]
  nameList <- strsplit(gBack, "#")[[1]][3]
  nameList <- strsplit(nameList, "[.]|[/]|[?]|[+]|[|]|[~]|\\^")[[1]]
  for (i in 1:length(nameList[[1]])) {
    nameList[[1]][i] <- gsub("HLA-","",nameList[[1]][i])
  }
  #ensuring no replicates
  nameList <- unique(nameList)

  if(expand) {
    for (i in 1:length(nameList)) {
      #list to put all name matches in for an allele
      nameOp <- vector("list", 1)
      rTU <- ""
      backBack <- nameList[i]
      #length of namelist value, will be used to shorten items in AlleleListHistory and to decide how to compare value
      LX <- length(strsplit(as.character(nameList[i]), "[:]|[*]")[[1]])
      #specific case when updating very old versions to recent version, must manually add ":" to get matches
      hope <- strsplit(strsplit(nameList[i], "[*]")[[1]][2], "")[[1]]
      #version 2.28.0 was last version not to use ":"
      if (outnum > 2280) {
        if (LX ==2) {
          #formatting input correctly so it may accurately be compared to the new version
          if (length(hope)==4) {
            xxx <- paste0(paste0(paste0(hope[1], hope[2]), ":"), paste0(hope[3]), hope[4])
            yyy <- paste0(paste0(strsplit(nameList[i], "[*]")[[1]][1], "*"), xxx)
            nameList[i] <- yyy
          }
          if (length(hope)==5) {
            xxx <- paste0(paste0(paste0(hope[1], hope[2]), ":"), paste0(hope[3]), hope[4])
            yyy <- paste0(paste0(strsplit(nameList[i], "[*]")[[1]][1], "*"), xxx)
            yyy <- paste0(paste0(yyy, ":0"), hope[5])
            nameList[i] <- yyy
          }
          if (length(hope)==6) {
            xxx <- paste0(paste0(paste0(paste0(hope[1], hope[2]), ":"), paste0(paste0(hope[3], hope[4]), ":")), paste0(hope[5], hope[6]))
            yyy <- paste0(paste0(strsplit(nameList[i], "[*]")[[1]][1], "*"), xxx)
            nameList[i] <- yyy
          }
        }
      }
      LX <- length(strsplit(as.character(nameList[i]), "[:]|[*]")[[1]])
      #search everything in column for match
      for (x in 1:length(alleleListHistory$AlleleListHistory[[outpos]])) {
        #original value in box being searched, will be added to list if altered allelename matches with input nameslist
        YORG <- alleleListHistory$AlleleListHistory[[outpos]][x]
        #creating altered allelename for specific item in AlleleListHistory
        YWAY <- ""
        YWA <- strsplit(as.character(alleleListHistory$AlleleListHistory[[outpos]][x]), "[:]|[*]")[[1]]
        try2 <- ""
        #this will only be used in older versions without decimals
        if (LX == 2) {
          hope2 <- strsplit(YWA[2], "")[[1]]
          if (length(hope2) > length(hope)) {
            te <- ""
            for (w in 1:length(hope)) {
              te <- paste0(te, hope2[w])
            }
            YWA[2] <- te
          }
          YWAY <- paste0(paste0(YWA[1], "*"), YWA[2])
        }
        #most inputs will go through this route
        if (LX > 2) {
          for (n in 2:(LX-1)) {
            YWAY <- paste0(paste0(YWAY, YWA[n]), ":")
          }
          YWAY <- paste0(paste0(YWA[1], "*"), YWAY)
          YWAY <- paste0(YWAY, YWA[LX])
        }
        #if input allelename and altered allelename from AllelelistHistory match, add original name from AlleleListHistory to the list of matches
        if(!is.na(YWAY)&!is.na(nameList[i])) {
          if(YWAY == nameList[i]){
            nameOp[[1]] <- append(nameOp[[1]], paste0("HLA-", YORG))
          }
        }

      }
      nameOp <- unique(nameOp)

      if (length(nameOp[[1]]) > 1) {
        for (w in 1:(length(nameOp[[1]])-1)) {
          rTU <- paste0(rTU, paste0(nameOp[[1]][w], "/"))
        }
        rTU <- paste0(rTU, nameOp[[1]][length(nameOp[[1]])])
      }
      if (length(nameOp[[1]]) ==1) {
        rTU <- nameOp[[1]][1]
      }

      if(is.null(nameOp[[1]])) {
        message(nameList[i], " has not been located in desired GL String Version.")
        gBack <- str_replace_all(gBack, fixed(backBack), "NA")
      }
      #editing string to return
      if (!is.null(nameOp[[1]])) {
        message("Returning all matches for ", nameList[i],".",sep="")
        gBack <- str_replace_all(gBack, fixed(backBack), rTU)
      }

    }
    #changing version to reflect the update
    gBack <- gsub(origV,Version,gBack, fixed = TRUE) ### v1.1.4 adding fixed = TRUE fixes the issue
    gBack #### this should be a return 1.1.4
  } else {
    for (i in 1:length(nameList)) {
      if(nameList[i] %in% alleleListHistory$AlleleListHistory[[inpos]]){
        #finding which row the allele is in
        ytemp <- which(alleleListHistory$AlleleListHistory[[inpos]] == nameList[i])
        #getting value of allele from desired version column
        valtemp <- alleleListHistory$AlleleListHistory[ytemp,outpos]
        #avoiding issues with NA
        if(is.na(valtemp)|is.null(valtemp)) {
          message(nameList[i], " has been removed from the Allele List or has had its name changed.")
          gBack <- str_replace_all(gBack, fixed(nameList[i]), "NA")
        }
        #if valtemp is not NA, replace
        else {
          #replacing words on gBack string to return
          gBack <- str_replace_all(gBack, fixed(nameList[i]), paste0("HLA-", valtemp))
        }
      }  else {
        message(nameList[i], " is not present in IPD-IMGT/HLA Database Version ", strsplit(GLString, "#")[[1]][[2]], ".")
        tempFind <- paste0(c("hla", origV, nameList[i]), collapse = "#")
        suppressMessages(lastEffort <- GLupdate(tempFind, Version = origV, expand = TRUE))
        cOne <- strsplit(lastEffort, split = "#",fixed=FALSE)[[1]][3] ### 1.1.4 fix - explicitly set fixed = FALSE
        cOne <- strsplit(cOne, "[/]",fixed = FALSE)[[1]][1] ### 1.1.4 fix - explicitly set fixed = FALSE
        if (cOne != "NA") {
          message("Returning first match for ", nameList[i],".",sep="")
          Q <- paste0(c("hla", origV, cOne), collapse = "#")
          temp <- GLupdate(Q, Version = Version)
          temp <- strsplit(temp, split = "#",fixed=FALSE)[[1]][3] ### 1.1.4 fix - explicitly set fixed = FALSE
          gBack <- str_replace_all(gBack, fixed(nameList[i]), temp)
        }
        else { gBack <- str_replace_all(gBack, fixed(nameList[i]), "NA")
        }
      }
    }
    gBack <- gsub(origV,Version,gBack, fixed = TRUE) ### v1.1.4 adding fixed = TRUE fixes the issue
    gBack

  }

}
################
##GLV
#'Retrieve version from input GL String.
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
#'Format GL String Code version number.
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
#'Locate matches for an incomplete IPD-IMGT/HLA Database version.
#'
#' Uses the provided version information to locate and return possible matches for an incompletely defined IPD-IMGT/HLA Database release version.
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
#'Reintroduce version decimals.
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
#'Validates a GL String Code.
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
