#GLupdater v0.4.0 7FEB24 R. Nickens

#library('dplyr')
#library('stringr')

################
##updateGL
#'Update a GL String Code to a specified IPD-IMGT/HLA Database version.
#'
#'A quality control wrapper for GLupdate, which updates a GL String Code to a desired reference database version.
#'
#'@param GLStringCode A string of HLA allele names and operators in GL String Code format, signifying their relation with one another and the associated IMGT/HLA Database release version.
#'@param Version A string identifying of the desired IPD-IMGT/HLA Database release version to which the alleles should be updated.
#'@param expand A logical that indicates whether user would like to return all allele names that contain the input allele name (TRUE), or if only the direct HLA ID match should be returned (FALSE).
#'
#'@return A version of the input GL String code (in the form of a string) updated to the desired version.
#'
#'
#'@export
#'
#'@examples
#'\dontrun{
#'updateGL("hla#3.36.0#HLA-B*15:35", "3.52.0")
#'
#'updateGL("hla#3.45.0#HLA-DPA1*02:01:01:05", "3.52.0")
#'
#'updateGL("hla#3.45.0#HLA-A*02:08", "3.52.0", expand = TRUE)
#'
#'updateGL("hla#1.05.0#HLA-DPA1*0106", "3.52.0", expand = TRUE)
#'updateGL("hla#1.05.0#HLA-DPA1*0106", "2.27.0", expand = TRUE)
#'updateGL("hla#1.05.0#HLA-DPA1*0106", "3.52.0")
#'}
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
#'@references Mack et al. HLA 05 July 2023 https://doi.org/10.1111/tan.15145
updateGL <- function(GLStringCode, Version, expand = FALSE) {
  #making sure desired version is possible
  outpos <- GLV2(Version)
  #checking GL String code input, making changes if necessary
  GLString <- GLvalidate(GLStringCode)
  #if code is not stopped by previous two lines of code, call GLupdate
  GLupdate(GLString = GLStringCode, Version = Version, expand = expand)
}

################
##multiUpdateGL
#'Update columns of GL String Code data to a desired IPD-IMGT/HLA Database version.
#'
#'Updates columns from a data frame in GL String Code format to a desired reference database version.
#'
#'@param GLstringArray An array of HLA allele names and operators in GL String code format identifying their relation with one another and the pertinent IPD-IMGT/HLA Database release version.
#'@param Version A text string identifying the desired version to which for the alleles should be updated.
#'@param expand A logical to determine whether to include only the direct HLA ID match, or all possible allele matches.
#'
#'@return A version of the input array of GL String Codes (in a data frame) updated to the desired version.
#'
#'@export
#'
#'@examples
#'\dontrun{
#'multiUpdateGL(GLSC.ex[[2]][1:5], Version = "3.53.0")
#'}
#'
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
#'@references Mack et al. HLA 05 July 2023 https://doi.org/10.1111/tan.15145
multiUpdateGL <- function(GLstringArray, Version, expand = FALSE) {
  #making sure desired version is possible
  outpos <- GLV2(Version)
  names <- colnames(GLstringArray)
  if (is.null(names)) {
    GLstringArray <-as.data.frame(GLstringArray)
  }
  names <- colnames(GLstringArray)
  if (length(GLstringArray) == 1) {
    nn <- names[1]
    for (i in 1:length(GLstringArray[[1]])) {
      GLstringArray[[nn]][i] <- updateGL(GLstringArray[[nn]][i], Version = Version, expand = expand)
    }
  }
  else {
    for (x in 2:length(GLstringArray)) {
      nn <- names[x]
      for (i in 1:length(GLstringArray[[1]])) {
        GLstringArray[[nn]][i] <- updateGL(GLstringArray[[nn]][i], Version = Version, expand = expand)
      }
    }
  }
  #View(GLstringArray)
  GLstringArray
}


################
##GLupdate
#'Update a GL String Code.
#'
#'Updates the allele names in a Genotype List (GL) String Code to the desired reference database version using the IPD-IMGT/HLA Database's Allele List History table.
#'
#'@param GLString A string of HLA allele names and operators in GL String Code format signifying their relation with one another and the pertinent HLA Allele List version.
#'@param Version A string of the desired version for the alleles to be updated to
#'@param expand A logical that determines whether to return only the direct HLA ID allele match or all possible HLA allele matches.
#'
#'@return An updated version the GL String Code input (in the form of a string) updated to the input desired version.
#'
#'@note For internal use only.
#'
#'@importFrom stringr fixed str_replace_all
#'
#'@export
#'
#'@examples
#'\dontrun{
#'GLupdate("hla#3.1.0#HLA-A*02:01+HLA-A*01:01:01:01", Version = "3.53.0")
#'GLupdate("hla#3.36.0#HLA-B*15:35", "3.52.0")
#'
#'GLupdate("hla#3.45.0#HLA-DPA1*02:01:01:05", "3.52.0")
#'GLupdate("hla#3.45.0#HLA-A*02:01:05", "3.52.0", expand = TRUE)
#'
#'GLupdate("hla#3.45.0#HLA-A*02:08", "3.52.0")
#'GLupdate("hla#3.45.0#HLA-A*02:08", "3.52.0", expand = TRUE)
#'}
#'
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
#'@references Mack et al. HLA 05 July 2023 https://doi.org/10.1111/tan.15145
#'@source https://github.com/ANHIG/IMGTHLA/blob/Latest/Allelelist_history.txt

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

  if (expand) {
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
        gBack <- str_replace_all(gBack, fixed(backBack), rTU)
      }

    }
    #changing version to reflect the update
    gBack <- gsub(origV,Version,gBack)
    gBack
  } else {
    for (i in 1:length(nameList)) {
      if(nameList[i] %in% alleleListHistory$AlleleListHistory[[inpos]]){
        #finding which row the allele is in
        ytemp <- which(alleleListHistory$AlleleListHistory[[inpos]] == nameList[i])
        #getting value of allele from desired version column
        valtemp <- alleleListHistory$AlleleListHistory[ytemp,outpos]
        #avoiding issues with NA
        if (is.na(valtemp)|is.null(valtemp)) {
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
        lastEffort <- GLupdate(tempFind, Version = origV, expand = TRUE)
        cOne <- strsplit(lastEffort, split = "#")[[1]][3]
        cOne <- strsplit(cOne, "[/]")[[1]][1]
        if (cOne != "NA") {
          message("returning match for ", cOne)
          Q <- paste0(c("hla", origV, cOne), collapse = "#")
          temp <- GLupdate(Q, Version = Version)
          temp <- strsplit(temp, split = "#")[[1]][3]
          gBack <- str_replace_all(gBack, fixed(nameList[i]), temp)
        }
        else { gBack <- str_replace_all(gBack, fixed(nameList[i]), "NA")
        }
      }
    }
    gBack <- gsub(origV,Version,gBack)
    gBack

  }

}
################
##GLV
#'Retrieve version from input GL String.
#'
#'Extracts the version data from an input GL String Code, or provides appropriate options if version input is not present in the IPD-IMGT/HLA Database.
#'
#'@param GLString A string of HLA allele names and operators in GL String Code format, signifying their relation with one another and the pertinent IPD-IMGT/HLA Database release version.
#'
#'@return Returns the version of the release as it appears in the pertinent alleleListHistory$AlleleListHistory-column-header.
#'
#'@export
#'
#'@note For internal use only.
#'
#'@examples
#'\dontrun{
#'GLV("hla#3.25.0#HLA-B15:35")
#'}
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
#'@references Mack et al. HLA 05 July 2023 https://doi.org/10.1111/tan.15145

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
#'Returns the AlleleListHistory-formatted version of a dot-delimited IPD-IMGT/HLA Database release, or provides appropriate options if version input is not present in the IPD-IMGT/HLA Database.
#'
#'@param Version  A string of the desired version for the alleles to be updated to.
#'
#'@return Returns alleleListHistory$AlleleListHistory-column-formatted version.
#'
#'@export
#'
#'@note For internal use only.
#'
#'@examples
#'\dontrun{
#'GLV2("3.34.0")
#'GLV2("3.0.0")
#'}
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
#'@references Mack et al. HLA 05 July 2023 https://doi.org/10.1111/tan.15145

GLV2 <- function(Version) {
  ve1 <-""
  ve2 <- ""
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
  }
  else {
    ve1
  }
}

################
##GLVhelper
#'Locate matches for an incomplete IPD-IMGT/HLA Database version.
#'
#' Uses the provided version information to locate and return possible matches for an incompletely defined IPD-IMGT/HLA Database release version.
#'
#'@param Version  A string of the desired version for the alleles to be updated to.
#'
#'@return Returns possible matches to a given incomplete input.
#'
#'@export
#'
#'@note For internal use only.
#'
#'@examples
#'\dontrun{
#'GLVhelper("2.25")
#'GLVhelper("3.9.0")
#'}
#'
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
#'@references Mack et al. HLA 05 July 2023 https://doi.org/10.1111/tan.15145

GLVhelper <- function(Version) {
  ve5 <-ve4 <-ve3 <-ve2 <- ""
  vte1 <- vte3 <- vector("list", 1)
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
      message("Returning all HLA Allele List versions with ", version, " as its first number, please re-enter one of these options into the function.")
      for (i in 1:length(vte1[[1]])) {
        vte1[[1]][i] <- redec(vte1[[1]][i])
      }
      print(vte1)
      stop()
    }
    else {
      message(version, " is not a valid first number for an HLA Allele List Version.")
      stop()
    }
  }
  #add 1 zero to the end, test for multiple other options
  if (length(version) == 2) {
    message("Version is too short, searching for matches by adding zero to the end.")
    version <- append(version, "0")
    for (i in 1:length(version)) {
      ve2 <- paste0(ve2,version[i])
    }
    ve2 <- paste0("X", ve2)
    ve2 <-redec(ve2)
    GLVhelper(ve2)
  }
  if (length(version) == 3) {
    tte3 <- vector("list", 1)
    vv3 <- ""
    for(i in 1:length(version)) {
      vv3 <- paste0(vv3,version[[i]])
    }
    #print(vv3)
    for(i in 1:length(version)) {
      xt <- append(version, "0", i)
      for (n in 1:length(xt)) {
        ve3 <- paste0(ve3,xt[n])
      }
      vte3[[1]][i] <- paste0("X",ve3)
      ve3 <- ""
    }
    for (i in 1: length(vte3[[1]])) {
      if(vte3[[1]][i] %in% names(alleleListHistory$AlleleListHistory)) {
        #message about what it is returning (its returning earliest version of input)
        tte3[[1]] <- append(tte3[[1]], vte3[[1]][i])
      }
    }
    for (i in 2:(length(names(alleleListHistory$AlleleListHistory)))) {
      yy <-names(alleleListHistory$AlleleListHistory)[i]
      yy <- substr(yy, start = 2, stop = 4)
      if(yy == vv3) {
        tte3[[1]] <- append(tte3[[1]], names(alleleListHistory$AlleleListHistory)[i])
      }
    }
    if (!is.null(tte3[[1]])) {
      message("No direct matches were found for an input containing the digits '", vv3, "', please re-enter one of these options into the function.")
      tte3[[1]] <- unique(tte3[[1]])
      for (i in 1:length(tte3[[1]])) {
        tte3[[1]][i] <- redec(tte3[[1]][i])
      }
      print(tte3)
      stop()
    }
    else {
      GLVhelper(version[1])
    }
  }

  #remove last element from list, and see if that works
  if (length(version) > 3) {
    xnames <- names(alleleListHistory$AlleleListHistory)[-1]
    for (i in 1:length(xnames)) {
      xnames[i] <- redec(xnames[i])
    }
    message("Please re-enter a version from this list of names. \n" )
    print(xnames)
    stop()
  }

}

################
##redec
#'Reintroduce version decimals.
#'
#'Correctly put decimals back into an AlleleListHistory column-formatted version name.
#'
#'@param Cname Version name in AlleleListHistory column format.
#'
#'@return Version name with decimals in appropriate places according to the IPD-IMGT/HLA Database.
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
#'Validate a GL String Code.
#'
#'A lightweight validator to protect from incorrectly formatted GL String Codes.
#'
#'@param GLString A string formatted to GL String Code specifications.
#'
#'@return A well-formatted GL String Code, or an error message.
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

GLvalidate <- function(GLString) {
  GLcopy <- GLString
  GLret <- ""
  if (length(strsplit(GLString, "#")[[1]]) > 3) {
    message("The input GL String code is formatted incorrectly, there are too many '#' delimited sections.")
    stop()
  }
  if (length(strsplit(GLString, "#")[[1]]) < 3) {
    message("The input GL String code is formatted incorrectly, there are too few '#' delimited sections.")
    stop()
  }
  if (length(strsplit(GLcopy, "#")[[1]]) == 3) {
    whole <- strsplit(GLcopy, '#')[[1]]
    first <- strsplit(GLcopy, '#')[[1]][1]
    se <- GLV(GLString)
    third <- strsplit(GLcopy, "#")[[1]][3]
    if (first != "hla") {
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
