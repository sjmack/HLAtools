#GLStoUNI (GL string to uniformat) v1.0.0 12APR2024

#library(stringr)

################
##GLStoUNI
#'Translate GL String to UNIFORMAT
#'
#'A wrapper function for GLtoUN, which translates strings from GL String format to UNIFORMAT format.
#'
#'@param GLstring A character string of HLA allele names and operators in the GL String format signifying their relation with one another.
#'@param prefix A character string of the pertinent gene-system prefix (default is "HLA-").
#'@param pre A logical that indicates whether returned allele names should contain 'prefix' (TRUE), or if 'prefix' should be excluded from returned names (FALSE).
#'
#'@return A version of 'GLstring' converted to UNIFORMAT format, or FALSE if 'GLstring' is invalid.
#'
#'@note This function does not function with the "?" operator, as the "?" operator has no cognate in UNIFORMAT.
#'
#'@export
#'
#'@examples
#'
#'GLStoUNI("HLA-A*02:01/HLA-A*02:02+HLA-A*03:01/HLA-A*03:02")
#'GLStoUNI("A+B/C~D^G|E^W+X/Y^Z+J")
#'
#'@references Nunes Tissue Antigens 2007;69(s1):203-205 https://doi.org/10.1111/j.1399-0039.2006.00808.x
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
GLStoUNI <- function(GLstring, prefix = "HLA-", pre = FALSE) {
  #using validator to better protect code
  if(validateGLstring(GLstring,"1.0")) {
          #if all is well, call GLtoUNI
        GLtoUN(GLstring = GLstring, prefix = prefix, pre = pre) } else {
          
          return(FALSE)
        }
}

################
##multiGLStoUNI
#'Translate Multiple GL Strings to UNIFORMAT
#'
#'Translate a data frame or a vector of GL Strings to UNIFORMAT strings.
#'
#'@param GLstringArray A data frame or a vector containing GL String formatted data. If 'GLstringArray' is a data frame with more than one column, the first column should contain only identifiers. If 'GLstringArray' is a vector, it should contain only GL Strings.
#'@param prefix A string of the desired locus prefix (default is "HLA-").
#'@param pre A logical. If 'pre' is TRUE, all allele names will be prefixed with 'prefix'. If 'pre' is FALSE, no allele names will be prefixed.
#'
#'@return A version of 'GLStringArray' in which the GL String data has been converted to UNIFORMAT format. If a 'GLstringArray' was a data frame, a data frame is returned. If 'GLstringArray' was a vector, a vector is returned. 
#'
#'@note GL Strings that include the "?" operator will not be translated, as the "?" operator has no cognate in UNIFORMAT.
#'
#'@export
#'
#'@examples
#'multiGLStoUNI(GLstring.ex[[2]][1:5],version) ## converting a vector of GL Strings to UNIFORMAT
#'multiGLStoUNI(GLstring.ex[1:5,1:2])  ## converting a data.frame of GL Strings to UNIFORMAT
#'
#'@references Nunes Tissue Antigens 2007;69(s1):203-205 https://doi.org/10.1111/j.1399-0039.2006.00808.x
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
multiGLStoUNI <- function(GLstringArray, prefix = "HLA-", pre = FALSE) {

    vec <- FALSE
    
      if(is.vector(GLstringArray)) {
        vec <- TRUE
        GLstringArray <-as.data.frame(GLstringArray)
      }
    
    names <- colnames(GLstringArray)

    for(x in ifelse(length(GLstringArray) == 1,1,2):length(GLstringArray)) {
      nn <- names[x]
      for(i in 1:length(GLstringArray[[1]])) {
        if(suppressMessages(validateGLstring(GLstringArray[[nn]][i],"1.0"))) {
          GLstringArray[[nn]][i] <- GLtoUN(GLstringArray[[nn]][i], prefix = prefix, pre = pre)
            } else {
                  message(paste("The GL String in row",i,"contains operators or characters that are not permitted in GL String version 1.0.",sep=" "))
                    }
               }
          }
  
  if(vec) { GLstringArray <- unname(unlist(GLstringArray)) }
    
  GLstringArray
}

################
##GLtoUN
#'Translate GL Strings to UNIFORMAT Strings
#'
#'A function that translates strings from GL String format to UNIFORMAT format.
#'
#'@param GLstring A character string of HLA allele names and operators in the GL String format signifying their relation with one another.
#'@param prefix A character string of the desired gene-system prefix (default is "HLA-").
#'@param pre A logical that indicates whether returned allele names should contain 'prefix' (TRUE), or if 'prefix' should be excluded from returned names (FALSE).
#'
#'@return A version of the input string converted to UNIFORMAT format.
#'
#'@note For internal use only.
#'@note This function does not function with the "?" operator, as the "?" operator has no cognate in UNIFORMAT.
#'
#'@importFrom stringr fixed str_replace_all
#'
#'@export
#'
#'@examples
#'GLtoUN("HLA-A*02:01/HLA-A*02:02+HLA-A*03:01/HLA-A*03:02")
#'GLtoUN("A+B/C~D^G|E^W+X/Y^Z+J")
#'
#'@references Nunes Tissue Antigens 2007;69(s1):203-205 https://doi.org/10.1111/j.1399-0039.2006.00808.x
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
GLtoUN <- function(GLstring, prefix = "HLA-", pre = FALSE) {

  UNIret <- GLstring
  
  UNIret <- gsub(prefix, "", UNIret)

  #splitting into large sections, which will be split further
  firstSplit <- strsplit(UNIret, "\\^|[~]|[|]")[[1]]
  
  for (i in 1: length(firstSplit)) {

    #note: SS = second split, TS = third split
    secondSplit <- strsplit(firstSplit[i], "[+]")[[1]]
    thirdSplit <- strsplit(firstSplit[i], "[/]|[+]")[[1]]
    nameSS <-unique(strsplit(firstSplit[i], "[/]|[+]")[[1]])
    lengthSS <- length(secondSplit)
    lengthTS <- length(thirdSplit)
    
    #only invoked for specific scenarios involving cartesian products
    if(lengthSS>=2 && lengthTS >=3) {
      firstSSplit <- strsplit(secondSplit[1], "[/]")[[1]]
      secondSSplit <- strsplit(secondSplit[2], "[/]")[[1]]
      list1 <- vector("list", 1)[[1]]
      tempS <- ""
      #for (n in 1:lengthSS) {
      for (n in 1:length(firstSSplit)) {
        #for (m in 1:lengthTS) {
        for (m in 1:length(secondSSplit)) {
          tempSS <- ""
          tempSS <- paste0(paste0(firstSSplit[n], ","), secondSSplit[m])
          list1 <- append(list1, fixed(tempSS), (length(list1)))
          list1 <- unique(list1)
        }
      }
      fixOne <- paste0(list1, collapse = "|")
      UNIret <- str_replace_all(UNIret, fixed(firstSplit[i]), fixOne)
    }

  }

  #substituting the one to one translations
  UNIret <- gsub("[+]",",",UNIret)

  UNIret <- gsub("\\^"," ",UNIret)

  if (pre) {
    baseSplit <- unique(strsplit(UNIret, "[,]|[ ]|[\t]|[|]|[~]|\\^")[[1]])
    finalSplit <- unique(strsplit(UNIret, "[,]|[ ]|[\t]|[|]|[~]|\\^")[[1]])
    for (i in 1:length(finalSplit)) {
      finalSplit[i] <- paste0(prefix, finalSplit[i])
      UNIret <- str_replace_all(UNIret, fixed(baseSplit[i]), fixed(finalSplit[i]))
    }
  }

  UNIret
}

################
##validateGLstring
#'Validate a GL String
#'
#'Evaluates a GL String to identify unsupported characters for a specified GL String version.
#'
#'@param GLstring A character string of allele names and operators in the GL String format.
#'
#'@param version A character string identifying the GL String version for evaluation. Options are "1.0" and "1.1".
#'
#'@return A logical. TRUE is returned when all characters in a 'version' GL String are permitted. FALSE is returned when forbidden characters are present, or when an incorrect 'version' is specified.
#'
#'@export
#'
#'@examples
#'validateGLstring("HLA-A*02:01/HLA-A*02:02+HLA-A*03:01/HLA-A*03:02", version = "1.0")
#'
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
validateGLstring <- function(GLstring, version) {
  #if there is anything not supposed to be there, badOperate will be false
  
  if(!version %in% c("1.0","1.1")) {return(FALSE)}
  
    if(version == "1.0") {
      badOperate <- all(strsplit(GLstring,split="")[[1]] %in% c("/","|","+","^","~",0:9,LETTERS, "b", "l", "a", "n", "k", "g","*","-",":"))
        } 
  
    if(version == "1.1") {
      badOperate <- all(strsplit(GLstring,split="")[[1]] %in% c("/","|","+","^","~","?",0:9,LETTERS, "b", "l", "a", "n", "k", "g","*","-",":"))
    }
    
     if (!badOperate) { message(paste(GLstring," contains characters or operators not permitted in GL String version ",version,".",sep="")) }

  badOperate
}
