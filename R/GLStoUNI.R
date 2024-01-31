#GLStoUNI (GL string to uniformat) v0.2.1 12AUG23

#library(stringr)

################
##GLStoUNI
#'Translates strings from GL String format to UNIFORMAT format.
#'
#'A wrapper function for GLtoUN, a function that translates strings from GL String format to UNIFORMAT format.
#'
#'@param GLstring A string of HLA allele names and operators in the GL String format, identifying their relation with one another.
#'@param prefix A string of the desired prefix (default is "HLA-"), this string should contain a "-" at the end.
#'@param pre A logical that indicates whether user would like all allele names to contain the prefix of their choice (TRUE), or if the prefix should not be appended to allele names (FALSE).
#'
#'@return An altered version of the input string, converted to UNIFORMAT format.
#'
#'@note This function does not function with the "?" operator, as the "?" operator has no cognate in UNIFORMAT.
#'
#'@export
#'
#'@examples
#'GLStoUNI("HLA-A*02:01/HLA-A*02:02+HLA-A*03:01/HLA-A*03:02")
#'GLStoUNI("A+B/C~D^G|E^W+X/Y^Z+J")
#'
#'@references Nunes Tissue Antigens 2007;69(s1):203-205 https://doi.org/10.1111/j.1399-0039.2006.00808.x
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
GLStoUNI <- function(GLstring, prefix = "HLA-", pre = FALSE) {
  #using validator to better protect code
  GLstring <- glstringValidate(GLstring)
  #if all is well, call GLtoUNI
  GLtoUN(GLstring = GLstring, prefix = prefix, pre = pre)
}

################
##multiGLStoUNI
#'Translates columns of GL Strings to UNIFORMAT strings.
#'
#'A function that translates columns of arrays containing strings in GL String format to UNIFORMAT format.
#'
#'@param GLstringArray A vector with GL String format. If the vector has more than one column, the first column should only contain labels/identifiers; the first column of a multi-column vector should NOT contain GL Strings. If the input vector only contains one column, the column should contain GL Strings.
#'@param prefix A string of the desired prefix (autoset to "HLA-"), this string should contain a "-" at the end.
#'@param pre A logical that indicates whether user would like all allele names to contain the prefix of their choice (TRUE), or if the prefix should not be appended to allele names (FALSE).
#'
#'@return An altered version of the input array, converted to UNIFORMAT format.
#'
#'@note This function does not function with the "?" operator, as the "?" operator has no cognate in UNIFORMAT.
#'
#'@export
#'
#'@examples
#'multiGLStoUNI(GLstring.ex[[2]][1:5])
#'
#'@references Nunes Tissue Antigens 2007;69(s1):203-205 https://doi.org/10.1111/j.1399-0039.2006.00808.x
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
multiGLStoUNI <- function(GLstringArray, prefix = "HLA-", pre = FALSE) {
  names <- colnames(GLstringArray)
  if (is.null(names)) {
    GLstringArray <-as.data.frame(GLstringArray)
  }
  names <- colnames(GLstringArray)
  if (length(GLstringArray) == 1) {
    #nn <- colnames(GLstringArray, do.NULL = FALSE)
    nn <- names[1]
    for (i in 1:length(GLstringArray[[1]])) {
      #if there is anything not supposed to be there, badOperate will be false
      badOperate <- all(strsplit(GLstringArray[[nn]][i],split="")[[1]] %in% c("/","|","+","^","~",0:9,LETTERS, "b", "l", "a", "n", "k", "g","*","-",":"))
      if (badOperate) {
        GLstringArray[[nn]][i] <- GLtoUN(GLstringArray[[nn]][i], prefix = prefix, pre = pre)
      }
      #populating return sentence with what is wrong
      if (badOperate == FALSE) {
        #simple return message
        message( paste( "In row", i, ",the input GL String contains operators or characters not present in GL String version 1.0."))
      }
    }
  }
  else {
    for (x in 2:length(GLstringArray)) {
      nn <- names[x]
      for (i in 1:length(GLstringArray[[1]])) {

        #if there is anything not supposed to be there, badOperate will be false
        badOperate <- all(strsplit(GLstringArray[[nn]][i],split="")[[1]] %in% c("/","|","+","^","~",0:9,LETTERS, "b", "l", "a", "n", "k", "g","*","-",":"))
        if (badOperate) {
          GLstringArray[[nn]][i] <- GLtoUN(GLstringArray[[nn]][i], prefix = prefix, pre = pre)
        }
        #populating return sentence with what is wrong
        if (badOperate == FALSE) {
          #simple return message
          message( paste( " In row ", i, ", the input GL String contains operators or characters not present in GL String version 1.0."))
        }
      }
    }
  }
  #View(GLstringArray)
  GLstringArray
}

################
##GLtoUN
#'Translates GL Strings to UNIFORMAT strings.
#'
#'A function that translates strings from GL String format to UNIFORMAT format.
#'
#'@param GLstring A string of HLA allele names and operators in the GL String format signifying their relation with one another.
#'@param prefix A string of the desired prefix (default is "HLA-"), this string should contain a "-" at the end.
#'@param pre A logical that indicates whether user would like all allele names to contain the prefix of their choice (TRUE), or if the prefix should not be appended to allele names (FALSE).
#'
#'@return An altered version of the input string, altered to be in UNIFORMAT format.
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
  #print(GLstring)
  #this string will be altered and returned
  UNIret <- GLstring
  preTemp <- prefix
  preTemp <- gsub("-", "", preTemp)
  #removing "HLA-" prefix, will be put back if prefix = TRUE
  UNIret <- gsub(prefix, "", UNIret)
  UNIret <- gsub("-","",UNIret)
  UNIret <- gsub(preTemp,"",UNIret)

  #splitting into large sections, which will be split further
  #firstSplit <- strsplit(GLstring, "\\^|[~]|[|]")[[1]]
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
  #UNIret <- gsub("\\^","\t",UNIret)
  UNIret <- gsub("\\^"," ",UNIret)

  if (pre) {
    baseSplit <- unique(strsplit(UNIret, "[,]|[ ]|[\t]|[|]|[~]|\\^")[[1]])
    finalSplit <- unique(strsplit(UNIret, "[,]|[ ]|[\t]|[|]|[~]|\\^")[[1]])
    for (i in 1:length(finalSplit)) {
      finalSplit[i] <- paste0("HLA-", finalSplit[i])
      #UNIret <- gsub(baseSplit[i], finalSplit[i], UNIret)
      UNIret <- str_replace_all(UNIret, fixed(baseSplit[i]), fixed(finalSplit[i]))
    }
  }

  UNIret
}

################
##glstringValidate
#'Validate a GL String.
#'
#'A function that validates that a GL String formatted string does not contain any letters or operators not supported in the GL String.
#'
#'@param GLstring A string of HLA allele names and operators in the GL String format signifying their relation with one another.
#'
#'@return If everything is in order, it will return the original string, if something is wrong, it will return a message and stop.
#'
#'@export
#'
#'@examples
#'glstringValidate("HLA-A*02:01/HLA-A*02:02+HLA-A*03:01/HLA-A*03:02")
#'
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
glstringValidate <- function(GLstring) {
  #if there is anything not supposed to be there, badOperate will be false
  badOperate <- all(strsplit(GLstring,split="")[[1]] %in% c("/","|","+","^","~",0:9,LETTERS, "b", "l", "a", "n", "k", "g","*","-",":"))

  #populating return sentence with what is wrong
  if (badOperate == FALSE) {
    #simple return message
    message("The input GL String contains operators or characters not present in GL String version 1.0.")
    stop()
  }

  GLstring

}
