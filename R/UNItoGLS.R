#UNItoGLS (uniformat to GL string) v1.0.0 7FEB2024

#library(stringr)

################
##UNItoGLS
#'Translate strings from UNIFORMAT format to GL String format.
#'
#'A wrapper function for UNtoGL(), which translates strings from UNIFORMAT format to GL String format.
#'
#'@param uniformat A string of HLA allele names and operators in the UNIFORMAT format signifying their relation with one another.
#'@param prefix A string of the desired prefix (autoset to "HLA-"). This string should contain a "-" at the end.
#'@param pre A logical that indicates whether user would like all allele names to contain the prefix of their choice (TRUE), or if the prefix should not be appended to allele names (FALSE).
#'
#'@return An altered version of the input string, converted to GL String format.
#'
#'@note This function does not return the "?" operator, as the "?" operator has no cognate in UNIFORMAT.
#'
#'@export
#'
#'@examples
#'UNItoGLS("A*02:01,A*03:01|A*02:01,A*03:02|A*02:02,A*03:01|A*02:02,A*03:02")
#'UNItoGLS("A,B|A,C~D G|E W,X|W,Y Z,J J,Z")
#'
#'@references Nunes Tissue Antigens 2007;69(s1):203-205 https://doi.org/10.1111/j.1399-0039.2006.00808.x
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
UNItoGLS <- function(uniformat, prefix = "HLA-", pre = TRUE) {
  #using validator to better protect code
  uniformat <- uniformatValidate(uniformat)
  #if all is well, call GLupdate
  UNtoGL(uniformat = uniformat, prefix = prefix, pre = pre)
}

################
##multiUNItoGLS
#'Translate columns of UNIFORMAT strings to GL Strings.
#'
#'A function that translates columns of arrays containing UNIFORMAT formatted strings to GL String format.
#'
#'@param UniformatArray A vector with UNIFORMAT formatted strings. if the vector has more than one column, the first column should only contain labels/identification data. The first column of a multi-column vector should NOT contain UNIFORMAT strings. if the input vector only contains one column, the column should contain UNIFORMAT strings.
#'@param prefix A string of the desired prefix (autoset to "HLA-"), this string should contain a "-" at the end.
#'@param pre A logical that indicates whether user would like all allele names to contain the prefix of their choice (TRUE), or if the prefix should not be appended to allele names (FALSE).
#'
#'@return An altered version of the input string, converted to GL String format.
#'
#'@note This function does not return the "?" operator, as the "?" operator has no cognate in UNIFORMAT.
#'
#'@export
#'
#'@examples
#'multiUNItoGLS(UNIFORMAT.example[1:2])
#'multiUNItoGLS(UNIFORMAT.example[2])
#'
#'@references Nunes Tissue Antigens 2007;69(s1):203-205 https://doi.org/10.1111/j.1399-0039.2006.00808.x
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
multiUNItoGLS <- function(UniformatArray, prefix = "HLA-", pre = TRUE) {
  names <- colnames(UniformatArray)
  if (is.null(names)) {
    UniformatArray <-as.data.frame(UniformatArray)
  }
  names <- colnames(UniformatArray)
  if (length(UniformatArray) == 1) {
    #nn <- colnames(GLstringArray, do.NULL = FALSE)
    nn <- names[1]
    for (i in 1:length(UniformatArray[[1]])) {
      #if there is anything not supposed to be there, badOperate will be false
      badOperate <- all(strsplit(UniformatArray[[nn]][i],split="")[[1]] %in% c(" ","|",",","\t","~",0:9,LETTERS,"b", "l", "a", "n", "k", "g","*","-",":"))
      if (badOperate) {
        UniformatArray[[nn]][i] <- UNtoGL(UniformatArray[[nn]][i], prefix = prefix, pre = pre)
      }
      #populating return sentence with what is wrong
      if (badOperate == FALSE) {
        #simple return message
        message( paste( "In row", i, ",the input uniformat string contains operators or characters not present in uniformat."))
      }
    }
  }
  else {
    for (x in 2:length(UniformatArray)) {
      nn <- names[x]
      for (i in 1:length(UniformatArray[[1]])) {

        #if there is anything not supposed to be there, badOperate will be false
        badOperate <- all(strsplit(UniformatArray[[nn]][i],split="")[[1]] %in% c(" ","|",",","\t","~",0:9,LETTERS,"b", "l", "a", "n", "k", "g","*","-",":"))
        if (badOperate) {
          UniformatArray[[nn]][i] <- UNtoGL(UniformatArray[[nn]][i], prefix = prefix, pre = pre)
        }
        #populating return sentence with what is wrong
        if (badOperate == FALSE) {
          #simple return message
          message( paste( " In row ", i, ", the input uniformat string contains operators or characters not present in uniformat."))
        }
      }
    }
  }
  #View(UniformatArray)
  UniformatArray
}

################
##UNtoGL
#'Translate UNIFORMAT strings to GL Strings.
#'
#'A function that translates strings from UNIFORMAT format to GL String format.
#'
#'@param uniformat A string of HLA allele names and operators in the UNIFORMAT format signifying their relation with one another.
#'@param prefix A string of the desired prefix (autoset to "HLA-"), this string should contain a "-" at the end.
#'@param pre A logical that indicates whether user would like all allele names to contain the prefix of their choice (TRUE), or if the prefix should not be appended to allele names (FALSE).
#'
#'@return An altered version of the input string, converted to GL String format.
#'
#'@note For internal use only.
#'@note This function does not return the "?" operator, as the "?" operator has no cognate in UNIFORMAT.
#'
#'@importFrom stringr fixed str_replace_all
#'
#'@export
#'
#'@examples
#'UNtoGL("A*02:01,A*03:01|A*02:01,A*03:02|A*02:02,A*03:01|A*02:02,A*03:02")
#'UNtoGL("A,B|A,C|D,B|D,C|E,B|E,C", pre = FALSE)
#'
#'@references Nunes Tissue Antigens 2007;69(s1):203-205 https://doi.org/10.1111/j.1399-0039.2006.00808.x
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
UNtoGL <- function(uniformat, prefix = "HLA-", pre = TRUE) {
  #print(uniformat)
  #this string will be altered and returned
  GLret <- uniformat
  preTemp <- prefix
  preTemp <- gsub("-", "", preTemp)
  #taking out all prefixes just in case they are present. this will be returned depending on prefix parameter
  GLret <- gsub(prefix, "", GLret)
  GLret <- gsub("-", "", GLret)
  GLret <- gsub(preTemp, "", GLret)

  #splitting across some operators, allowing more difficult translations to be carried out on segments
  firstSplit <- strsplit(GLret, "[ ]|[~]|[\t]")[[1]]
  for (i in 1: length(firstSplit)) {
    #divide alleles amongst these lists
    list1 <- vector("list", 1)
    list2 <- vector("list", 1)
    #more refined splits and other information, to be used when creating "/" denoted lists of ambiguous alleles
    #note: SS = second split, TS = third split
    secondSplit <- strsplit(firstSplit[i], "[|]")[[1]]
    nameSS <-unique(strsplit(firstSplit[i], "[|]|[,]")[[1]])
    #print(nameSS)
    lengthSS <- length(secondSplit)
    lengthTS <- length(strsplit(firstSplit[i], "[|]|[,]")[[1]])
    #only invoked for specific scenarios involving cartesian products
    if(lengthSS>=2 && lengthTS>2) {
      fixOne <- fixTwo <- ""
      #populating lists with first pair
      initialPopulate <- strsplit(secondSplit[1], "[,]")[[1]]
      list1 <- append(list1[[1]], initialPopulate[1])
      list2 <- append(list2[[1]], initialPopulate[2])
      #removing initial populaters names from nameSS
      nameSS <- nameSS[nameSS != initialPopulate[1]]
      nameSS <- nameSS[nameSS != initialPopulate[2]]
      #secondSplit <- secondSplit[secondSplit != secondSplit[1]]
      #print(secondSplit)
      #print(length(nameSS))
      for (n in 1:length(nameSS)) {
        #sorting to either list, if two items are in same section (e.g. |a,c|) they will be placed on different lists
        FirstList <- TRUE
        #specific scenario involving the "blank" allele
        isBlank <- FALSE
        for (w in 1:lengthSS){
          SmallSub <- strsplit(secondSplit[w], "[,]")[[1]]
          #if there are no names, this means there are only repeats, so set isBlank to true and break
          if (length(nameSS) == 0) {
            isBlank <- TRUE
            break()
          }
          if (nameSS[n] %in% SmallSub && list1[1] %in% SmallSub){
            FirstList = FALSE
            break()
          }
          if (nameSS[n] %in% SmallSub && list2[1] %in% SmallSub){
            FirstList = TRUE
            break()
          }
        }
        #handling if there is a blank allele possibility
        if (isBlank) {
          if ("blank" %in% list1) {
            list1 <- append(list1, list2, (length(list1)))
            list1 <- unique(list1)
          }
          if ("blank" %in% list2) {
            list2 <- append(list2, list1, (length(list1)))
            list2 <- unique(list2)
          }
          if (length(secondSplit) == 2) {
            list1 <- list1[list1 %in% NA == FALSE]
            list2 <- list2[list2 %in% NA == FALSE]
          }
        }
        if (FirstList) {
          #list1[[length(list1)+1]] = nameSS[n]
          list1 <- append(list1, nameSS[n], (length(list1)))
        }
        #if (!FirstList) {
        else {
          #list2[[length(list2)+1]] = nameSS[n]
          list2 <- append(list2, nameSS[n], (length(list2)))
        }
      }
      #creating single string with "/" denoted lists of ambiguous alleles
      fixOne <- paste0(list1, collapse="/")
      fixOne <- paste0(fixOne, "+")
      fixTwo <- paste0(list2, collapse="/")
      fixOne <- paste0(fixOne, fixTwo)
      #print(fixOne)
      GLret <- str_replace_all(GLret, fixed(firstSplit[i]), fixOne)
    }

  }
  #substituting the one to one translations
  GLret <- gsub(",","+",GLret)
  GLret <- gsub(" ","^",GLret)
  GLret <- gsub("\t","^",GLret)

  if (pre) {
    baseSplit <- unique(strsplit(GLret, "[.]|[/]|[+]|[|]|[~]|\\^")[[1]])
    finalSplit <- unique(strsplit(GLret, "[.]|[/]|[+]|[|]|[~]|\\^")[[1]])
    #print(finalSplit)
    for (i in 1:length(finalSplit)) {
      #print(finalSplit[i])
      finalSplit[i] <- paste0(prefix, fixed(finalSplit[i]))
      GLret <- str_replace_all(GLret, fixed(baseSplit[i]), fixed(finalSplit[i]))
    }
  }

  GLret
}


################
##uniformatValidate
#'Validates a UNIFORMAT string.
#'
#'A function that validates that a UNIFORMAT formatted string, ensuring that it does not contain any characters not supported in UNIFORMAT.
#'
#'@param uniformat A string of HLA allele names and operators in the UNIFORMAT format signifying their relation with one another.
#'
#'@return If the input string is valid, the original string is returned. If the input string is invalid, an error message is generated.
#'
#'@export
#'
#'@examples
#'uniformatValidate("A*02:01,A*03:01|A*02:01,A*03:02|A*02:02,A*03:01|A*02:02,A*03:02")
#'
#'@references Nunes Tissue Antigens 2007;69(s1):203-205 https://doi.org/10.1111/j.1399-0039.2006.00808.x
uniformatValidate <- function(uniformat) {
  #if there are any of these, badOperate will be true
  badOperate <- all(strsplit(uniformat,split="")[[1]] %in% c(" ","|",",","\t","~",0:9,LETTERS,"b", "l", "a", "n", "k", "g","*","-",":"))
  #populating return sentence with what is wrong
  if (badOperate == FALSE) {
    message("The input uniformat string contains operators or characters not present in UNIFORMAT.")
    stop()
  }

  uniformat

}
