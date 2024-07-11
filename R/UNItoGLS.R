#UNItoGLS (uniformat to GL string) v1.0.0 7FEB2024

#library(stringr)

################
##UNItoGLS
#'Translate UNIFORMAT to GL String
#'
#'A wrapper function for UNtoGL() which translates strings from UNIFORMAT format to GL String format.
#'
#'@param uniformat A character string of HLA allele names and operators in the UNIFORMAT format signifying their relation with one another.
#'@param prefix A character string of the desired gene-system prefix (default is "HLA-").
#'@param pre A logical that indicates whether returned allele names should contain 'prefix' (TRUE), or if 'prefix' should be excluded from returned names (FALSE).
#'
#'@return  An version of 'uniformat' converted to GL String format, or FALSE if 'uniformat' is invalid.
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

  if(validateUniformat(uniformat)){
    
          UNtoGL(uniformat = uniformat, prefix = prefix, pre = pre) } else { 
              FALSE}
}

################
##multiUNItoGLS
#'Translate Multiple UNIFORMAT Strings to GL Strings
#'
#'Translate a data frame or vector of UNIFORMAT strings to GL Strings.
#'
#'@param UniformatArray A data frame or vector of UNIFORMAT formatted strings. If 'UniformatArray' is a data frame with more than one column, the first column should contain only identifiers. If 'UniformatgArray' is a vector, it should contain only UNIFORMAT strings.
#'@param prefix A character string of the desired locus prefix (default is "HLA-").
#'@param pre A logical. If 'pre' is TRUE, all allele names will be prefixed with 'prefix'. If 'pre' is FALSE, no allele names will be prefixed.
#'
#'@return A version of 'UniformatArray' in which the UNIFORMAT data have been converted to GL String format. If a 'UniformatArray' was a data frame, a data frame is returned. If 'UniformatArray' was a vector, a vector is returned. 
#'
#'@note This function does not return the GL String "?" operator, as the "?" operator has no cognate in UNIFORMAT.
#'
#'@export
#'
#'@examples
#'multiUNItoGLS(unname(as.vector(UNIFORMAT.example[2]))[[1]])
#'multiUNItoGLS(UNIFORMAT.example[1:2])
#'
#'@references Nunes Tissue Antigens 2007;69(s1):203-205 https://doi.org/10.1111/j.1399-0039.2006.00808.x
#'@references Mack et al. HLA 2023;102(2):206-212 https://doi.org/10.1111/tan.15126
multiUNItoGLS <- function(UniformatArray, prefix = "HLA-", pre = TRUE) {
 
  vec <- FALSE
  if (is.vector(UniformatArray)) {
      vec <- TRUE
      UniformatArray <-as.data.frame(UniformatArray)
    }
  
  names <- colnames(UniformatArray)

    for(x in ifelse(length(UniformatArray) == 1,1,2):length(UniformatArray)) {
 
      nn <- names[x]
      for (i in 1:length(UniformatArray[[x]])) {
 
                if (suppressMessages(validateUniformat(UniformatArray[[nn]][i]))) {
          UniformatArray[[nn]][i] <- UNtoGL(UniformatArray[[nn]][i], prefix = prefix, pre = pre)
        } else {

          message( paste( "The UNIFORMAT string in row ",i,"contains operators or characters that are not permitted in UNIFORMAT."))
        }
      }
    }

  if(vec) { UniformatArray <- unname(unlist(UniformatArray)) }
  
  UniformatArray
}

################
##UNtoGL
#'Translate UNIFORMAT Strings to GL Strings
#'
#'A function that translates strings from UNIFORMAT format to GL String format.
#'
#'@param uniformat A character string of HLA allele names and operators in the UNIFORMAT format signifying their relation with one another.
#'@param prefix A character string of the desired prefix (default is "HLA-").
#'@param pre A logical that indicates whether user would like all allele names to contain the prefix of their choice (TRUE), or if the prefix should not be appended to allele names (FALSE).
#'
#'@return An altered version of the input character string, converted to GL String format.
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
 
  GLret <- uniformat
  
  preTemp <- prefix

  GLret <- gsub(prefix, "", GLret)

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
##validateUniformat
#'Validate a UNIFORMAT String
#'
#'Evaluates a UNIFORMAT string to identify unsupported characters.
#'
#'@param uniformat A character string of allele names and operators in the UNIFORMAT format.
#'
#'@return A logical. TRUE is returned when all characters in 'uniformat' are permitted. FALSE is returned when forbidden characters are present.
#'
#'@export
#'
#'@examples
#'validateUniformat("A*02:01,A*03:01|A*02:01,A*03:02|A*02:02,A*03:01|A*02:02,A*03:02")
#'
#'@references Nunes Tissue Antigens 2007;69(s1):203-205 https://doi.org/10.1111/j.1399-0039.2006.00808.x
validateUniformat <- function(uniformat) {
 
  badOperate <- all(strsplit(uniformat,split="")[[1]] %in% c(" ","|",",","\t","~",0:9,LETTERS,"b", "l", "a", "n", "k", "g","*","-",":"))
  
  #populating return sentence with what is wrong
  if (!badOperate) {  message(paste(uniformat, "contains operators or characters not permitted in UNIFORMAT.",sep=" "))  }

  badOperate

}
