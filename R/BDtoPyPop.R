## BDtoPyPop v1.0.1 10 Jul 2024

################
##FormatHead
#'Format PyPop Data Frame Headers
#'  
#'Format the header of a PyPop-formatted data frame.
#'
#'@param colHead A vector of column names. PyPop format requires that paired locus/gene names should end in '_1' and '_2', respectively.
#'
#'@return A vector in which the locus names are suffixed with '_1' and '_2'.
#'
#'@note This function assumes that the first two elements the 'colHead' vector are not locus/gene names.
#'
#'@references Lancaster et al. Front Immunol. 2024 Apr 2;15:1378512. https://pubmed.ncbi.nlm.nih.gov/38629078/
#'@references Pappas et al. Hum Immunol. 2016 Mar 77(3):283-287. https://pubmed.ncbi.nlm.nih.gov/26708359/
#'
#'@export
#'
#'@examples
#'formatHead(colHead = colnames(sHLAdata))
#'
formatHead <- function(colHead) {
  
  headLength <- length(colHead) 
  
  colHead[seq(3,headLength,by=2)] <- paste(colHead[seq(3,headLength,by=2)],"_1",sep="")
  colHead[seq(4,headLength,by=2)] <- paste(colHead[seq(4,headLength,by=2)],"_2",sep="")
  
  colHead
}

################
##convertAny
#' Convert Values Across an Entire Data Frame or Vector
#'  
#' Converts all instances of a value in a data frame or vector to a specified value.
#'
#' @param dataset A data frame or vector of values.
#' @param change.from The value present in the data frame or vector that should be changed. The default value is NA.
#' @param change.to The value to which 'change.from' should be changed to. The default value is "****".
#'
#' @return A data frame or vector in which 'change.from' has been converted to 'change.to'.
#'
#' @import dplyr
#'
#' @export
#' 
#' @examples
#' dataset <- convertAny(sHLAdata)
#'
convertAny <- function(dataset,change.from=NA,change.to="****") {
  
  if(is.na(change.from)) {
    change.from <- is.na(dataset)
  } else {
    change.from <- (dataset == change.from)
  }
  dataset %>% replace(change.from,change.to)
}

################
##pypopHeaders
#'Convert BIGDAWG File Headers to PyPop Format 
#'
#'Convert the header of a BIGDAWG-formatted data frame into a Python for Population Genomics (PyPop) formatted header.
#'
#'@param colHead A vector of column names. BIGDAWG format requires that the first two header values are the sample identifier character string and the subject status (0 or 1). All other header values should be paired locus/gene names.
#'
#'@return A PyPop-formatted column header vector.
#'
#'@export
#'
#'@references Lancaster et al. Front Immunol. 2024 Apr 2;15:1378512. https://pubmed.ncbi.nlm.nih.gov/38629078/
#'@references Pappas et al. Hum Immunol. 2016 Mar 77(3):283-287. https://pubmed.ncbi.nlm.nih.gov/26708359/
#'
#'@note This function returns unique locus names for each column. For example, "A_1" and "A_2" headers will be returned for two "A" locus columns. PyPop configuration (.ini) files allow this format to be customized. 
#'@seealso [PyPop Configuration File](http://pypop.org/docs/guide-chapter-usage.html#a-minimal-configuration-file)
#'
#'@examples
#'pyHead <- pypopHeaders(colnames(sHLAdata))
#'
pypopHeaders <- function(colHead) {
  
  if(length(grep("[_][1]|[_][2]",colHead)) == 0) { ## no underscores
    
          if(length(grep("[.][1]|[.[2]",colHead)) == 0) { # no periods
           
            colHead <- formatHead(colHead)
      
            } else { ## there are *some* periods
              
              targets <- grep("[.][1]|[.[2]",colHead)
              if(length(targets) == length(colHead)-2) { # everything has a period
                
                colHead <- gsub(".","_",colHead,fixed=TRUE)
                
                    } else { ## guessing here all of the periods are ".1" on the second locus name
                
                      if(all(grep(".1",colHead,fixed = TRUE) == targets)) { ## confirm that guess
                      
                          colHead <- gsub(".1","",colHead,fixed=TRUE)
                          
                          colHead <- formatHead(colHead)
                      
                      } else { ## the mysterious fix -- strip everything 
                        
                        colhead <- formatHead(gsub("[.][1]|[.][2]","",colHead))
                        
                      }
                  }
             }
    
        }
  
  colHead
}

################
## BDtoPyPop
#'Convert BIGDAWG datasets to PyPop datasets
#'  
#'Converts a BIGDAWG-formatted data frame into a pair of PyPop-formatted case and control data frames.
#'
#'@param dataset A data frame containing a BIGDAWG-formatted case-control dataset
#'@param filename A character string identifying the desired path and name for the generated files or the elements of the returned list object.
#'@param save.file A logical that determines if a pair of files should be written (TRUE) or if a list object should be returned. The default value is TRUE.
#'@param save.path A character string identifying the path in which to write the pair of files when save.file is TRUE. The default value is tempdir().
#'
#'@return When save.file = TRUE, a pair of files named "'filename'.positive.pop" and "'filename'.negative.pop" are generated in the working directory. When save.file = FALSE, a list of two-elements, named "'filename'.positive" and "'filename'.negative", is returned. 
#'
#'@note Files generated by BDtoPyPop are intended for analysis using Pypop version 1.#.#. PyPop is not an R package, and must be installed separately. 
#'
#'@references Lancaster et al. Front Immunol. 2024 Apr 2;15:1378512. https://pubmed.ncbi.nlm.nih.gov/38629078/
#'@references Pappas et al. Hum Immunol. 2016 Mar 77(3):283-287. https://pubmed.ncbi.nlm.nih.gov/26708359/
#'
#'@export
#'
#'@examples
#' HLAdata.PP <- BDtoPyPop(sHLAdata,"BDHLA",FALSE)
#'
BDtoPyPop <- function(dataset, filename, save.file=TRUE,save.path = tempdir()) {
  
  colnames(dataset) <- pypopHeaders(colnames(dataset))
  dataset <- convertAny(dataset)
  
  posSet <- dataset[dataset[,2] %in% c(1,"1"),]
  negSet <- dataset[dataset[,2] %in% c(0,"0"),]
  
  if(save.file) {
    write.table(posSet,paste(save.path,paste(filename,"positive","pop",sep="."),sep=.Platform$file.sep),append = FALSE,quote = FALSE,sep ="\t",na = "****",row.names = FALSE,col.names = TRUE)
    write.table(negSet,paste(save.path,paste(filename,"negative","pop",sep="."),sep=.Platform$file.sep),append = FALSE,quote = FALSE,sep ="\t",na = "****",row.names = FALSE,col.names = TRUE)
  } else {
    
    retList <- list(posSet,negSet)
    names(retList) <- c(paste(filename,"positive",sep="."), paste(filename,"neagtive",sep="."))
    
    retList
    
  }
}
