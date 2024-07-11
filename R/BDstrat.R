#### BDstrat -- multi-allele stratification in BIGDAWG formatted datasets  
####            v2.1.1 SJM July 10, 2024 

####################
##BDstrat
#'Stratify BIGDAWG Datasets by Specific Alleles
#'
#'Divides a BIGDAWG-formatted dataset into two subsets (strata) that do and do not include specified alleles.
#'
#'@param dataset A BIGDAWG-formatted data frame or a path to a BIGDAWG-formatted, tab-delimited text file.
#'@param alleles A vector of allele names in the locus-asterisk-allele_name format (e.g., "A*01:01:01:01").
#'@param warnBelow An integer that defines a low number of subjects in a stratum, generating a warning message. The default value is 21.
#'
#'@return A list-object of two BIGDAWG-formatted data frames titled dataset$`<alleles>-positive` and dataset$`<alleles>-negative`. The positive list element includes all subjects with the specified alleles, and the negative list element includes all subjects without those specified alleles.
#'
#'@examples
#' HLA_data.multi.strat <- BDstrat(sHLAdata,c("DRB1*16:02:01:01","DRB1*04:07:01:01","A*25:01:01:01"))
#' HLA_data.single.strat <- BDstrat(sHLAdata,"DRB1*16:02:01:01")
#'
#'@export
#'
#'@references \href{https://CRAN.R-project.org/package=BIGDAWG/vignettes/BIGDAWG.html#input-data}{BIGDAWG Data Format}
#' 
BDstrat <- function(dataset,alleles,warnBelow = 21){  
  
   ## loading the dataset if it is a path
  if(!is.data.frame(dataset)) { 
    dataset <- read.table(dataset,header = TRUE,sep="\t",check.names = TRUE,stringsAsFactors = FALSE)
    # Since all loci are duplicated, check.names = TRUE generates "locus","Locus.1" name pairs if they are not unique
  } 
  
  ## checking for locus suffixes in the dataset's column header
  suffix <- NA
  suffVec <- grepl("_",colnames(dataset)[3:length(colnames(dataset))],fixed=TRUE)
  
  #no mixture of suffixes allowed
  if(length(unique(suffVec)) == 2) { stop("Each locus name must either be sufffixed with '_1' and '_2, or have no suffix.") }
  
  if(unique(suffVec) == FALSE) {suffix = FALSE} else {suffix = TRUE} ## set suffixed or not 
  
  # and make.names(x,unique=TRUE) does the same for data frames
  if(!suffix) { colnames(dataset) <- make.names(colnames(dataset),unique = TRUE) }
   
  datafalse <- as.data.frame(matrix(data=FALSE,nrow=nrow(dataset),ncol = ncol(dataset)))
  colnames(datafalse) <- colnames(dataset)
  allcols <- c()
  
  for(s in 1:length(alleles)) { # main loop for all of the stratified alleles
    
        allele <- c(strsplit(alleles[s],"*",fixed = TRUE)[[1]])
                                                                                                                  
        if(suffix) {
          currLoc <- (c(paste(allele[1],"1",sep="_"),paste(allele[1],"2",sep="_"))) } else {
            currLoc <- make.names(c(allele[1],allele[1]),unique=TRUE) }
                                                                                                                  
        allcols <- append(allcols,currLoc) ## the targeted column names
                                                                 
        for(t in 1:nrow(datafalse)) {               
          
          if(!is.na(dataset[t,colnames(dataset) %in% currLoc[1]])) {
                if(dataset[t,colnames(dataset) %in% currLoc[1]] == allele[2]) { datafalse[t,colnames(dataset) %in% currLoc[1]] <- TRUE}
                      } #end NA check
          if(!is.na(dataset[t,colnames(dataset) %in% currLoc[2]])) {
                if(dataset[t,colnames(dataset) %in% currLoc[2]] == allele[2]) { datafalse[t,colnames(dataset) %in% currLoc[2]] <- TRUE}
                      } #end NA check
                } # end of t loop
            
          }  # end of s loop
  
  datafalse <- datafalse[,colnames(datafalse) %in% unique(allcols)]
  
  sumVec <- rep(0,nrow(datafalse))
  
      for(i in 1:nrow(datafalse)){
          sumVec[i] <- sum(datafalse[i,])
        }

  posStrat <- dataset[sumVec[] > 0,]
  negStrat <- dataset[sumVec[] == 0,]
 
  # Add them as elements of the stratPair list, named with the selected or excluded alleles
  stratPair <- list()
  stratPair[[paste(paste(alleles,collapse="+"),"-positive",sep="")]] <- posStrat
  stratPair[[paste(paste(alleles,collapse="+"),"-negative",sep="")]] <- negStrat
  
  for(i in 1:2){
    # Strip suffixes from column names if they have any
    colnames(stratPair[[i]]) <- gsub(".1","",colnames(stratPair[[i]]),fixed=TRUE)
    
    # Silently eliminate empty rows
    stratPair[[i]][is.na(stratPair[[i]])] <- ""
    stratPair[[i]] <- stratPair[[i]][!substr(rownames(stratPair[[i]]),1,2) == "NA",]
    
    # Return Phenotype to integer
    stratPair[[i]][,2] <- as.integer(stratPair[[i]][,2])
    
    # Make blanks into NAs
    stratPair[[i]][stratPair[[i]][,] == ""] <- NA
    
    ## alert if number of rows in a stratum is < 20 (an arbitrary cutoff, to be adjusted)
    if(nrow(stratPair[[i]])<warnBelow){warning(paste("The",names(stratPair[i]),"stratum contains",nrow(stratPair[[i]]),ifelse(nrow(stratPair[[i]])==1,"subject.","subjects."),sep=" "))}
  }
  #return object
  stratPair
  
}
