#### BDstrat -- single allele stratification in BIGDAWG formatted datasets  
####            v1.0.0 SJM March 24,2024

####################
##BDstrat
#'Allele stratification for BIGDAWG datasets
#'
#'Divides a BIGDAWG-formatted dataset into two subsets (strata) that do and do not include a specific allele.
#'
#'@param dataset A BIGDAWG-formatted data.frame or a path to a BIGDAWG-formatted tab-delimited text file.
#'@param allele A single allele in the <locus>*<allele name> format (e.g., "A\*01:01:01:01" for HLA alleles or "KIR2DL1\*001" for KIR alleles).
#'@param warnBelow An integer that defines a low number of subjects in a stratum, generating a warning message. The default value is 21.
#'
#'@return A list-object of two BIGDAWG-formatted data frames titled dataset$`allele-positive` and dataset$`allele-negative`. The positive list element includes all subjects with a specified allele, and the negative list element includes all subjects without that specified allele.
#'
#'@examples
#'\dontrun{
#' HLA_data.strat <- BDstrat(BIGDAWG::HLA_data,"DRB1*15:01:01:01",Run.Tests = "L",HLA=TRUE)
#' for(i in 1:2) {BIGDAWG(HLA_data.strat[[i]],HLA = TRUE,Run.Tests = "L")}
#'}
#'
#'@export
#'
BDstrat <- function(dataset,allele,warnBelow = 21){
  
  ## parsing allele name 
  tempName <- strsplit(allele,"*",fixed=TRUE)
  locus <- tempName[[1]][1]
  alleles <- tempName[[1]][2]
  loaded <- FALSE
  
  ## loading the dataset if it is a path
  if(!is.data.frame(dataset)) { 
    dataset <- read.table(dataset,header = TRUE,sep="\t",check.names = TRUE,stringsAsFactors = FALSE)
    # Since all loci are duplicated, check.names = TRUE generates "locus","Locus.1" name pairs if they are not unique
    loaded <- TRUE
    } 
  
  ## checking for locus suffixes in the dataset's column header
  suffix <- NA
  suffVec <- grepl("_",colnames(dataset)[3:length(colnames(dataset))],fixed=TRUE)
  
  #no mixture of suffixes allowed
  if(length(unique(suffVec)) == 2) { stop("Each locus name must either be sufffixed with '_1' and '_2, or have no suffix.") }
  
  if(unique(suffVec) == FALSE) {suffix = FALSE} else {suffix = TRUE} ## set suffixed or not 
  
  # suffxing the selected locus to match the dataset
  if(!suffix) {
      locus <- c(locus,paste(locus,"1",sep=".")) } else {
      locus <- c(paste(locus,"1",sep="_"),paste(locus,"2",sep="_"))
      }
  
  # and make.names(x,unique=TRUE) does the same for data frames
  if(!suffix) { colnames(dataset) <- make.names(colnames(dataset),unique = TRUE) }
  
  # make sure the locus is in the dataset
  if(!locus[1] %in% colnames(dataset)) {
    stop(paste(locus[1],"is not found in the dataset.",sep=" "))
  }
  
  #  make sure the allele is in the dataset
  locus.alleles <-  unique(unique(dataset[,colnames(dataset) %in% locus[1]]),unique(dataset[,colnames(dataset) %in% locus[2]]))
  if(!alleles %in% locus.alleles) {stop(paste(allele,"is not really found in the dataset.",sep=" "))}
  
  # Everything is stored in the strataSet list 
  stratSet <- data.frame(rep(NA,nrow(dataset)))

  # Identify the rows containing each target allele for each locus column
  for(i in 1:length(alleles)){
      stratSet <- cbind(stratSet,dataset[[locus[1]]]==alleles[i],dataset[[locus[2]]]==alleles[i])
      }

  # Identify the rows containing any target allele
  for(i in 1:nrow(stratSet)){
      stratSet[i,1] <- any(unlist(stratSet[i,2:((2*length(alleles))+1)]))
      }

  # Split the parent dataset into two stratified subsets
  posStrat <- dataset[stratSet[,1]==TRUE,]
  negStrat <- dataset[stratSet[,1]==FALSE,]

  # Add them as elements of the stratPair list, named with the selected or excluded alleles
  stratPair <- list()
  stratPair[[paste(strsplit(locus[1],"_",fixed = TRUE)[[1]][1],"*",paste(unlist(alleles),collapse=paste("+",locus,"*",sep="")),"-positive",sep="")]] <- posStrat
  stratPair[[paste(strsplit(locus[1],"_",fixed = TRUE)[[1]][1],"*",paste(unlist(alleles),collapse=paste("+",locus,"*",sep="")),"-negative",sep="")]] <- negStrat

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
