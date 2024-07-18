#### GIANT -- version 1.0.0 July 17, 2024

################
##GIANT
#'GLupdate-Integrated Allele Name Translation
#'
#'Calls GLupdate() to translate the HLA allele names in a vector or data frame between two IPD-IMGT/HLA Database release versions.
#'
#'@param data A data frame or vector of HLA allele name character string data. Allele names in vectors must include locus prefixes. Data frames must include column headers identifying a single locus for the allele names in that column. 
#'@param transFrom A dot-formatted character string identifying the IPD-IMGT/HLA Database version (e.g. 3.57.0) under which the HLA data were generated.
#'@param transTo A dot-formatted character string identifying the IPD-IMGT/HLA Database version (e.g. 2.25.0) to which which the HLA data are to be translated.
#'
#'@return A data frame or vector of HLA allele name data in the 'transTo' release. 
#'
#'@note GIANT() will ignore the first two columns of BIGDAWG-formatted data frames, which include non-HLA data. All columns in non-BIGDAWG-formatted data frames are expected to contain HLA allele names.
#'
#'@export
#'
#'@examples
#'GIANT(c("A*01:01:01:01","DOA*01:01:01:01"),"3.56.0","2.20.0")
#'updsHLAdata <- GIANT(sHLAdata,"3.56.0","2.25.0")
#'
GIANT <- function(data,transFrom,transTo) {

if(!validateVersion(transFrom)) {return(paste(transFrom,"is not a valid release version. Please provide a valid release version for 'transFrom'.",sep=" "))}
if(!validateVersion(transTo)) {return(paste(transTo,"is not a valid release version. Please provide a valid release version for 'transTo'.",sep=" "))}

### figure out data parameters 
        if(is.data.frame(data)) { ## determine type of data 
                  df.data <- TRUE # is a data frame 
                  is.vector <- FALSE
                  startCol <- 1 # assume all columns of a vector contain data 
                  prefix = ifelse(identical(grep("*",data[,3],fixed=TRUE),integer(0)),FALSE,TRUE) # determine if the alleles are prefixed
                  if(any(c(0,1) %in% data[,2])) {
                                bd.data <- TRUE
                                startCol <- 3 # ignore the first two columns of BD data
                                } 
                            } else {
                                  df.data <- FALSE
                                  is.vector <- TRUE} # is a vector
 ### for data.frames
          if(df.data){
                  
            fullframe <- data.frame(matrix(NA,nrow=nrow(data),ncol=ncol(data))) ## this dataframe contains full (suffixed) allele names 
            fullframe[,1:2] <- data[,1:2] ## if this is a BIGDAWG formatted dataset, otherwise it gets overwritten anyway when startCol = 1, below
                  
                  for(i in startCol:ncol(data)){
                    
                        if(!prefix) {   # Add prefixes if necessary
                            fullframe[,i] <- paste(colnames(data)[i],data[,i],sep="*")
                                } else {
                           fullframe[,i] <- data[,i]
                                  }
                    
                            updVec <- GIANT(fullframe[,i],transFrom,transTo)
                            
                            if(!prefix) {
                              for(j in 1:length(updVec)) {
                                      ifelse(is.na(updVec[j]),fullframe[j,i] <- NA,fullframe[j,i] <- strsplit(updVec[j],"*",TRUE)[[1]][2])
                                }
                            } else {
                              fullframe[,i] <- updVec
                            }
                                      
                        colnames(fullframe) <- colnames(data)
                     
                    } 
            
            return(fullframe)
            
            # end of data frame section
        
 ### for vectors        
                } else { 
                        uniqueAlleles <- unique(data) 

          ignore <- grep("NA",uniqueAlleles,fixed=TRUE)

          uniqueAlleles <- uniqueAlleles[!1:length(uniqueAlleles) %in% ignore]

          #### translate the unique alleles

          GLScode <- paste("hla",transFrom,paste(uniqueAlleles,collapse="~"),sep="#")

          updAlleles <- gsub("HLA-","",strsplit(strsplit(updateGL(GLScode,transTo),"#",fixed = TRUE)[[1]][3],"~",fixed=TRUE)[[1]],fixed=TRUE)

                    if(length(updAlleles) != length(data)){
                      
                      retVec <- rep(NA,length(data))
                            
                      transTab <- rbind(uniqueAlleles,updAlleles)
                    
                          for(i in 1:length(data)) {
                        
                            if(any(transTab[1,] == data[i])) { 
                              retVec[i] <- transTab[2,transTab[1,] == data[i]]
                                  } else {
                              retVec[i] <- NA }
                          }
                      
                          return(retVec)
            
                            } else {
            
                        return(updAlleles) }
            
              } # end of vector section
          } # end
