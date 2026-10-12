## build WMDA Serological Relationships version 1.0 11AUG2026

################
#buildWMDASerologicallRelationships
#'Build local copies of the rel_ser_ser and rel_dna_ser tables.
#'
#'Consumes the ANHIG/IMGTHLA GitHub repository's rel_dna_ser.txt and rel_ser_ser.txt files to create the WMDASerologicalRelationships object.
#'
#'@param version A character string describing the desired release version (branch) of the ANHIG/IMGTHLA Github repository (e.g. '3.53.0'). The default value ('Latest') returns a list of 'rel_ser_ser', 'rel_dna_ser' and 'version' objects for the most recent release version.
#'@param asterisk A boolean indicating if the asterisk in the first column of the rel_dna_ser.txt file should be included in (TRUE) or excluded from (FALSE) the generated rel_dna_ser object. The default value is TRUE.
#'
#'@returns A list object of three elements, the rel_ser_ser table, the rel_dna_ser table, and a version vector identyfying the IPD-IMGT/HLA Database release version under which this object was generated. When 'version' is 3.63.0 or lower, the rel_dna_ser table will include six columns. When 'version' us 3.64.0 or higher, the rel_dna_ser table will include a seventh 'HATS' column.
#'
#'@importFrom dplyr case_when
#'@importFrom utils read.csv
#'
#'@export
#'
#'@source The rel_ser_ser object is derived from the wmda folder in the specified version of the IPD-IMGT/HLA GitHub repository at https://github.com/ANHIG/IMGTHLA/tree/Latest/wmda.
#'@note The names of columns 3-7 in the $rel_dna_ser table are abbreviated: 'USA' stands for 'Unambiguous Serological Antigen', 'PSA' stands for 'Possible Serological Antigen', 'ASA' stands for 'Assumed Serological Antigen', 'Expert' stands for 'Expert assigned exceptions in search determinants of some registries', and HATS stands for 'HLA allele to serotype'.  
#'@note Details about the HATS application can be found at https://github.com/kosoegawa/HATS and https://onlinelibrary.wiley.com/doi/10.1111/tan.14662.
#'@note The column names in the $rel_ser_ser table identify the HLA locus (Locus), HLA antigen name (Antigen), split antigen (Split), and associated antigen (associated).
#'
buildWMDASerRel <- function(version = "Latest",asterisk = TRUE) {        
        
          if(version == "Latest") { version <- getLatestVersion()}
  
          if(!validateVersion(version)) {
             message(paste("Version",version,"does not exist.",sep=" "))
              return(NA)
            }

          version <- squashVersion(version)
          
    ### Version-specific parameters
    ## version 3.00.0 to 3.63.0
        # rel dna ser has 6 columns -- 4 line header 3.0.0-3.21.0; 5 line header 3.22.0-3.24.0; 6 line 3.25.0-3.63.0
        # rel ser ser has 4 columns -- 4 line header 3.0.0-3.21.0; 5 line header 3.22.0-3.24.0; 6 line 3.25.0-3.63.0
        # release version not present in releases 3.0.0 to 3.21.0; present on line 3 in release 3.22.0 onward
  
    ## version 3.64.0 and up 
        # rel dna ser has 7 columns -- 6 line header 3.64.0-3.65.0 (current)
        # rel ser ser has 4 columns -- 6 line header 3.64.0-3.65.0 (current)
 
               num.version <- as.numeric(squashVersion(version))
               
               version.present <- num.version > 3210 ## determine if the release version appears in the header
           
               body.start <- case_when(
                 num.version < 3220 ~ 5,
                 num.version >= 3220 && num.version < 3250 ~ 6,
                 num.version >= 3250 ~ 7
               )
               
               num.columns <- ifelse(num.version > 3630,7,6) ## for rel dna ser only
  
            relDNAser <- readLines(paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/refs/heads/",version,"/wmda/rel_dna_ser.txt",sep=""),-1,warn = F)
                DNAversion <- ifelse(version.present,strsplit(relDNAser[3],split = " ",fixed = T)[[1]][4],version) ## pull version from the document when available
                relDNAser.tab <- read.csv(text=relDNAser[body.start:length(relDNAser)],sep = ";",header = F)
                if(num.columns == 7) {
                    colnames(relDNAser.tab) <- c("Locus","Allele","USA","PSA","ASA","Expert","HATS")
                        } else {
                          colnames(relDNAser.tab) <- c("Locus","Allele","USA","PSA","ASA","Expert")
                                }
                if(!asterisk) {relDNAser.tab$Locus <- sub(".$","",relDNAser.tab$Locus)}

            relSERser <- readLines(paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/refs/heads/",version,"/wmda/rel_ser_ser.txt",sep=""),-1,warn = F)
                relSERser.tab <- read.csv(text=relSERser[body.start:length(relSERser)],sep = ";",header = F)
                colnames(relSERser.tab) <- c("Locus","Antigen","Split","Associated")

                  wSR <- vector(mode = "list",length = 3)

                    wSR[[1]] <- relSERser.tab
                    wSR[[2]] <- relDNAser.tab
                    wSR[[3]] <- DNAversion
      
                  names(wSR) <- c("rel_ser_ser","rel_dna_ser","version")
                  
                  wSR
      }
