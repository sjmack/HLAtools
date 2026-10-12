 ######### GenomicReferences Functions v.1.0 4OCT2026 SJM

 ################
 #genomicReferenceLookup
 #Search the IMGT Genomic References Table to Retrieve Genomic Reference Data for Specific IMGT/HLA Genes
 #'
 #'Searches the IMGTGenomicReferences$GenomicReferences table to return genomic reference data for specific IPD-IMGT/HLA genes.
 #'
 #'@param locus A string identifying the gene or genes for which genomic reference data is desired. All loci can be specified by setting 'locus' to "all". The default 'locus' value is "all".
 #'@param column A vector of the desired column or columns to look-up for 'locus'. All columns can be specified by setting 'column' to "all". Column names have been simplified for searching, as detailed in the note. The default 'column' value is "all". 
 #'
 #'@return A data.frame identifying the values in 'column' for each gene in 'locus'. 
 #'
 #'@export
 #'
 #'@note The names of column's 2-6 in the GenomicReferences table have been shortened to facilitate searches -- "Ensemblversion" is "Ensembl.ver", "GRCh38 Location" is "GRCh38.loc", "Start of ATG (Initiator Met)" is "ATG", "GRCh38 reference allele" is "GRCh38.ref", and "IPD-IMGT/HLA reference allele" is "IMGT.ref". The full-length column names are included in the returned data.frame.
 #'@note The IMGTGenomicReferences table is available at https://www.ebi.ac.uk/ipd/imgt/hla/help/genomics.html.
 #'
 genomicReferencesLookup <- function(locus="all",column="all"){
  
            locus <- unique(locus)
            column <- unique(column)
            
            # append a locus column and shorten native column names
            splitNames <- strsplit(IMGTGenomicReferences$GenomicReferences$Gene,split = "-",fixed = T)
            singleNames <- rep("",length(splitNames))
        
            for(i in 1:length(splitNames)) { # append locus column
                 if(length(splitNames[[i]])==1) {
                   singleNames[i] <- splitNames[[i]][1]
                     } else {singleNames[i] <- splitNames[[i]][2]
                         } 
                            }
            
            RGRefs <- cbind(singleNames,IMGTGenomicReferences$GenomicReferences)
            colnames(RGRefs)[1] <- "Locus"
            orig.colnames_3_7 <- colnames(RGRefs)[3:7] 
            colnames(RGRefs)[3:7] <- c("Ensembl.ver","GRCh38.loc","ATG","GRChg38.ref","IMGT.ref")
            
            if(any(c("all","") %in% locus)) {
                  locus <- RGRefs$Locus
                  }
            
            if(any(c("all","") %in% column)) {
                  column <- colnames(RGRefs)[3:7]  
                }
            
            # remove HLA- prefixes for locus look-up
            hla.prefix <- which(substr(locus,1,4) == "HLA-")
            
                if(!length(hla.prefix) == 0) {
                  
                  for(i in 1:length(hla.prefix)) {
                          locus[hla.prefix[i]] <- strsplit(locus[hla.prefix[i]],"-",TRUE)[[1]][2]
                      }
                }

            # remove invalid locus and column values 
            if(!all(locus %in% RGRefs$Locus)) {
              
              loc.err <- locus[!locus %in% RGRefs$Locus]
                
              if(length(loc.err) > 0) {
                      message("The ", ifelse(length(loc.err) == 1, paste(loc.err, "gene is not in the 'GenomicReferences' table.", collapse = " "),
                                             paste(paste(loc.err[1:length(loc.err) -1],collapse = ", "), " and ", loc.err[length(loc.err)], " loci are not found in the 'GenomicReferences' table.",sep="")))
             
                      locus <- locus[!locus %in% loc.err]
                      } 
                   }
            
              col.err <- column[!column %in% colnames(RGRefs)]
            
              if(length(col.err) > 0) {
                      message("The ",ifelse(length(col.err) == 1, paste(col.err, "column is not in the 'GenomicReferences' table.", collapse = " "),
                                    paste(paste(col.err[1:length(col.err)-1],collapse = ", ")," and ", col.err[length(col.err)]," columns are not found in genomicReferneces.",sep="")))
                
                    column <- column[!column %in% col.err]
               }
              
              # generate the output 
                retDF <- RGRefs[RGRefs$Locus %in% locus,colnames(RGRefs) %in% c("Gene",column),drop=FALSE]
            
                colnames(retDF)[2:ncol(retDF)] <- orig.colnames_3_7[colnames(RGRefs)[3:7] %in% colnames(retDF)]
                
                rownames(retDF) <- 1:nrow(retDF)
                
                retDF
        }
   
   
   