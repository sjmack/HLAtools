### Extract IPD-IMGT/HLA Database-Supported Gene Types v3.0.0 16 Apr 2024

####################
## BuildIMGTHLAGeneTypes
#'Describe IPD-IMGT/HLA Database Genes, Identifying Pseudogenes and Gene Fragments
#'
#'buildIMGTHLAGeneTypes() scrapes 'hla.alleles.org/pages/genes/genes_list' and generates a data frame identifying each gene. It requires internet access to function. As such, this function always returns data for the current IPD-IMGT/HLA Database release.
#'
#'@return A list object of two elements -- 'version' and 'GeneTypes'. The 'version' element identifies the date on which the source table was generated. The 'GeneTypes' element is a data frame of three columns, identifying each gene supported by the IPD-IMGT/HLA Database, along with its molecular characteristics and its status as either a pseudogene or a gene fragment.
#'
#'@importFrom rvest read_html html_table html_nodes
#'@importFrom stats setNames
#'@importFrom stringr %>%
#'
#'@export
#'
#'@note For internal HLAtools use.
#'
buildIMGTHLAGeneTypes <- function(){
  url <-  "https://hla.alleles.org/pages/genes/genes_list"
  
    # Setup and version extraction
    rawPage <- readLines(url,-1,warn = FALSE)
    updateLine <- which(rawPage[] == "        <p class=\"last-updated text-right content\">")+1 ## the line with the update information on it. 
    version <- strsplit(rawPage[updateLine],": ")[[1]][2]

          ## Extract the molecular characteristics table from the page HTML
          df <- url |> 
          read_html() |> 
          html_nodes("table") |> 
          html_table(fill = T) %>%
          lapply(., function(x) setNames(x, c("Names", "Equivalent", "Alignments", 
                                        "Molecular Characteristics","Function","Links")))
  
          ## Capture columns 1 and 4, and ad a pseudogene/genefragment column
          baseTab <- as.data.frame(df[[1]][,c(1,4)])
          baseTab <- cbind(baseTab,rep("",nrow(baseTab)))
          colnames(baseTab)[3] <- "Pseudogene/Fragment"  
  
          ## Identify Pseudogenes and Gene Fragments
          pseudoRows <- grep("pseudogene",baseTab[,2],fixed = FALSE)
          fragmentRows <- grep("gene fragment",baseTab[,2],fixed = FALSE)

          baseTab[,3][pseudoRows] <- "Pseudogene"
          baseTab[,3][fragmentRows] <- "Fragment"
  
    geneTypes <- list(baseTab,version)
    names(geneTypes) <- c("GeneTypes","version")
    
    return(geneTypes)
}
