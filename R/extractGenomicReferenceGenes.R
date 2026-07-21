## ExtractGenomicReferences v1.0.0 July 17, 2026 

# BuildGenomicReferences
#'Generate a Table of Genomic Reference Alleles and Characteristics of HLA Genes Defined by the IPD-IMGT/HLA Database 
#'
#'@description
#'buildGenomicReferences() extracts the table of genomic HLA reference alleles at 'www.ebi.ac.uk/ipd/imgt/hla/about/help/genomics/' and generates the IMGTGenomicReferences data object. 
#'
#'@return A list object of two elements -- GenomicReferences and Version. 
#'
#'@importFrom rvest read_html html_elements html_table
#'@importFrom stringr %>%
#'
#'@export
#'
#'@note This function requires internet access to function.
#'@note For internal HLAtools use.
#'
#'@references Barker et al. NAR (2025) 54:D1152–D1158 https://doi.org/10.1093/nar/gkaf1218
buildGenomicReferences <- function() {      
    
    url <- "https://www.ebi.ac.uk/ipd/imgt/hla/about/help/genomics/"
    
             page <- read_html(url)
    
          # Extract all HTML tables (currently #2 is the genome reference table)
            tables <- page %>%
             html_elements("table") %>%
                 html_table(fill = TRUE)
    
    genomicRefTab <- as.data.frame(tables[[2]])
    
    generationDate <- strsplit(as.character(Sys.time()),fixed = T," ")[[1]][1]
    
    IGR <- list(genomicRefTab,generationDate)
    names(IGR) <- c("GenomicReferences","version")
    
    IGR
}
