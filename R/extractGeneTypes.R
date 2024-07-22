### Extract IPD-IMGT/HLA Database-Supported Gene Types v2.0.1 20 March 2024

####################
## BuildIMGTHLAGeneTypes
#'Describe IPD-IMGT/HLA Database Genes, Identifying Pseudogenes and Gene Fragments
#'
#'buildIMGTHLAGeneTypes() scrapes 'hla.alleles.org/genes/index.html' and generates a data frame identifying each gene. It requires internet access to function. As such, this function always returns data for the current IPD-IMGT/HLA Database release.
#'
#'@return A list object of two elements -- 'version' and 'GeneTypes'. The 'version' element identifies the date that the source table was generated. The 'GeneTypes' element is a data frame of three columns, identifying each gene supported by the IPD-IMGT/HLA Database, along with its molecular characteristics and its status as either a pseudogene or a gene fragment.
#'
#'@export
#'
#'@note For internal HLAtools use.
#'
buildIMGTHLAGeneTypes <- function(){
  rawPage <- readLines("https://hla.alleles.org/genes/index.html",-1,warn = FALSE)
  start <- which(rawPage[] == "          <th>Molecular characteristics</th>") +3
  end <- which(rawPage[] == "      </table>") -2
  dateMade <- rawPage[grep("#BeginDate",rawPage,fixed=TRUE)][1]
  dateMade <- strsplit(strsplit(dateMade,"-->",fixed = TRUE)[[1]][2],"<!--",fixed=TRUE)[[1]][1]
  rawPage <- rawPage[start:end]
  trs <- which(rawPage[] %in% "        <tr>")
  rawPageTable <- data.frame(rawPage[!1:length(rawPage) %in% c(trs, trs-1)])
  PageTable <- as.data.frame(matrix(NA,nrow(rawPageTable)/3,3))

  j <- 1
  for(i in 1:nrow(rawPageTable)) {
    PageTable[ceiling(i/3),j] <- rawPageTable[i,1]
    j <- j+1
    if(j == 4){j <- 1}
  }

  PageTable <- cbind(PageTable[,1],PageTable[,3]) ## Dump the middle column of old names
  PageTable[,1] <- str_replace(PageTable[,1],pattern = "          <td><em>",replacement = "")
  PageTable[,1] <- str_replace(PageTable[,1],pattern = "</span></em></td>",replacement = "")
  PageTable[,2] <- str_replace(PageTable[,2],pattern = "          <td>",replacement = "")
  PageTable[,2] <- str_replace(PageTable[,2],pattern = "</td>",replacement = "")
  PageTable[,2] <- str_replace(PageTable[,2],pattern = "&alpha;",replacement = "Alpha")
  PageTable[,2] <- str_replace(PageTable[,2],pattern = "&beta;",replacement = "Beta")
  colnames(PageTable) <- c("Names","Molecular Characteristics")

  PageTable <- cbind(PageTable,(rep("",nrow(PageTable))))
  pseudoCols <- grep("pseudogene",PageTable[,2],fixed = FALSE)
  fragmentCols <- grep("gene fragment",PageTable[,2],fixed = FALSE)

  PageTable[,3][pseudoCols] <- "Pseudogene"
  PageTable[,3][fragmentCols] <- "Fragment"

  colnames(PageTable)[3] <- "Pseudogene/Fragment"

  PageTable <- as.data.frame(PageTable)
  
  typeList <- list(firstList <- dateMade,  
                  secondList <- PageTable)
  
  names(typeList) <- c("version","GeneTypes")
  
  typeList
}
