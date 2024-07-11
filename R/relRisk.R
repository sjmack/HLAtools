## HLA RRisk Calculation -- Steven J Mack July 10, 2024 v 1.4.1

#######################
#' Calculate Relative Risk for Individual Alleles and Genotypes in BIGDAWG-formatted Non-Case-Control Datasets
#' 
#' This function returns a list object containing relative risk, confidence interval and p-value data for the individual alleles and individual genotypes at each locus in a BIGDAWG-formatted non-case-control genotype data frame or file. 
#' 
#' @param dataset A character string describing the name of a non-case-control genotype dataset using the BIGDAWG format. Here, "non-case-control" means that while two subject categories are required, the categories should not be patients and controls; instead, the categories may be, e.g., for a dataset of patients, either of two disease states, where one disease state is coded as 0 and the other is coded as 1 in the second column of the dataset. Either a tab-delimited file or a data frame can be specified. 
#' @param return A logical identifying if the list object should be returned (return=TRUE), or if pairs of tab-delimited text files of results (one for alleles and one for genotypes) should be written to the working directory for each locus.
#' @param save.path A character string identifying the path in which to write the pair of files when return is FALSE. The default value is tempdir().
#' 
#' @keywords relative risk genotype allele
#' 
#' @importFrom fmsb riskratio
#' 
#' @importFrom utils read.table capture.output write.table
#' 
#' @return A list object of two lists ("alleles" and "genotypes"), each of which contains a list of nine-column data frames containing results for each unique allele or genotype (in rows) at each locus. Column headers in each dataframe are, *Locus*, *Variant*, *Status_1*, *Status_0*, *RelativeRisk*, *CI.low*, *CI.high*, *p.value*, and *Significant*.
#' 
#' @export
#' 
#' @examples
#' rr <- relRisk(sHLAdata)
#' 
#' @references \href{https://CRAN.R-project.org/package=BIGDAWG/vignettes/BIGDAWG.html#input-data}{BIGDAWG Data Format}
#' 
relRisk <- function(dataset,return=TRUE,save.path = tempdir()){ ## if return == TRUE, a list object is returned; if return==FALSE, two tab-delimited text files are written for each locus; one for alleles and one for genotypes

      if(!is.data.frame(dataset)){
      # load in HLA BIGDAWG formatted data
        HLAgenoData <- read.table(dataset,header=TRUE,sep="\t",quote = "",stringsAsFactors = FALSE,colClasses = "character",check.names = FALSE, na.strings = c("****","NA","","-","na","Na"))
      } else {
              HLAgenoData <- dataset }
  
        colnames(HLAgenoData) <- sub("\\.\\d|\\_\\d","",colnames(HLAgenoData))
  
        loci <- unique(colnames(HLAgenoData)[3:ncol(HLAgenoData)])

        ## sorting genotypes in ascending order
        for(j in 1:nrow(HLAgenoData)) {
              for(i in 1:length(loci)){
                  if(!NA %in% HLAgenoData[j,colnames(HLAgenoData) %in% loci[i]]) { #ignore NAs
              HLAgenoData[j,colnames(HLAgenoData) %in% loci[i]] <- sort(unlist(HLAgenoData[j,colnames(HLAgenoData) %in% loci[i]]),decreasing = FALSE)
              }
        }
    
    }
        cases <- HLAgenoData[HLAgenoData[,2]=="1",] # status_1
        controls <- HLAgenoData[HLAgenoData[,2]=="0",] # status_0

        ## variant level RRs 
        RRtab <- masterTab <- list("alleles"=list(),
                           "genotypes"=list())
        for(k in 1:2){ ## looping through alleles RRtab[[1]], and genotypes RRtab[[2]]

            for(i in 1:length(loci)){ # looping through loci

            if(k == 2){
            ## Genotype evaluations
            masterTab[[k]][[i]] <- merge(as.data.frame(table(apply(cases[,colnames(cases) %in% loci[i]],1,paste,collapse="+")),stringsAsFactors = FALSE),as.data.frame(table(apply(controls[,colnames(controls) %in% loci[i]],1,paste,collapse="+")),stringsAsFactors = FALSE),by.x="Var1",by.y="Var1",all=TRUE)
                } else {
                    ## Allele evaluations
                      masterTab[[k]][[i]] <- merge(as.data.frame(table(unlist(cases[,colnames(cases) %in% loci[i]])),stringsAsFactors = FALSE),as.data.frame(table(unlist(controls[,colnames(controls) %in% loci[i]])),stringsAsFactors = FALSE),by.x="Var1",by.y="Var1",all=TRUE)
                    }

                colnames(masterTab[[k]][[i]]) <- c(loci[i],"Cases","Controls")
                masterTab[[k]][[i]][is.na(masterTab[[k]][[i]])] <- 0
                masterTab[[k]][[i]] <- masterTab[[k]][[i]][!masterTab[[k]][[i]][,1] %in% "NA+NA",] # Remove any NA+NA genotypes
                masterTab[[k]][[i]] <- masterTab[[k]][[i]][substr(masterTab[[k]][[i]][,1],1,3) != "NA+",] # Remove any NA+X genotypes (should't be any ideally)

                totCases <- sum(masterTab[[k]][[i]]$Cases)
                totControls <- sum(masterTab[[k]][[i]]$Controls)

# Matrix for Relative Risk via fmsb:riskratio()
#  Allele | # cases | # controls
# --------+---------+-----------
#    A    |    a    |     b
#  not-A  |    c    |     d
#
# riskratio() parameters
#    X = a
#    Y = c = (totCases - a)
#   m1 = a+b 
#   m2 = c+d = (totCases - a) + (totControls - b) 

            RRtab[[k]][[i]] <- data.frame(Locus=character(),
                 Variant=character(),
                 Status_1=character(),
                 Status_0=character(),
                 RelativeRisk=numeric(), 
                 CI.low=numeric(), 
                 CI.high=numeric(),
                 p.value=numeric(),
                 Significant=character(),
                 stringsAsFactors=FALSE)

            for(j in 1:nrow(masterTab[[k]][[i]])){ # looping through variants
                  ## This is sort of kludgy -- in some instances when these values go directly into riskratio, the p-value is not calculated, and a warning:
                  ## "NAs produced by integer overflow" is produced, but by converting these to characters and then back into integers, that goes away.
                  x <- as.character(masterTab[[k]][[i]][j,2])
                  y <- as.character(totCases-masterTab[[k]][[i]][j,2])
                  M1 <- as.character(masterTab[[k]][[i]][j,2]+masterTab[[k]][[i]][j,3])
                  M2 <- as.character((totCases-masterTab[[k]][[i]][j,2])+(totControls-masterTab[[k]][[i]][j,3]))
  
              silenceRR <- capture.output(RR <- riskratio(as.numeric(x),as.numeric(y),as.numeric(M1),as.numeric(M2)))
              RRtab[[k]][[i]][j,1:9] <- c(loci[i],masterTab[[k]][[i]][j,1],masterTab[[k]][[i]][j,2],masterTab[[k]][[i]][j,3],RR$estimate,RR$conf.int[1],RR$conf.int[2],RR$p.value,ifelse(RR$p.value<0.05,"*",""))
                    } # end of variant loop j  

              if(!return){write.table(RRtab[[k]][[i]],file=paste(save.path,paste(loci[i],ifelse(k == 2,"genotype","allele"),"relative_risk.txt",sep="."),sep=.Platform$file.sep),append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE,col.names = TRUE)
                    message(loci[i],ifelse(k == 2,".genotype",".allele"),".relative_risk.txt file written.")
                      }
                } # end of Loci loop i 
            } # end of allele/genotype loop k
        if(return){return(RRtab)}
      } # end of Function
