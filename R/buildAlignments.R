#BuildAlignments v1.5.0 21JUN2024 - LT/SM

################
##BuildAlignments
#'Build Amino Acid, cDNA and gDNA Alignments
#'
#'Returns a list of data frames of amino acid, codon, nucleotide and genomic alignments, with version information.
#'
#'@param loci A character vector of HLA gene names (e.g., "DRB1", c("A","C")).
#'@param source A character vector of alignment types. The allowed values are "AA", "cDNA", and "gDNA". If 'source' is "cDNA", both codon and cDNA nucleotide alignments are generated. If source is 'AA' or 'gDNA', a single peptide or genomic nucleotide alignment is generated. Up to four alignments will be returned for a locus, as determined by its ability to be transcribed or translated.
#'@param version A character string describing the desired release version (branch) of the ANHIG/IMGTHLA Github repository (e.g. '3.53.0'). The default value ('Latest') returns alignments for the most recent release.
#'
#'@return A list object with a data frame of all allele names (and trimmed allele names) and their corresponding sequences (Amino Acid, codon, cDNA, or gDNA) for a specific locus, as well as version details for the returned information. These alignments identify locations of feature boundaries in relation to amino acid, codon, cDNA, and gDNA sequences.
#'
#'@importFrom stringr str_squish
#'@importFrom tibble add_column
#'@importFrom dplyr filter %>%
#'@importFrom utils head tail capture.output
#'
#'@note For internal HLAtools use.
#'
#'@export
#'
buildAlignments<-function(loci, source, version = "Latest"){
  
  if(version != "Latest"){ #
    if(!validateVersion(version)){stop(paste(version," is not a valid IPD-IMGT/HLA Database release version."))}
  }else{ version <- getLatestVersion()}
  
  source <- checkSource(source)
  
  #checks if input locus is present in version 3.38.0 HLA loci
  #skip name checks for DRB1/3/4/5, as they are a part of the DRB alignment
  for(j in 1:length(loci)){
    if(loci[j]=="DRB1"|loci[j]=="DRB3"|loci[j]=="DRB4"|loci[j]=="DRB5") next
    for(x in 1:length(source)) {
      if(source[x] == "cDNA") {
        if(loci[j]%in% HLAgazeteer$nuc == FALSE) {
          return(warning(paste(loci[j], "is not currently supported for", source[x])))
        }
      }
      if(source[x] == "gDNA") {
        if(loci[j]%in% HLAgazeteer$gen == FALSE) {
          return(warning(paste(loci[j], "is not currently supported for", source[x])))
        }
      }
      if(source[x] == "AA") {
        if(loci[j] %in% HLAgazeteer$prot == FALSE) {
          return(warning(paste(loci[j], "is not currently supported for", source[x])))
        }
      }
    }
  }


  #creates empty variables for future for loops
  comb_list<-start<-end<-alignment<-list()

  #empty variables for correspondence table
  alignmentVersion<-pepsplit<-refexon<-aligned<-final_alignment<-AA_codon_alignments<-DNAalignments<-HLAalignments<-exonB<-inDels<-AA_atlas<-gDNA_atlas<-cDNA_atlas<-atlas<-corr_table<-cols<-downloaded_segments<-w<-alignment_positionsx3<-alignment_positions<-alignment_length<-DNA_start<-alignment_start<-space_diff<-prot_extractions<-refblock_number<-sapply(loci, function(x) NULL)

  #sub out periods in version, if there are any
  # version <- gsub(".", "", version, fixed=T) ## change to 
  version <- repoVersion(version)

  #if version is not latest, turn into numeric object
  if(version != "Latest"){
    version <- as.numeric(version)
  }
  
  for(i in 1:length(loci)){


    for(j in 1:length(source)){
      HLAalignments<-sapply(loci, function(x) NULL)
      #expressions unique to each table source (AA or cDNA)
      if(source[j]=="AA"){
        suffix <- "_prot.txt"
        type <- "Prot"
        delete_lines <- c(1,2)
        divide <- 1
        sequence_name <- "AAsequence" 
        } else if(source[j]=="cDNA"){
        suffix <- "_nuc.txt"
        type <- "cDNA"

        #change delete lines to first and second lines if locus does not have protein sequence ## LT
        if(loci[i] %in% HLAgazeteer$nuc[!HLAgazeteer$nuc %in% HLAgazeteer$prot]){
          delete_lines<-c(1,2)
                              #### Fix: DPA2, DPB2 and HLA-N cDNA alignments included an 'AA codon' row in several releases, through they are both pseudogenes
                              if(loci[i] == "DPA2" && version %in% c(3530,3520,3510,3500,3490,3480,3470,
                                                 3460,3450,3440,3430,3420,3410,3400,
                                                 3390,3380,3370,3360,3350,3340,3330,
                                                 3320,3310,3300,3290,3280,3270)) {
                                delete_lines <- c(1,2,3)} # 3.53.0 - 3.27.0
          
                              if(loci[i] == "DPB2" && version %in% c(3530,3520,3510,3500,3490,3480,3470,
                                                 3460,3450,3440,3430,3420,3410,3400,
                                                 3390,3380,3370,3360,3350,3340,3330,
                                                 3320,3310,3300,3290,3280,3270,3260,
                                                 3250,3240)) {
                                delete_lines <- c(1,2,3)} # 3.53.0 - 3.24.0
          
                              if(loci[i] %in% c("N","S","U","Y") && version %in% c(3350,3340,3330)) {
                                delete_lines <- c(1,2,3)} # 3.35.0 - 3.33.0
          
                              if(loci[i] =="Y" && version %in% c(3320,3310,3300,3320,3310,
                                                                 3300,3290,3280,3270,3260,
                                                                 3250,3240,3230,3220,3210,3200)) {
                               delete_lines <- c(1,2,3)} # 3.32.0 - 3.20.0
          
                              if(loci[i] %in% c("W","T") && version == 3270) {
                                delete_lines <- c(1,2,3)} # 3.35.0 - 3.33.0
  
                        } else{
                  delete_lines <- c(1,2,3)}
        divide <- 3
        sequence_name <- "cDNAsequence"
      } else if(source[j]=="gDNA"){
        suffix <- "_gen.txt"
        type <- "gDNA"
        delete_lines <- c(1,2)
        divide <- 1
        sequence_name <- "gDNAsequence"
      }
      
      if(version == 3131) {version <- 3130} ## Fix for downloading the alignments, as the repo URL includes "3130", but the files use "3131".
      
      #downloads relevant locus alignment file -- readLines allows for space preservation, which is important in
      #finding where the alignment sequence starts
      if(source[j] == "AA"|source[j] == "cDNA"){
        alignment[[loci[i]]] <- readLines(paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/", repoVersion(version), "/alignments/",paste(ifelse(loci[i]=="DRB1"|loci[i]=="DRB3"|loci[i]=="DRB4"|loci[i]=="DRB5","DRB",loci[i]),suffix,sep=""),sep=""),-1,ok=TRUE,skipNul = FALSE)
      }  else if(source[j] == "gDNA"){
        alignment[[loci[i]]] <- readLines(paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/", repoVersion(version), "/alignments/",paste(loci[i],suffix,sep=""),sep=""),-1,ok=TRUE,skipNul = FALSE)
      }

      #if version is equal to latest, or if version is numeric and > 3310, obtain
      #version number from line 3, and skip first 7 rows and last 3 rows
      if((version == "Latest") | (is.numeric(version) & version > 3310)){
        # alignmentVersion[[loci[i]]] <-alignment[[loci[i]]][3]
        alignmentVersion[[loci[i]]] <-substr(alignment[[loci[i]]][3],12,nchar(alignment[[loci[i]]][3]))

        #alters alignment file by cutting out non-pertinent information in beginning
        #and end of alignment file
        alignment[[loci[i]]] <- head(alignment[[loci[i]]],-3)
        alignment[[loci[i]]] <- tail(alignment[[loci[i]]],-7)
      }
          ### Fix for extraneous block of allele name rows in HLA-B cDNA alignment in 3.44.0 and 3.43.0.
      if(loci[i] == "B" && version == 3440 && source == "cDNA"){ alignment[[loci[i]]] <- alignment[[loci[i]]][-c(127540:135510)] } # 3.44.0
      if(loci[i] == "B" && version == 3430 && source == "cDNA"){ alignment[[loci[i]]] <- alignment[[loci[i]]][-c(124129:131886)] } # 3.43.0
         ### Fix for extraneous HLA-C allele name rows in version 3.2.0 
      if(loci[i] == "C" && version == 320 && source == "cDNA"){ alignment[[loci[i]]] <- alignment[[loci[i]]][-c(14467:15430)] } # 3.43.0
      
         ### Fix for missing "AA codon" lines in HFE versions 3.27.0 to 3.22.0.
      if(loci[i] == "HFE" && version %in% c(3270,3260,3250,3240,3230,3220) && source == "cDNA") {
        
              firstPos <- -25
              alignment[[loci[i]]] <- addCodonLine(alignment[[loci[i]]],firstPos)
        
            }
      
         ### Fix for missing "AA codon" lines in DPA, DPB, TAP1 and TAP2 in version 0.00.0
      if(version == 300 && source == "cDNA" && loci[i] %in% c("DPA","DPB","TAP1","TAP2")) {

                  if(loci[i] == "DPA") { firstPos <- -31 }
                  if(loci[i] == "DPB") { firstPos <- -29 }
                  if(loci[i] %in% c("TAP1","TAP2")) { firstPos <- 1 }
        
                  alignment[[loci[i]]] <- addCodonLine(alignment[[loci[i]]],firstPos)
                }
        
      ### Fix for missing carriage-return between lines 2 and 3 in HLA-V versions 3.14.0, and converting "3.15.0" to "3.14.0"    
       if(loci[i] == "V" && version == 3140) { 
                  alignment[[loci[i]]] <- append(alignment[[loci[i]]],"Sequences Aligned: 2014 January 17",after=2)
                  alignment[[loci[i]]][2] <- "IMGT/HLA Release: 3.14.0"
                  }
      
      #if version is numeric and <= 3310, obtain version number from line 2, and
      #skip first 6 rows and last 2 rows
      if((is.numeric(version) & version <= 3310)){ ## this format is 'IPD-IMGT/HLA Release: 3.31.0'. I'd like it to be 'IPD-IMGT/HLA 3.31.0' to match the post 3.31.0 version structure ****
        #alignmentVersion[[loci[i]]] <-alignment[[loci[i]]][2]
        alignmentVersion[[loci[i]]] <- paste(substr(alignment[[loci[i]]][2],1,12),substr(alignment[[loci[i]]][2],gregexpr(pattern =':',alignment[[loci[i]]][2])[[1]][1]+2,nchar(alignment[[loci[i]]][2])),sep = " ")

        #alters alignment file by cutting out non-pertinent information in beginning
        #and end of alignment file
        alignment[[loci[i]]] <- head(alignment[[loci[i]]],-2)
        alignment[[loci[i]]] <- tail(alignment[[loci[i]]],-6)
      }

      #see countSpaces function at beginning of script (AA table only)
      #Counts difference between Prot to -30 and beginning of Prot to -30 + 1 due to zero number indexing to find where
      #the alignment sequence actually starts
      if(source[j]=="AA"){space_diff[[loci[i]]]<-(nchar(strsplit(alignment[[loci[i]]][3], " ")[[1]][2])+countSpaces(alignment[[loci[i]]][3])[2]+1)-countSpaces(alignment[[loci[i]]][2])[1]}

      #reduces repeated whitespace in alignment file and removes rows with empty values for proper
      #start and stop subsetting
      alignment[[loci[i]]] <-str_squish(alignment[[loci[i]]])
      alignment[[loci[i]]] <-alignment[[loci[i]]][-which(alignment[[loci[i]]] == "")]

      #determines positions of "cDNA" or "Prot" and the end of that reference block segment
      start[[loci[i]]] <-as.numeric(grep(type, alignment[[loci[i]]]))
      
      ## Fix for DOA in 3.01.0 - Only one sequence block
      if(length(start[[loci[i]]]) == 1) {
        end[[loci[i]]] <- length(alignment[[loci[i]]]) ## the only end is the end of the sequence block
          } else {
      end[[loci[i]]] <- as.numeric(c(start[[loci[i]]][2:length(start[[loci[i]]])]-1,length(alignment[[loci[i]]]))) }
      
      if(version == 3480 && source[j] == "gDNA" && loci[i] == "DRB1"){ ## Fix for DRB1*15:200:01:01N and DRB1*15:200:01:02N in 3.48.0 gDNA alignment
          alignment[[loci[i]]] <- gsub("DRB1*15:200:01:01N","DRB1*15:200:01:01N ",alignment[[loci[i]]], fixed=TRUE)
          alignment[[loci[i]]] <- gsub("DRB1*15:200:01:02N","DRB1*15:200:01:02N ",alignment[[loci[i]]], fixed=TRUE)
      }

      if(source[j]=="AA"){
        #extracts rows with "Prot" and reference sequence position information
        #extracts only relevant reference sequence positions
        #NOTE: the first row containing "Prot" contains two numbers -- -30 and 1 -- where only -30, is extracted,
        #as the actual sequence start will always be 1
        for (k in 1:length(start[[loci[i]]])){

          prot_extractions[[loci[i]]][k]<-strsplit(alignment[[loci[i]]][start[[loci[i]]][k]], " ")

          refblock_number[[loci[i]]][k]<-as.numeric(sapply(prot_extractions[[loci[i]]][k], "[", 2))

          #determines the alignment start by adding -30 to the difference between white spaces found above
          alignment_start[[loci[i]]]<-refblock_number[[loci[i]]][1]+space_diff[[loci[i]]]
        }
      }  else if(source[j]=="cDNA"){
        #these loci do not have protein sequences; set alignment start to 1
        if(loci[i] %in% HLAgazeteer$nuc[!HLAgazeteer$nuc %in% HLAgazeteer$prot]){
          alignment_start[[loci[i]]] <- 1
        } else{
        #determines the alignment start by finding the second vector in second list and removing "codon"
        alignment_start[[loci[i]]] <-as.numeric(sub("AA codon ", "", alignment[[loci[i]]][2]))
        }
        #print(alignment_start[[loci[i]]])
        DNA_start[[loci[i]]] <-as.numeric(sub("cDNA ", "", alignment[[loci[i]]][1]))
        #print(DNA_start[[loci[i]]])
      }  else if(source[j]=="gDNA"){
        DNA_start[[loci[i]]] <-as.numeric(sub("gDNA ", "", alignment[[loci[i]]][1]))
      }

      #closes all white space in the alignment file, except for the white space separating the allele and peptide sequence
      alignment[[loci[i]]] <-paste(substr(alignment[[loci[i]]],1,regexpr(" ",text = alignment[[loci[i]]],fixed = TRUE)), gsub(" ","",substr(alignment[[loci[i]]],regexpr(" ",text = alignment[[loci[i]]],fixed = TRUE),nchar(alignment[[loci[i]]]))),sep = "")

      #string splits at white spaces to yield allele and peptide sequences
      alignment[[loci[i]]]  <- strsplit(alignment[[loci[i]]]," ", fixed=T)

      #binds the previously split strings by row
      alignment[[loci[i]]] <- do.call(rbind,alignment[[loci[i]]])

      #if the seq column is equal to the allele column due to premature peptide termination,
      #insert a blank in place of the allele in the seq column (AA table only)
      if(source[j]=="AA"){alignment[[loci[i]]][which(alignment[[loci[i]]][,1]==alignment[[loci[i]]][,2]),2] <- ""}

      #renames columns to "alleles" and "seq"
      colnames(alignment[[loci[i]]])<-c(paste(loci[i], "alleles", sep="_"), "seq")

      #due to ANHIG formatting, cases where an allele contains newly reference peptide sequences will not
      #contain the same number of rows as previous reference peptide blocks
      #this for loop is invoked to add "."for all other alleles for each character in the newly reference peptide
      #to preserve structural integrity
      #changes 10/9/19 to accommodate if there is more than one extraneous allele with an extended amino acid sequence

      if(source[j]=="cDNA"|source[j]=="AA"&loci[i]=="TAP2"){
        for(l in 1:length(start[[loci[i]]])){
          if(nrow(alignment[[loci[i]]][start[[loci[i]]][l]:end[[loci[i]]][l],])!=nrow(alignment[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],])){
            x<-as.data.frame(alignment[[loci[i]]][,1][start[[loci[i]]][1]:end[[loci[i]]][1]][-c(1,2)], stringsAsFactors = F)
            colnames(x)<-paste(loci[i], "alleles", sep="_")
            x<-cbind.data.frame(x, seq=as.character(paste(rep(".", nchar(tail(alignment[[loci[i]]][,2], 1))), collapse = "")), stringsAsFactors=FALSE)
            y<-data.frame(tail(alignment[[loci[i]]], (nrow(alignment[[loci[i]]][start[[loci[i]]][l]:end[[loci[i]]][l],][nrow(alignment[[loci[i]]][start[[loci[i]]][l]:end[[loci[i]]][l],])!=nrow(alignment[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],]),])-2)), stringsAsFactors = F)
            x$seq[match(y[,1], x[,1])]<-y$seq
            alignment[[loci[i]]]<-as.matrix(rbind(head(alignment[[loci[i]]], -(nrow(alignment[[loci[i]]][start[[loci[i]]][l]:end[[loci[i]]][l],][nrow(alignment[[loci[i]]][start[[loci[i]]][l]:end[[loci[i]]][l],])!=nrow(alignment[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],]),])-2)), x))
            start[[loci[i]]]<-as.numeric(grep(type, alignment[[loci[i]]]))
            
            end[[loci[i]]] <- as.numeric(c(start[[loci[i]]][2:length(start[[loci[i]]])]-1,nrow(alignment[[loci[i]]])))
            }
        }
      }

      if(source[j]=="AA"|source[j]=="gDNA"|source[j]!="AA"&loci[i]!="TAP2"){
        for(l in 1:length(start[[loci[i]]])){
          if(nrow(alignment[[loci[i]]][start[[loci[i]]][l]:end[[loci[i]]][l],])!=nrow(alignment[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],])){
            x<-as.data.frame(alignment[[loci[i]]][,1][start[[loci[i]]][1]:end[[loci[i]]][1]][-c(1,2)], stringsAsFactors = F)
            colnames(x)<-paste(loci[i], "alleles", sep="_")
            temp_vec<-alignment[[loci[i]]][start[[loci[i]]][l]:end[[loci[i]]][l],][-c(1,2),]

            #if there is only one allele with an extended amino acid sequence, the conversion from a named
            #vector to a data frame is not properly
            #if the character vector is = 2, there is only one allele with an extended amino acid sequence, so use
            #as.list and then data.frame so data frame is created correctly
            if(length(temp_vec) == 2){
              temp_filter<-data.frame(as.list(temp_vec), stringsAsFactors = F)
            } else{
              temp_filter<-data.frame(temp_vec, stringsAsFactors = F)
            }

            #find the greatest number of peptides in extraneous alleles to determine
            #how many "." to add on
            max_nchar<-(temp_filter %>%
                          add_column(nchar = nchar(.$seq)) %>%
                          filter(nchar == max(nchar)))$nchar
            x<-cbind.data.frame(x, seq=as.character(paste(rep(".", max_nchar[1]), collapse = "")), stringsAsFactors=FALSE)
            y<-data.frame(tail(alignment[[loci[i]]], (nrow(alignment[[loci[i]]][start[[loci[i]]][l]:end[[loci[i]]][l],][nrow(alignment[[loci[i]]][start[[loci[i]]][l]:end[[loci[i]]][l],])!=nrow(alignment[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],]),])-2)), stringsAsFactors = F)
            x$seq[match(y[,1], x[,1])]<-y$seq
            alignment[[loci[i]]]<-as.matrix(rbind(head(alignment[[loci[i]]], -(nrow(alignment[[loci[i]]][start[[loci[i]]][l]:end[[loci[i]]][l],][nrow(alignment[[loci[i]]][start[[loci[i]]][l]:end[[loci[i]]][l],])!=nrow(alignment[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],]),])-2)), x))
            start[[loci[i]]]<-as.numeric(grep("Prot", alignment[[loci[i]]]))
            end[[loci[i]]] <- as.numeric(c(start[[loci[i]]][2:length(start[[loci[i]]])]-1,nrow(alignment[[loci[i]]])))}
        }
      }

      #if a locus has extra formatting, resulting in unequal rows, start and end will be updated to reflect subsetting
      #if a locus has no extra formatting, start and end will remain the same, as procured by earlier code
      if(length(start[[loci[i]]]) > 1) { # Fix #2a for DOA in 3.01.0 -- only need to bind together multiple sequence blocks
      for(m in 1:length(start[[loci[i]]])){
        HLAalignments[[loci[i]]]<-cbind(HLAalignments[[loci[i]]], alignment[[loci[i]]][start[[loci[i]]][m]:end[[loci[i]]][m],])
              }
          } else { HLAalignments[[loci[i]]] <- alignment[[loci[i]]]}

      #removes rows containing "cDNA", "AA codon", or "Prot"
      HLAalignments[[loci[i]]] <- HLAalignments[[loci[i]]][-delete_lines,]

      #designates columns to be combined as every other so allele names are not included
      #in pasting all the amino acid sequences together 
      
      ### Fix #2b for DOA in 3.01.0 -- only need to bind together multiple sequence blocks 
     if(length(start[[loci[i]]]) > 1) {
            cols<-seq(0, ncol(HLAalignments[[loci[i]]]), by=2)
            HLAalignments[[loci[i]]]<-cbind(HLAalignments[[loci[i]]][,1], apply(HLAalignments[[loci[i]]][,cols], 1 ,paste, collapse = "")) 
            }                     

      #creates a new matrix with the number of columns equal to the number of characters in the reference sequence
      corr_table[[loci[i]]]<-matrix(, nrow = 3, ncol = as.numeric(nchar(HLAalignments[[loci[i]]][,2][1])))

      if(source[j]=="AA"|source[j]=="cDNA"){
        #if the first position enumeration is negative (i.e. has a leader peptide sequence), determines alignment length based on the total number of characters plus the alignment start (which is negative)
        if(grepl("-", alignment_start[[loci[i]]][[1]])==TRUE){
          alignment_length[[loci[i]]]<-(as.numeric(nchar(HLAalignments[[loci[i]]][,2][1]))/divide)+(alignment_start[[loci[i]]])
        } else{ #if there is no leader peptide (i.e sequence starts at 1), determines alignment length based on total number of characters
          alignment_length[[loci[i]]]<-as.numeric(nchar(HLAalignments[[loci[i]]][,2][1]))
        }

        #pastes alignment_start to alignment_length together in sequential order
        #captures output as "w"
        w[[loci[i]]] <- capture.output(cat(alignment_start[[loci[i]]]:alignment_length[[loci[i]]]))

        #splits string formed by cat for separate character variables
        alignment_positions[[loci[i]]]<-as.character(unlist(strsplit(w[[loci[i]]], " ")))

        #eliminates "0" if present in the alignment positions, as the alignment sequence from ANHIG does not contain 0
        if("0" %in% alignment_positions[[loci[i]]]==TRUE){
          alignment_positions[[loci[i]]]<-alignment_positions[[loci[i]]][-which(alignment_positions[[loci[i]]] == 0)]
        }
      }

      #triples alignment to account for codons (cDNA table only)
      if(source[j]=="cDNA"){
        for(n in 1:length(alignment_positions[[loci[i]]])){
          alignment_positionsx3[[loci[i]]]<-c(alignment_positionsx3[[loci[i]]],alignment_positions[[loci[i]]][n],alignment_positions[[loci[i]]][n],alignment_positions[[loci[i]]][n])
        }
        alignment_positions[[loci[i]]]<-alignment_positionsx3[[loci[i]]]
      }

      #string splits to extract locus in the allele name
      #assigns to new variable "aligned"
      aligned[[loci[i]]]<- as.matrix(do.call(rbind,strsplit(HLAalignments[[loci[i]]][,1],"[*]")))

      #adds a new column of pasted locus and trimmed two field alleles to aligned
      aligned[[loci[i]]]<- cbind(aligned[[loci[i]]], paste(aligned[[loci[i]]][,1], apply(aligned[[loci[i]]],MARGIN=c(1,2),FUN=getField,res=2)[,2], sep="*"))

      #binds aligned and HLAalignments -- renames columns
      HLAalignments[[loci[i]]] <- cbind(aligned[[loci[i]]], HLAalignments[[loci[i]]])
 
      colnames(HLAalignments[[loci[i]]]) <- c("locus", "full_allele", "trimmed_allele", "allele_name", sequence_name)

      if(source[j]=="AA"){
        #if locus is DRB3/4/5, use DRB line 1 as reference sequence
        if(loci[i]=="DRB3"|loci[i]=="DRB4"|loci[i]=="DRB5"){
          #sets refexon to a reference peptide for each HLA locus based on the reference sequences in HLAalignments
          refexon[[loci[i]]] <- rbind(HLAalignments[[loci[i]]][1,])[which(rbind(HLAalignments[[loci[i]]][1,])[,"locus"]=="DRB1"),sequence_name]
        } else{       
                ## Fix error in version 3.23.0 with DPA1 cDNA alignments, where the filename is "DPA_prot.txt, and DPB1 with filenane "DPB_prot.txt".
                if(loci[i] %in% c("DPA","DPB","DQA","DQB") && source == "AA" && version %in% c(3230,3220,3210,3200,3190,3180,
                                                                                               3170,3160,3150,3140,3130,3120,
                                                                                               3110,3100,390,380,370,360,
                                                                                               350,340,330,320,310,300)) {
                  
                 if(loci[i] == "DPA") {refexon[[loci[i]]] <- rbind(HLAalignments[[loci[i]]][1,])[which(rbind(HLAalignments[[loci[i]]][1,])[,"locus"]=="DPA1"),sequence_name]}
                 if(loci[i] == "DPB") {refexon[[loci[i]]] <- rbind(HLAalignments[[loci[i]]][1,])[which(rbind(HLAalignments[[loci[i]]][1,])[,"locus"]=="DPB1"),sequence_name]}
                 if(loci[i] == "DQA") {refexon[[loci[i]]] <- rbind(HLAalignments[[loci[i]]][1,])[which(rbind(HLAalignments[[loci[i]]][1,])[,"locus"]=="DQA1"),sequence_name]}
                 if(loci[i] == "DQB") {refexon[[loci[i]]] <- rbind(HLAalignments[[loci[i]]][1,])[which(rbind(HLAalignments[[loci[i]]][1,])[,"locus"]=="DQB1"),sequence_name]}
                  
                  } else {
                      refexon[[loci[i]]] <- rbind(HLAalignments[[loci[i]]][1,])[which(rbind(HLAalignments[[loci[i]]][1,])[,"locus"]==loci[i]),sequence_name]
                          }
                      }
        #if input loci is DRB, use grep statement to match input loci to loci in HLAalignments
        if(loci[i]=="DRB"){
          refexon[[loci[i]]]<-rbind(HLAalignments[[loci[i]]][1,])[grepl(loci[i], rbind(HLAalignments[[loci[i]]][1,])[,"locus"]), sequence_name]}
          }

      #splits sequence column at every amino acid or nucleotide, resulting in a split list of each for each row
      pepsplit[[loci[i]]] <- sapply(HLAalignments[[loci[i]]][,sequence_name],strsplit,split="*")

      #fills in spaces with '.' for alleles with premature termination to make it the same number of characters
      #as the reference sequence (AA table only)
      if(source[j]=="AA"){pepsplit[[loci[i]]]<- lapply(pepsplit[[loci[i]]],function(x) c(x,rep(".",nchar(refexon[[loci[i]]])-length(x))))}

      #nullifies variable names
      names(pepsplit[[loci[i]]]) <- NULL

      #binds pep_split together by element in its previous list form by row ### need to skip this for DPA2 and DPB2??
      pepsplit[[loci[i]]]<- do.call(rbind,pepsplit[[loci[i]]])

      #binds all columns together to form desired output, as described above
      HLAalignments[[loci[i]]] <- cbind.data.frame(HLAalignments[[loci[i]]][,1:4],pepsplit[[loci[i]]], stringsAsFactors=FALSE)

      if(source[j]=="AA"|source[j]=="cDNA"){
        #get reference sequence from DRB alignment
        if(loci[i]=="DRB3"|loci[i]=="DRB4"|loci[i]=="DRB5"){
          DRBref<-HLAalignments[[loci[i]]][1,]}

        #if the locus is DRB1/3/4/5, subset the specific locus from DRB alignment
        if(loci[i]=="DRB1"|loci[i]=="DRB3"|loci[i]=="DRB4"|loci[i]=="DRB5"){
          HLAalignments[[loci[i]]]<-HLAalignments[[loci[i]]][HLAalignments[[loci[i]]]$locus==loci[i],]}

        #add reference sequence to DRB3/4/5, reset row names to numerical order
        if(loci[i]=="DRB3"|loci[i]=="DRB4"|loci[i]=="DRB5"){
          HLAalignments[[loci[i]]]<- rbind(DRBref, HLAalignments[[loci[i]]])
          rownames(HLAalignments[[loci[i]]])<-NULL}
      }

      #inputs HLAalignments alignment sequence into the corr_table with "InDel" still present
      corr_table[[loci[i]]][1,]<-names(HLAalignments[[loci[i]]][5:ncol(HLAalignments[[loci[i]]])])

      #finds positions in HLAalignments that have ".", indicating an inDel
      inDels[[loci[i]]]<-colnames(HLAalignments[[loci[i]]][1, 5:ncol(HLAalignments[[loci[i]]])][HLAalignments[[loci[i]]][1, 5:ncol(HLAalignments[[loci[i]]])] %in% "."])

      #finds positions in HLAalignments that have "|", indicating a exon boundary
      exonB[[loci[i]]]<-colnames(HLAalignments[[loci[i]]][1, 5:ncol(HLAalignments[[loci[i]]])][HLAalignments[[loci[i]]][1, 5:ncol(HLAalignments[[loci[i]]])] %in% "|"])

      #inDel inclusion if there are inDels present
      if(length(inDels[[loci[i]]])!=0){
        for(o in 1:length(inDels[[loci[i]]])){
          corr_table[[loci[i]]][2,][inDels[[loci[i]]][[o]]==corr_table[[loci[i]]][1,]]<-paste("INDEL", o, sep="-")
          corr_table[[loci[i]]][3,][inDels[[loci[i]]][[o]]==corr_table[[loci[i]]][1,]]<-paste("INDEL", o, sep="-")
        }
      }

      #exonB inclusion if there are exon boundaries present
      if(length(exonB[[loci[i]]])!=0){
        for(p in 1:length(exonB[[loci[i]]])){
          corr_table[[loci[i]]][2,][exonB[[loci[i]]][[p]]==corr_table[[loci[i]]][1,]]<-paste("EXONB", p, sep="-")
          corr_table[[loci[i]]][3,][exonB[[loci[i]]][[p]]==corr_table[[loci[i]]][1,]]<-paste("EXONB", p, sep="-")
        }
      }

      if(source[j]=="AA"|source[j]=="cDNA"){
        #pastes alignment_positions into corr_table accounting for InDels and exon boundaries
        codon_count<-as.vector(1)
        for(q in 1:length(corr_table[[loci[i]]][2,])){
          if(is.na(corr_table[[loci[i]]][2,q])){
            corr_table[[loci[i]]][2,q]<-alignment_positions[[loci[i]]][codon_count]
            codon_count<-codon_count+1
          }
        }
      }

      # create number sequence for alignment (cDNA and gDNA tables only)
      if(source[j]=="cDNA"){ alignment_positions[[loci[i]]]<-as.character(unlist(strsplit(capture.output(cat(DNA_start[[loci[i]]]:length(corr_table[[loci[i]]][3,]))), " ")))
      } else if(source[j]=="gDNA"){ if(DNA_start[[loci[i]]] < 0) {alignment_positions[[loci[i]]]<-as.character(unlist(strsplit(capture.output(cat(c(DNA_start[[loci[i]]]:-1,1:length(corr_table[[loci[i]]][3,])))), " "))) ## most loci have negative start positions
      } else {alignment_positions[[loci[i]]]<-as.character(unlist(strsplit(capture.output(cat(c(DNA_start[[loci[i]]]:length(corr_table[[loci[i]]][3,])))), " ")))  #### HLA-P starts at position 1, but in theory some other locus could start at a different positive position
      }
      } # In 3.53.0, only the HLA-P alignment started with 1, so, using the original version above, the first four alignment_positions values are "1"   "0"   "-1"  "1"
      # pastes number sequence into corr_table accounting for InDels and exon boundaries (cDNA and gDNA tables only)

      #pastes number sequence into corr_table accounting for InDels and exon boundaries (cDNA and gDNA tables only)
      if(source[j]=="cDNA"|source[j]=="gDNA"){
        cDNAcount<-as.vector(1)
        for(r in 1:length(corr_table[[loci[i]]][3,])){
          if(is.na(corr_table[[loci[i]]][3,r])){
            corr_table[[loci[i]]][3,r]<-alignment_positions[[loci[i]]][cDNAcount]
            cDNAcount<-cDNAcount+1
          }
        }
      }

      #run if InDels or exon boundaries are present in alignment
      if(any(grepl("INDEL|EXONB", corr_table[[loci[i]]][2,]))){
        #finds which positions have InDels or exon boundaries
        v<-as.numeric(corr_table[[loci[i]]][1,][which(grepl("INDEL", corr_table[[loci[i]]][2,]))])
        w<-as.numeric(corr_table[[loci[i]]][1,][which(grepl("EXONB", corr_table[[loci[i]]][2,]))])

        #splits cumulative sequences
        vsplit<-split(v, cumsum(c(TRUE, diff(v) != 1L)))
        wsplit<-split(w, cumsum(c(TRUE, diff(w) != 1L)))


        if(any(grepl("INDEL", corr_table[[loci[i]]][2,]))){
          #adds increments of .1 to any InDels in row 2
          if(source[j]=="AA"|source[j]=="cDNA"){
            for(s in 1:length(vsplit)){
              if(grepl("-", corr_table[[loci[i]]][2,][vsplit[[s]][[1]]-1])){
                corr_table[[loci[i]]][2,][vsplit[[s]]]<-paste(corr_table[[loci[i]]][2,][vsplit[[s]][[1]]-1], gsub(0, "", seq(1, length(vsplit[[s]]))/10), sep="")
             # }
     #         else{corr_table[[loci[i]]][2,][vsplit[[s]]]<-as.numeric(corr_table[[loci[i]]][2,][vsplit[[s]][[1]]-1])+as.numeric(paste(0,".",seq(1, length(vsplit[[s]])),sep = ""))}
             }  else{corr_table[[loci[i]]][2,][vsplit[[s]]]<-paste(corr_table[[loci[i]]][2,][vsplit[[s]][[1]]-1], paste(".",seq(1,length(vsplit[[s]])),sep=""),sep="")}
            }
          }


          #adds increments of .1 to any InDels in row 3
          if(source[j]=="cDNA"|source[j]=="gDNA"){
            for(t in 1:length(vsplit)){
              corr_table[[loci[i]]][3,][vsplit[[t]]]<-paste(corr_table[[loci[i]]][3,][vsplit[[t]][[1]]-1],".",seq(1, length(vsplit[[t]])),sep = "")
            }
          }
        }


        if(any(grepl("EXONB", corr_table[[loci[i]]][2,]))){

          #changes EXONB to E.(# of last exon)-(# of next exon)
          if(source[j]=="cDNA"){
            for(u in 1:length(wsplit)){
              corr_table[[loci[i]]][2,][wsplit[[u]]]<-paste("E.",u,"-",u+1, sep = "")
              corr_table[[loci[i]]][3,][wsplit[[u]]]<-paste("E.",u,"-",u+1, sep = "")
            }
          }

          #changes EXONB to corresponding boundary e.g. UTR-E.1
          #also creates table (atlas) of gDNA boundaries and their locus
          if(source[j]=="gDNA"){
            if(loci[i] %in% names(fragmentFeatureNames)) { # simplified gDNA atlas builder for genes in the fragemtnFeatureNames object

              gDNA_atlas[[loci[i]]]<-matrix(, nrow = 1, ncol = length(wsplit))
              gDNA_atlas[[loci[i]]]<- as.numeric(corr_table[[i]][3,w+1])
              newNames <- rep(NA,length(fragmentFeatureNames[names(fragmentFeatureNames) %in% loci[i] == TRUE][[1]][[1]])-1)
              for(f in 1:(length(newNames))) {
                newNames[f] <- paste(fragmentFeatureNames[names(fragmentFeatureNames) %in% loci[i] == TRUE][[1]][[1]][f],"-",fragmentFeatureNames[names(fragmentFeatureNames) %in% loci[i] == TRUE][[1]][[1]][f+1], sep="")
              }
              gDNA_atlas[[loci[i]]] <- as.data.frame(t(gDNA_atlas[[loci[i]]]))
              colnames(gDNA_atlas[[loci[i]]]) <- newNames
              rownames(gDNA_atlas[[loci[i]]]) <- "gDNA"

            } else{
              intron<-as.vector(paste("I.",seq(1,((length(wsplit)-1)/2)-0.5),sep = ""))
              exon<-as.vector(paste("E.",seq(1,((length(wsplit)-1)/2)+0.5),sep = ""))
              a<-unlist(strsplit(paste0(exon," ",intron), " "))
              a<-c("UTR",head(a, length(a)-1),"UTR")
              b<-vector()
              for(y in 1:(length(a)-1)){
                b<-c(b,paste(a[[y]],"-",a[[y+1]],sep = ""))}
              gDNA_atlas[[loci[i]]]<-matrix(, nrow = 1, ncol = length(wsplit))
              #covers 3 loci
              if(length(b) != length(gDNA_atlas[[loci[i]]])) {
                gDNA_atlas[[loci[i]]]<-matrix(, nrow = 1, ncol = length(b))
              }
              colnames(gDNA_atlas[[loci[i]]])<-b
              rownames(gDNA_atlas[[loci[i]]])<-"gDNA"
              for(z in 1:length(wsplit)){
                gDNA_atlas[[loci[i]]][1,z]<-corr_table[[loci[i]]][3,][wsplit[[z]]+1]
                corr_table[[loci[i]]][3,][wsplit[[z]]]<-b[[z]]
              }
            }
          }
          #creates table (atlas) of cDNA boundaries and their locus
          if(source[j]=="cDNA"){
            cDNA_atlas[[loci[i]]]<-matrix(, nrow = 2, ncol = length(wsplit))
            colnames(cDNA_atlas[[loci[i]]])<-corr_table[[loci[i]]][2,][unlist(wsplit)]
            rownames(cDNA_atlas[[loci[i]]])<-c("cDNA","AA")
            for(z in 1:length(wsplit)){
              cDNA_atlas[[loci[i]]][1,z]<-corr_table[[loci[i]]][3,][wsplit[[z]]+1]
              cDNA_atlas[[loci[i]]][2,z]<-corr_table[[loci[i]]][2,][wsplit[[z]]+1]
            }
          }

          #creates table (atlas) of exon boundaries and their locus -- PROVISIONAL SJM
          if(source[j]=="AA"){
            AA_atlas[[loci[i]]]<-matrix(, nrow = 2, ncol = length(wsplit))
            colnames(AA_atlas[[loci[i]]])<-corr_table[[loci[i]]][2,][unlist(wsplit)]
            rownames(AA_atlas[[loci[i]]])<-c("cDNA","AA")
            for(z in 1:length(wsplit)){
              AA_atlas[[loci[i]]][1,z]<-corr_table[[loci[i]]][3,][wsplit[[z]]+1]
              AA_atlas[[loci[i]]][2,z]<-corr_table[[loci[i]]][2,][wsplit[[z]]+1]
            }
          }   ## End PROVISIONAL

        }

      }

      #Creates separate tables for AA/codon and DNA positions
      AA_codon_alignments[[loci[i]]]<-HLAalignments[[loci[i]]]
      DNAalignments[[loci[i]]]<-HLAalignments[[loci[i]]]

      if(source[j]=="cDNA"|source[j]=="AA"){
        #renames columns in AA_codon_alignments to codon names
        colnames(AA_codon_alignments[[loci[i]]]) <- c("locus","allele","trimmed_allele","allele_name", corr_table[[loci[i]]][2,])
        #distributes  reference sequence from row 1
        #into all other rows, if they contain a "-"
        #amino acids with changes will not be impacted
        for(w in 5:ncol(AA_codon_alignments[[loci[i]]])) {
          AA_codon_alignments[[loci[i]]][,w][which(AA_codon_alignments[[loci[i]]][,w]=="-")] <- AA_codon_alignments[[loci[i]]][,w][1]}
      }

      if(source[j]=="cDNA"|source[j]=="gDNA"){
        #renames columns in DNAalignments to cDNA names
        colnames(DNAalignments[[loci[i]]]) <- c("locus","allele","trimmed_allele","allele_name", corr_table[[loci[i]]][3,])
        #distributes  reference sequence from row 1
        #into all other rows, if they contain a "-"
        #amino acids with changes will not be impacted
        for(x in 5:ncol(DNAalignments[[loci[i]]])) {
          DNAalignments[[loci[i]]][,x][which(DNAalignments[[loci[i]]][,x]=="-")] <- DNAalignments[[loci[i]]][,x][1]}
      }

      if(source[j]=="AA"){
        final_alignment[[loci[i]]]<-c(final_alignment[[loci[i]]],AA = list(AA_codon_alignments[[loci[i]]]))
      }else if(source[j]=="cDNA"){
        comb_list[[loci[i]]]<-list(codon = AA_codon_alignments[[loci[i]]],cDNA = DNAalignments[[loci[i]]])
        final_alignment[[loci[i]]]<-c(final_alignment[[loci[i]]],comb_list[[loci[i]]])
      }else if(source[j]=="gDNA"){
        final_alignment[[loci[i]]]<-c(final_alignment[[loci[i]]],gDNA = list(DNAalignments[[loci[i]]]))
      }

      #Adds version number
      #     final_alignment[[loci[i]]]<-c(final_alignment[[loci[i]]],`ANHIG/IMGTHLA Alignments Version` = alignmentVersion[[loci[i]]])
      final_alignment[[loci[i]]]<-c(final_alignment[[loci[i]]],'Version' = alignmentVersion[[loci[i]]])

    }

  }
  # final_alignment<-c(final_alignment, atlas = list(atlas)) ## No need to add atlases as there is a separate atlas maker
  return(final_alignment)
}
