#atlasMaker v1.6.0 25 July 2024

################
##atlasMaker
#'Identify the Gene-Feature Boundaries in HLA Region Genes
#'
#'Builds an 'atlas' of the gene-feature (exon, intron and UTR) boundaries for IPD-IMGT/HLA loci, using IPD-IMGT/HLA Database alignments.
#'
#'@param loci A character string identifying a specific HLA gene (ex. "DRB1"). More than one gene can be specified (e.g., c("DRB1","DQB1")).
#'@param source A character string identifying the type of alignment used to build the atlas. The allowed options are "AA", for protein alignments, "cDNA", for nucleotide alignments, and "gDNA", for genomic alignments. More than one alignment type can be specified (e.g., c("AA","gDNA")).
#'@param version A character string identifying the release version of the IPD-IMGT/HLA Github repository for which alignments should be generated.
#'
#'@return A list object of 'atlas' dataframes and a 'release version' character string for each locus and alignment type.
#'
#'@description
#'The 'AA' atlas identifies the amino acid (AA) positions encoded by codons that span each exon (E) boundary. The 'cDNA' atlas describes the 3' nucleotide position in the cDNA sequence that follows each exon boundary, as well as the codon that spans each exon boundary. The 'gDNA' atlas describes first 3' nucleotide position following each UTR (U), exon or intron (I) boundary. Each feature is followed by its rank (e.g U3 is the 3' UTR). Non-standard gene-features (H, J, N, and S) are detailed in the documentation for the ffN() function. 
#'
#'@note For internal HLAtools use.
#'@note Nucleotide atlases for pseudogenes will include a 'codon' row populated with NA values.
#'@note Some HLA-DQB1*05 and *06 alleles include a 5th exon that is not present in the DQB1*05:01:01:01 reference allele. In this case, the Exon 4 to Exon 5 boundary is defined as 'Ins' (insertion). For all other alleles there is no E.4-5 boundary.
#'@note Boundaries for non-standard hybrid (H), join (J), novel (N) and segment (S) features may be included in gene fragment and pseudogene atlases.
#'
#'@importFrom DescTools StrLeft IsOdd
#'@importFrom dplyr filter
#'@importFrom stringr str_replace str_squish
#'@importFrom tibble add_column
#'@importFrom utils capture.output head tail
#'
#'@export
#'
atlasMaker<-function(loci, source, version = "Latest"){

  if(version != "Latest"){ #
    if(!validateVersion(version)){stop(paste(version," is not a valid IPD-IMGT/HLA Database release version."))}
      }else{ version <- getLatestVersion()}
  
  source <- (checkSource(source))
  
  loci <- multiLocusValidation(loci, source) # Added as a check to make sure that HLAgazeteer#version is enforced; otherwise validateLocus could return FALSE
  
  if(validateLocus(loci=loci,source=source)) { ## primary 'check'.

  subSource<- expressed<-genestructure<-nucVersion<-pep_start<-nuc<-nuc_df<-extract_ref<-nuc_extract<-start<-end<-pipe_split<-column_names<-boundaries<-boundary_split<-AA_atlas<-sapply(loci, function(x) NULL)

  #creates empty variables for future for loops
  comb_list<-start<-end<-alignment<-list()

  #empty variables for correspondence table
  alignmentVersion<-pepsplit<-refexon<-aligned<-final_alignment<-AA_codon_alignments<-DNAalignments<-HLAalignments<-exonB<-inDels<-gDNA_atlas<-cDNA_atlas<-atlas<-corr_table<-cols<-downloaded_segments<-w<-alignment_positionsx3<-alignment_positions<-alignment_length<-DNA_start<-alignment_start<-space_diff<-prot_extractions<-refblock_number<-sapply(loci, function(x) NULL)

  #sub out periods in version, if there are any
  version <- gsub(".", "", version, fixed=T)

  #if version is not latest, turn into numeric object
 # if(version != "Latest"){

    version <- as.numeric(version)

 # }

  for(i in 1:length(loci)){ # main loop of i values

    expressed[[i]] <- !loci[i] %in% HLAgazeteer$nuc[HLAgazeteer$nuc %in% c(HLAgazeteer$pseudo,HLAgazeteer$frag)] # Exclude pseudogenes and gene fragments

    for(j in 1:length(source)){ #### Source Loop of j values

      subSource[[j]] <- source[[j]] ## identify the submitted source so that that the cDNA code can be used for both AA and cDNA
      if(subSource[[j]] == "AA") { source[[j]] <- "cDNA" ## Bypass all of the AA code for source = AA, in favor of the cDNA code
        } else {subSource[[j]] <- 0}

      if(source[j]=="AA"){ ## start of the AA routine
        # download nucleotide alignment from IPD-IMGT/HLA Github repository
        nuc[[loci[[i]]]]<-readLines(paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/", repoVersion(version),"/alignments/",paste(ifelse(loci[[i]] %in% c("DRB1", "DRB3", "DRB4", "DRB5"),"DRB",loci[[i]]),"_nuc.txt",sep=""),sep=""),-1,ok=TRUE,skipNul = FALSE)

        if((version == "Latest") | (is.numeric(version) & version > 3310)){ ## Account for the changes in file header structures

          nucVersion[[loci[[i]]]]<-substr(nuc[[loci[i]]][3],12,nchar(nuc[[loci[i]]][3])) } else {
          nucVersion[[loci[[i]]]] <- paste(substr(nuc[[loci[i]]][2],1,12),substr(nuc[[loci[i]]][2],gregexpr(pattern =':',nuc[[loci[i]]][2])[[1]][1]+2,nchar(nuc[[loci[i]]][2])),sep = " ")
        } # closes two lines above

        # reduces repeated whitespace between allele and nucleotide
        nuc[[loci[[i]]]]<-str_squish(nuc[[loci[[i]]]])

        # remove impertinent header/footer information
        nuc[[loci[[i]]]] <- head(nuc[[loci[[i]]]],-2)
        nuc[[loci[[i]]]] <- tail(nuc[[loci[[i]]]],-6)

        # remove all whitespace, except for the whitespace between the allele and nucleotide sequence
        nuc[[loci[[i]]]]<-paste(substr(nuc[[loci[i]]],1,regexpr(" ",text = nuc[[loci[i]]],fixed = TRUE)), gsub(" ","",substr(nuc[[loci[i]]],regexpr(" ",text = nuc[[loci[i]]],fixed = TRUE),nchar(nuc[[loci[i]]]))),sep = "")

        # split at white spaces to yield allele and nucleotide sequences
        nuc[[loci[i]]]  <- strsplit(nuc[[loci[i]]]," ", fixed=T)

        # bind the previously split strings by row
        nuc[[loci[i]]] <- do.call(rbind,nuc[[loci[i]]])

        #extracts beginning alignment enumeration
        pep_start[[loci[[i]]]]<-as.numeric(gsub("codon", "", nuc[[loci[[i]]]][[2,2]]))
        ## correction for nuc alignments that start at codon 1 -- currently TAP1 and TAP2

        colnames(nuc[[loci[[i]]]])<-c(paste(loci[[i]], "alleles", sep="_"), "pepseq")

        # find start and end of each "cDNA" block
        start[[loci[i]]]<-as.numeric(grep("cDNA", nuc[[loci[i]]]))
        end[[loci[i]]] <- as.numeric(c(start[[loci[i]]][2:length(start[[loci[i]]])]-1,nrow(nuc[[loci[i]]])))

        # Due to ANHIG formatting, cases where an allele contains newly reference nucleotide sequences will not
        # contain the same number of rows as previous reference peptide blocks
        # this for loop is invoked to add "."for all other alleles for each character in the new reference nucleotide sequence
        # to preserve structural integrity
        for(k in 1:length(start[[loci[i]]])){ ## pad out ... loop
          if(nrow(nuc[[loci[[i]]]][start[[loci[i]]][k]:end[[loci[i]]][k],])!=nrow(nuc[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],])){
            x<-as.data.frame(nuc[[loci[i]]][,1][start[[loci[i]]][1]:end[[loci[i]]][1]][-c(1,2)], stringsAsFactors = F)
            colnames(x)<-paste(loci[[i]], "alleles", sep="_")
            x<-cbind.data.frame(x, pepseq=as.character(paste(rep(".", nchar(tail(nuc[[loci[i]]][,2], 1))), collapse = "")), stringsAsFactors=FALSE)
            y<-data.frame(tail(nuc[[loci[i]]], (nrow(nuc[[i]][start[[loci[i]]][k]:end[[loci[i]]][k],][nrow(nuc[[i]][start[[loci[i]]][k]:end[[loci[i]]][k],])!=nrow(nuc[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],]),])-2)), stringsAsFactors = F)
            colnames(y)<-c(paste(loci[[i]], "alleles", sep="_"), "pepseq")
            x$pepseq[match(y[,1], x[,1])]<-y$pepseq
            nuc[[loci[i]]]<-as.matrix(rbind(head(nuc[[loci[i]]], -(nrow(nuc[[i]][start[[loci[i]]][k]:end[[loci[i]]][k],][nrow(nuc[[i]][start[[loci[i]]][k]:end[[loci[i]]][k],])!=nrow(nuc[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],]),])-2)), x))
            start[[loci[i]]]<-as.numeric(grep("cDNA", nuc[[loci[i]]]))
            end[[loci[i]]] <- as.numeric(c(start[[loci[i]]][2:length(start[[loci[i]]])]-1,nrow(nuc[[loci[i]]])))}
        } # close ... loop

        # subset nucleotide alignment based on start and end blocks
        for(e in 1:length(start[[loci[[i]]]])){
          nuc_extract[[loci[i]]]<-cbind(nuc_extract[[loci[i]]], nuc[[loci[i]]][start[[loci[i]]][e]:end[[loci[i]]][e],])}

        # remove first two rows containing AA position and "Prot"
        nuc_extract[[loci[i]]] <- nuc_extract[[loci[i]]][-c(1,2,3),]

        # extract only the nucleotide reference sequence
        extract_ref[[loci[i]]]<-rbind(nuc_extract[[loci[[i]]]][1,])

        # paste every two columns together to get full nucleotide reference sequence without allele repetition
        cols<-seq(0, ncol(extract_ref[[loci[i]]]), by=2)
        extract_ref[[loci[i]]]<-cbind(extract_ref[[loci[i]]][,1], apply(extract_ref[[loci[i]]][,cols, drop=F], 1 ,paste, collapse = ""))

        # split reference sequences at pipes
        pipe_split[[loci[[i]]]]<-strsplit(extract_ref[[loci[[i]]]][,2], "|", fixed = T) ### Re: TAP1/TAP2 boundary error -- this pipe_spilt is correct

        for(z in 1:(length(pipe_split[[loci[[i]]]][[1]])-1)){ ## z-loop
          if(StrLeft(pipe_split[[loci[[i]]]][[1]][z], 2)!="..") {
            pipe_split[[loci[[i]]]][[1]][z]<-paste(pipe_split[[loci[[i]]]][[1]][z], StrLeft(pipe_split[[loci[[i]]]][[1]][z+1], 2), sep="")
          }

          if(StrLeft(pipe_split[[loci[[i]]]][[1]][z], 2)=="..") {
            pipe_split[[loci[[i]]]][[1]][z]<-pipe_split[[loci[[i]]]][[1]][z]
          }

          if(StrLeft(pipe_split[[loci[[i]]]][[1]][z+1], 2)=="..") {next}
          pipe_split[[loci[[i]]]][[1]][z+1]<-str_replace(pipe_split[[loci[[i]]]][[1]][z+1], substr(pipe_split[[loci[[i]]]][[1]][z+1],1, 2), "")
        } ## close z-loop
        # split every nucleotide within each boundary, and remove InDels
        for(k in 1:length(pipe_split[[loci[[i]]]][[1]])){ ## k-loop
          boundary_split[[loci[[i]]]][[k]]<-strsplit(pipe_split[[loci[[i]]]][[1]][[k]], "*")[[1]]
          if(boundary_split[[loci[[i]]]][[k]][[length(boundary_split[[loci[[i]]]][[k]])]]=="."){next}
          if(boundary_split[[loci[[i]]]][[k]][1]=="."){next}
          boundary_split[[loci[[i]]]][[k]]<-boundary_split[[loci[[i]]]][[k]][boundary_split[[loci[[i]]]][[k]]!="."]
        } ## close k-loop

        # form a blank dataframe with the number of rows equal to the number of exon boundaries
        AA_atlas[[loci[[i]]]]<-data.frame(matrix("", ncol=(length(boundary_split[[loci[[i]]]])-1)), stringsAsFactors = F)

        rownames(AA_atlas[[loci[[i]]]])<-"AA"

        # paste together individual nucleotides to form peptides; count number of peptides (codons) present between boundaries
        for(q in 1:(length(boundary_split[[loci[[i]]]])-1)){
          d <- seq.int(1L,length(boundary_split[[loci[[i]]]][[q]]),by = 3L)
          boundaries[[loci[[i]]]][[q]]<-length(paste0(boundary_split[[loci[[i]]]][[q]][d],boundary_split[[loci[[i]]]][[q]][d+1], boundary_split[[loci[[i]]]][[q]][d+2]))
        }

        # break out of for loop to add alignment start enumeration to actual start
        boundaries[[loci[[i]]]][[1]]<-boundaries[[loci[[i]]]][[1]]+pep_start[[loci[[i]]]] ### This seems to be the source of the discrepancy with TAP1 and TAP2 their pep_start values should be 0

        positions<-sapply(length(boundary_split[[loci[[i]]]])-1, function(x) NULL)

        tapLengths <- c(-1,-2,-2,-2,-2,-2,-2,-2,-3,-3) ## an ad hoc solution for the TAP1 and TAP2 genes. I'd prefer something programmatic, but I'm assuming that no more genes will be added
        # fill in AA_atlas information
        # boundaries obtained by adding up cumulative lengths
        for(q in 1:(length(boundary_split[[loci[[i]]]])-1)){
          column_names<-c(if(q==1){}else{column_names},column_names[[q]]<-paste("E", paste(seq(1, length(boundary_split[[loci[[i]]]]))[[q]], seq(1, length(boundary_split[[loci[[i]]]]))[[q+1]], sep="-"), sep="."))
          positions[[q]]<-as.numeric((ifelse(pep_start[[loci[[i]]]] == 1, tapLengths[q],0) + cumsum(boundaries[[loci[[i]]]])[[q]]))
                  AA_atlas[[loci[[i]]]][1,q]<-positions[[q]]
        }
        colnames(AA_atlas[[loci[[i]]]])<-column_names

        atlas[[loci[i]]] <-c(atlas[[loci[i]]], AA = list(AA_atlas[[loci[i]]]),`Version` = nucVersion[[loci[[1]]]])
        ## end of the AA routine
        } else if(source[j]=="cDNA" | source[j]=="gDNA"){ ### cDNA/gDNA routine
        HLAalignments<-sapply(loci, function(x) NULL)
        # expressions unique to each table source (AA or cDNA)
        if(source[j]=="AA"){
          suffix <- "_prot.txt"
          type <- "Prot"
          delete_lines <- c(1,2)
          divide <- 1
          sequence_name <- "AAsequence"
        } else if(source[j]=="cDNA"){
          suffix <- "_nuc.txt"
          type <- "cDNA"
          delete_lines <- c(1,2,3) # cDNA alignments for unexpressed genes (e.g., HLA-H, J, K, L, N, S, T, U, V, W and Y) only have 2 lines
          if(!expressed[[i]]) {delete_lines <- c(1,2)} ## prior to release 3.54.0, DPA2 and DPB2 had 3 lines, but are pseudogenes
          
          ### Some cDNAs for pseudogenes still have 3 lines. These lines fix that. 
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
          
          divide <- 3
          sequence_name <- "cDNAsequence"
        } else { if(source[j]=="gDNA"){
          suffix <- "_gen.txt"
          type <- "gDNA"
          delete_lines <- c(1,2)
          divide <- 1
          sequence_name <- "gDNAsequence"
        } }

        if(version == 3131) {version <- 3130} ## Fix for downloading the alignments, as the repo URL includes "3130", but the files use "3131".
        
        # download relevant locus alignment file -- readLines allows for space preservation, which is important in
        # finding where the alignment sequence starts
        if(source[j] == "AA"|source[j] == "cDNA"){
          alignment[[loci[i]]] <- readLines(paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/", repoVersion(version), "/alignments/",paste(ifelse(loci[[i]]=="DRB1"|loci[[i]]=="DRB3"|loci[[i]]=="DRB4"|loci[[i]]=="DRB5","DRB",loci[[i]]),suffix,sep=""),sep=""),-1,ok=TRUE,skipNul = FALSE)
        } else { if(source[j] == "gDNA"){
          alignment[[loci[i]]] <- readLines(paste("https://raw.githubusercontent.com/ANHIG/IMGTHLA/", repoVersion(version), "/alignments/",paste(loci[[i]],suffix,sep=""),sep=""),-1,ok=TRUE,skipNul = FALSE)
        } }

        ### Fix for extraneous block of allele name rows in HLA-B cDNA alignment in 3.44.0 and 3.43.0.
        if(loci[i] == "B" && version == 3440 && source == "cDNA"){ alignment[[loci[i]]] <- alignment[[loci[i]]][-c(127540:135510)] } # 3.44.0
        if(loci[i] == "B" && version == 3430 && source == "cDNA"){ alignment[[loci[i]]] <- alignment[[loci[i]]][-c(124129:131886)] } # 3.43.0
        ### Fix for extraneous HLA-C allele name rows in version 3.2.0 
        if(loci[i] == "C" && repoVersion(version) == 320 && source == "cDNA"){ alignment[[loci[i]]] <- alignment[[loci[i]]][-c(14467:15430)] } # 3.43.0
        
        ### Fix for missing "AA codon" lines in HFE versions 3.27.0 to 3.22.0.
        if(loci[i] == "HFE" && version %in% c(3270,3260,3250,3240,3230,3220) && source == "cDNA") {
          
          firstPos <- -25
          alignment[[loci[i]]] <- addCodonLine(alignment[[loci[i]]],firstPos)
          
        }
        
        ### Fix for missing "AA codon" lines in DPA, DPB, TAP1 and TAP2 in version 0.00.0
        if(repoVersion(version) == 300 && source == "cDNA" && loci[i] %in% c("DPA","DPB","TAP1","TAP2")) {
          
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
        
        
        # if version is equal to latest, or if version is numeric and > 3310, obtain
        # version number from line 3, and skip first 7 rows and last 3 rows
        if((version == "Latest") | (is.numeric(version) & version > 3310)){
          alignmentVersion[[loci[i]]] <-substr(alignment[[loci[i]]][3],12,nchar(alignment[[loci[i]]][3]))

          # alter alignment file by cutting out non-pertinent information in beginning
          # and end of alignment file
          alignment[[loci[i]]] <- head(alignment[[loci[i]]],-3)
          alignment[[loci[i]]] <- tail(alignment[[loci[i]]],-7)
        }

        # if version is numeric and <= 3310, obtain version number from line 2, and
        # skip first 6 rows and last 2 rows
        if((is.numeric(version) & version <= 3310)){
           alignmentVersion[[loci[i]]] <- paste(substr(alignment[[loci[i]]][2],1,12),substr(alignment[[loci[i]]][2],gregexpr(pattern =':',alignment[[loci[i]]][2])[[1]][1]+2,nchar(alignment[[loci[i]]][2])),sep = " ")

          # alters alignment file by cutting out non-pertinent information in beginning
          # and end of alignment file
          alignment[[loci[i]]] <- head(alignment[[loci[i]]],-2)
          alignment[[loci[i]]] <- tail(alignment[[loci[i]]],-6)
        }

        if(loci[[i]] %in%  c("DPA2","DPB2","Y")) { # In release version 3.53.0 and earlier these pseudogenes had 'AA codon' positions in their cDNA alignments, but in 3.54.0 those rows have been removed; similarly the HLA-Y pseudogene had an AA codon line prior to version 3.36.0 
          alignment[[loci[i]]] <- alignment[[loci[[i]]]][!substr(alignment[[loci[[i]]]][],1,9) == " AA codon"]
        }

        if(version == 3480 && source[j] == "gDNA" && loci[i] == "DRB1"){ ## Fix for DRB1*15:200:01:01N and DRB1*15:200:01:02N in 3.48.0 gDNA alignment
          alignment[[loci[i]]] <- gsub("DRB1*15:200:01:01N","DRB1*15:200:01:01N ",alignment[[loci[i]]], fixed=TRUE)
          alignment[[loci[i]]] <- gsub("DRB1*15:200:01:02N","DRB1*15:200:01:02N ",alignment[[loci[i]]], fixed=TRUE)
        }
        
        
        # See countSpaces function (AA table only)
        # Counts difference between Prot to -30 and beginning of Prot to -30 + 1 due to zero number indexing to find where
        # the alignment sequence actually starts
        if(source[j]=="AA"){space_diff[[loci[i]]]<-(nchar(strsplit(alignment[[loci[i]]][3], " ")[[1]][2])+countSpaces(alignment[[loci[i]]][3])[2]+1)-countSpaces(alignment[[loci[i]]][2])[1]}

        # reduce repeated whitespace in alignment file and removes rows with empty values for proper
        # start and stop subsetting
        alignment[[loci[i]]] <-str_squish(alignment[[loci[i]]])
        alignment[[loci[i]]] <-alignment[[loci[i]]][-which(alignment[[loci[i]]] == "")]

        # determine positions of "cDNA" or "Prot" and the end of that reference block segment
        start[[loci[i]]] <-as.numeric(grep(type, alignment[[loci[i]]]))
        end[[loci[i]]] <- as.numeric(c(start[[loci[i]]][2:length(start[[loci[i]]])]-1,length(alignment[[loci[i]]])))

        if(source[j]=="AA"){
          # extract rows with "Prot" and reference sequence position information
          # extract only relevant reference sequence positions
          # NOTE: the first row containing "Prot" contains two numbers -- -30 and 1 -- where only -30, is extracted,
          # as the actual sequence start will always be 1
          for (k in 1:length(start[[loci[i]]])){

            prot_extractions[[loci[i]]][k]<-strsplit(alignment[[loci[i]]][start[[loci[i]]][k]], " ")

            refblock_number[[loci[i]]][k]<-as.numeric(sapply(prot_extractions[[loci[i]]][k], "[", 2))


            # determine the alignment start by adding -30 to the difference between white spaces found above
            alignment_start[[loci[i]]]<-refblock_number[[loci[i]]][1]+space_diff[[loci[i]]]
          }
        } else if(source[j]=="cDNA"){
          # determine the alignment start by finding the second vector in second list and removing "codon"
          if(expressed[[i]]) { # exclude loci with cDNA alignment but no AA alignment (e.g., not expressed)
          alignment_start[[loci[i]]] <-as.numeric(sub("AA codon ", "", alignment[[loci[i]]][2]))
            }
          if(is.null(alignment_start[[loci[i]]])) {
          alignment_start[[loci[i]]] <- 1
            }
          DNA_start[[loci[i]]] <-as.numeric(sub("cDNA ", "", alignment[[loci[i]]][1]))
                 } else if(source[j]=="gDNA"){
                     DNA_start[[loci[i]]] <-as.numeric(sub("gDNA ", "", alignment[[loci[i]]][1]))
                                        }

        # close all white space in the alignment file, except for the white space separating the allele and peptide sequence
        alignment[[loci[i]]] <-paste(substr(alignment[[loci[i]]],1,regexpr(" ",text = alignment[[loci[i]]],fixed = TRUE)), gsub(" ","",substr(alignment[[loci[i]]],regexpr(" ",text = alignment[[loci[i]]],fixed = TRUE),nchar(alignment[[loci[i]]]))),sep = "")

        # string split at white spaces to yield allele and peptide sequences
        alignment[[loci[i]]]  <- strsplit(alignment[[loci[i]]]," ", fixed=T)

        # bind the previously split strings by row
        alignment[[loci[i]]] <- do.call(rbind,alignment[[loci[i]]])

        # if the seq column is equal to the allele column due to premature peptide termination,
        # insert a blank in place of the allele in the seq column (AA table only)
        if(source[j]=="AA"){alignment[[loci[i]]][which(alignment[[loci[i]]][,1]==alignment[[loci[i]]][,2]),2] <- ""} ##### I think this can be eliminated, since it is not in the AA block

        #renames columns to "alleles" and "seq"
        colnames(alignment[[loci[i]]])<-c(paste(loci[[i]], "alleles", sep="_"), "seq")

        # Due to ANHIG formatting, cases where an allele contains newly reference peptide sequences will not
        # contain the same number of rows as previous reference peptide blocks
        # this for loop is invoked to add "." for all other alleles for each character in the newly reference peptide
        # to preserve structural integrity
        # changes (10/9/19) to accommodate if there is more than one extraneous allele with an extended amino acid sequence

        if(source[j]=="cDNA"|source[j]=="AA"&loci[i]=="TAP2"){ ## the logic in the first expression is odd, but it works
          for(l in 1:length(start[[loci[i]]])){
            if(nrow(alignment[[i]][start[[loci[i]]][l]:end[[loci[i]]][l],])!=nrow(alignment[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],])){
              x<-as.data.frame(alignment[[loci[i]]][,1][start[[loci[i]]][1]:end[[loci[i]]][1]][-c(1,2)], stringsAsFactors = F)
              colnames(x)<-paste(loci[[i]], "alleles", sep="_")
              x<-cbind.data.frame(x, seq=as.character(paste(rep(".", nchar(tail(alignment[[loci[i]]][,2], 1))), collapse = "")), stringsAsFactors=FALSE)
              y<-data.frame(tail(alignment[[loci[i]]], (nrow(alignment[[i]][start[[loci[i]]][l]:end[[loci[i]]][l],][nrow(alignment[[i]][start[[loci[i]]][l]:end[[loci[i]]][l],])!=nrow(alignment[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],]),])-2)), stringsAsFactors = F)
              x$seq[match(y[,1], x[,1])]<-y$seq
              alignment[[loci[i]]]<-as.matrix(rbind(head(alignment[[loci[i]]], -(nrow(alignment[[i]][start[[loci[i]]][l]:end[[loci[i]]][l],][nrow(alignment[[i]][start[[loci[i]]][l]:end[[loci[i]]][l],])!=nrow(alignment[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],]),])-2)), x))
              start[[loci[i]]]<-as.numeric(grep(type, alignment[[loci[i]]]))
              end[[loci[i]]] <- as.numeric(c(start[[loci[i]]][2:length(start[[loci[i]]])]-1,nrow(alignment[[loci[i]]])))}
          }
        }
         # Like HLA-C, TAP2 has some alleles with extra nuc and prot 'tails' at the end of the alignment page.
        if(source[j]=="AA"|source[j]=="gDNA"|source[j]!="AA"&loci[i]!="TAP2"){
          for(l in 1:length(start[[loci[i]]])){
            if(nrow(alignment[[i]][start[[loci[i]]][l]:end[[loci[i]]][l],])!=nrow(alignment[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],])){ ####<<<< for HLA-C this is broken
              x<-as.data.frame(alignment[[loci[i]]][,1][start[[loci[i]]][1]:end[[loci[i]]][1]][-c(1,2)], stringsAsFactors = F)
              colnames(x)<-paste(loci[[i]], "alleles", sep="_")
              temp_vec<-alignment[[i]][start[[loci[i]]][l]:end[[loci[i]]][l],][-c(1,2),]
              cat(l,"\n")
              # If there is only one allele with an extended amino acid sequence, the conversion from a named
              # vector to a data frame is not properly done
              # If the character vector is = 2, there is only one allele with an extended amino acid sequence, so use
              # as.list and then data.frame so data frame is created correctly
              if(length(temp_vec) == 2){
                temp_filter<-data.frame(as.list(temp_vec), stringsAsFactors = F)
              } else{
                temp_filter<-data.frame(temp_vec, stringsAsFactors = F)
              }

              # find the greatest number of peptides in extraneous alleles to determine
              # how many "." to add on
              max_nchar<-(temp_filter %>%
                            add_column(nchar = nchar(.$seq)) %>%
                            filter(nchar == max(nchar)))$nchar
              x<-cbind.data.frame(x, seq=as.character(paste(rep(".", max_nchar), collapse = "")), stringsAsFactors=FALSE)
              y<-data.frame(tail(alignment[[loci[i]]], (nrow(alignment[[i]][start[[loci[i]]][l]:end[[loci[i]]][l],][nrow(alignment[[i]][start[[loci[i]]][l]:end[[loci[i]]][l],])!=nrow(alignment[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],]),])-2)), stringsAsFactors = F)
              x$seq[match(y[,1], x[,1])]<-y$seq
              alignment[[loci[i]]]<-as.matrix(rbind(head(alignment[[loci[i]]], -(nrow(alignment[[i]][start[[loci[i]]][l]:end[[loci[i]]][l],][nrow(alignment[[i]][start[[loci[i]]][l]:end[[loci[i]]][l],])!=nrow(alignment[[loci[i]]][start[[loci[i]]][1]:end[[loci[i]]][1],]),])-2)), x))
              start[[loci[i]]]<-as.numeric(grep("Prot", alignment[[loci[i]]]))
              end[[loci[i]]] <- as.numeric(c(start[[loci[i]]][2:length(start[[loci[i]]])]-1,nrow(alignment[[loci[i]]])))
              }
          }
        } ## end of not-TAP2 block

        # if a locus has extra formatting, resulting in unqequal rows, start and end will be updated to reflect subsetting
        # if a locus has no extra formatting, start and end will remain the same, as procured by earlier code
        for(m in 1:length(start[[loci[i]]])){
          HLAalignments[[loci[i]]]<-cbind(HLAalignments[[loci[i]]], alignment[[loci[i]]][start[[loci[i]]][m]:end[[loci[i]]][m],])
        }

        # remove rows containing "cDNA", "AA codon", or "Prot"
        HLAalignments[[loci[i]]] <- HLAalignments[[loci[i]]][-delete_lines,] ## this was the pseudo/frag problem; With no AAcodon row, it takes one too many lines out. Now corrected.

        # designate columns to be combined as every other so allele names are not included
        # in pasting all the amino acid sequences together
        cols<-seq(0, ncol(HLAalignments[[loci[i]]]), by=2)
        HLAalignments[[loci[i]]]<-cbind(HLAalignments[[loci[i]]][,1], apply(HLAalignments[[loci[i]]][,cols], 1 ,paste, collapse = "")) # The full-length alignment, with 1 row per allele

        # create a new matrix with the number of columns equal to the number of characters in the reference sequence
        corr_table[[loci[i]]]<-matrix(, nrow = 3, ncol = as.numeric(nchar(HLAalignments[[loci[i]]][,2][1])))

        if(source[j]=="AA"|source[j]=="cDNA"){  # closes at line 384(381)
          # if the first position enumeration is negative (i.e. has a leader peptide sequence), determines alignment length based on the total number of characters plus the alignment start (which is negative)
          if(grepl("-", alignment_start[[loci[i]]][[1]])==TRUE){
            alignment_length[[loci[i]]]<-(as.numeric(nchar(HLAalignments[[loci[i]]][,2][1]))/divide)+(alignment_start[[loci[[i]]]])
              } else { # if there is no leader peptide (i.e sequence starts at 1), determines alignment length based on total number of characters
            alignment_length[[loci[i]]]<-as.numeric(nchar(HLAalignments[[loci[i]]][,2][1]))
          }
          if(expressed[[i]]){   # pertinent to cDNA alignments with expressed proteins (codons are nucleotide triplets)
          # paste alignment_start to alignment_length together in sequential order
          # capture output as "w"
          w[[loci[i]]] <- capture.output(cat(alignment_start[[loci[i]]]:alignment_length[[loci[i]]]))

          # split string formed by cat for separate character variables
          alignment_positions[[loci[i]]]<-as.character(unlist(strsplit(w[[loci[i]]], " ")))

          # eliminate "0" if present in the alignment positions, as the alignment sequence from ANHIG does not contain 0
          if("0" %in% alignment_positions[[loci[i]]]==TRUE){
            alignment_positions[[loci[i]]]<-alignment_positions[[loci[i]]][-which(alignment_positions[[loci[i]]] == 0)]
                }
            } # closing the if not a pseudogene or fragment (i.e., expressed)
        } # closing AA/cDNA

        if(expressed[[i]]) {
        # triple alignment to account for codons (cDNA table only) # Does not pertain to pseduogenes and fragments
        if(source[j]=="cDNA"){
          for(n in 1:length(alignment_positions[[loci[i]]])){
            alignment_positionsx3[[loci[i]]]<-c(alignment_positionsx3[[loci[i]]],alignment_positions[[loci[i]]][n],alignment_positions[[loci[i]]][n],alignment_positions[[loci[i]]][n])
          }
          alignment_positions[[loci[i]]]<-alignment_positionsx3[[loci[i]]]
        }

        } ## closing the expressed block

        # string splits to extract locus in the allele name
        # assign to new variable "aligned"
        aligned[[loci[i]]]<- as.matrix(do.call(rbind,strsplit(HLAalignments[[loci[i]]][,1],"[*]")))

        # add a new column of pasted locus and trimmed two field alleles to aligned
        aligned[[loci[i]]]<- cbind(aligned[[loci[i]]], paste(aligned[[loci[i]]][,1], apply(aligned[[loci[i]]],MARGIN=c(1,2),FUN=getField,res=2)[,2], sep="*"))

        # bind aligned and HLAalignments -- rename columns
        HLAalignments[[loci[i]]] <- cbind(aligned[[loci[i]]], HLAalignments[[loci[i]]])
        colnames(HLAalignments[[loci[i]]]) <- c("locus", "full_allele", "trimmed_allele", "allele_name", sequence_name)

        if(source[j]=="AA"){
          # if locus is DRB3/4/5, use DRB line 1 as reference sequence
          if(loci[[i]]=="DRB3"|loci[[i]]=="DRB4"|loci[[i]]=="DRB5"){
            # set refexon to a reference peptide for each HLA locus based on the reference sequences in HLAalignments
            refexon[[loci[i]]] <- rbind(HLAalignments[[loci[i]]][1,])[which(rbind(HLAalignments[[loci[i]]][1,])[,"locus"]=="DRB1"),sequence_name]
          } else{
            refexon[[loci[i]]] <- rbind(HLAalignments[[loci[i]]][1,])[which(rbind(HLAalignments[[loci[i]]][1,])[,"locus"]==loci[[i]]),sequence_name]
          }

          # if input loci is DRB, use grep statement to match input loci to loci in HLAalignments
          if(loci[[i]]=="DRB"){
            refexon[[loci[i]]]<-rbind(HLAalignments[[loci[i]]][1,])[grepl(loci[[i]], rbind(HLAalignments[[loci[i]]][1,])[,"locus"]), sequence_name]}
        } # close the AA specific block

        # split sequence column at every amino acid or nucleotide, resulting in a split list of each for each row
        pepsplit[[loci[i]]] <- sapply(HLAalignments[[loci[i]]][,sequence_name],strsplit,split="*")

        # fill in space with NA for alleles with premature termination to make it the same number of characters
        # as the reference sequence (AA table only)
        if(source[j]=="AA"){pepsplit[[loci[i]]]<- lapply(pepsplit[[loci[i]]],function(x) c(x,rep("NA",nchar(refexon[[loci[i]]])-length(x))))}

        # nullify variable names
        names(pepsplit[[loci[i]]]) <- NULL

        # bind pep_split together by element in its previous list form by row
        pepsplit[[loci[i]]]<- do.call(rbind,pepsplit[[loci[i]]])

        # bind all columns together to form desired output, as described above
        HLAalignments[[loci[i]]] <- cbind.data.frame(HLAalignments[[loci[i]]][,1:4],pepsplit[[loci[i]]], stringsAsFactors=FALSE)

        if(source[j]=="AA"|source[j]=="cDNA"){
          # get reference sequence from DRB alignment
          if(loci[[i]] %in% c("DRB3","DRB4","DRB5")){
            DRBref<-HLAalignments[[loci[[i]]]][1,] #}

           #if the locus is DRB1/3/4/5, subset the specific locus from DRB alignment
             HLAalignments[[loci[i]]]<-HLAalignments[[loci[i]]][HLAalignments[[loci[[i]]]]$locus==loci[[i]],] #}

          # add reference sequence to DRB3/4/5, reset row names to numerical order
            HLAalignments[[loci[[i]]]]<- rbind(DRBref, HLAalignments[[loci[[i]]]])
            rownames(HLAalignments[[loci[[i]]]])<-NULL}
        } ## close DRB/DRB1 specific AA/cDNA block

        # input HLAalignments alignment sequence into the corr_table with "InDel" still present
        corr_table[[loci[[i]]]][1,]<-names(HLAalignments[[loci[[i]]]][5:ncol(HLAalignments[[loci[[i]]]])])

        # find positions in HLAalignments that have ".", indicating an inDel
        inDels[[loci[[i]]]]<-colnames(HLAalignments[[loci[[i]]]][1, 5:ncol(HLAalignments[[loci[[i]]]])][HLAalignments[[loci[[i]]]][1, 5:ncol(HLAalignments[[loci[[i]]]])] %in% "."])

        # find positions in HLAalignments that have "|", indicating a exon boundary
        exonB[[loci[[i]]]]<-colnames(HLAalignments[[loci[[i]]]][1, 5:ncol(HLAalignments[[loci[[i]]]])][HLAalignments[[loci[[i]]]][1, 5:ncol(HLAalignments[[loci[[i]]]])] %in% "|"])

        # inDel inclusion if there are inDels present
        if(length(inDels[[loci[[i]]]])!=0){
          for(o in 1:length(inDels[[loci[[i]]]])){
            corr_table[[loci[[i]]]][2,][inDels[[loci[[i]]]][[o]]==corr_table[[loci[[i]]]][1,]]<-paste("INDEL", o, sep="-")
            corr_table[[loci[[i]]]][3,][inDels[[loci[[i]]]][[o]]==corr_table[[loci[[i]]]][1,]]<-paste("INDEL", o, sep="-")
          }
        } # close indel block

        # exonB inclusion if there are exon boundaries present
        if(length(exonB[[loci[[i]]]])!=0){
          for(p in 1:length(exonB[[loci[[i]]]])){
            corr_table[[loci[[i]]]][2,][exonB[[loci[[i]]]][[p]]==corr_table[[loci[[i]]]][1,]]<-paste("EXONB", p, sep="-")
            corr_table[[loci[[i]]]][3,][exonB[[loci[[i]]]][[p]]==corr_table[[loci[[i]]]][1,]]<-paste("EXONB", p, sep="-")
          }
        }

        if(expressed[[i]]) {
        if(source[j]=="AA"|source[j]=="cDNA"){
          # paste alignment_positions into corr_table accounting for InDels and exon boundaries
          codon_count<-as.vector(1)
          for(q in 1:length(corr_table[[loci[[i]]]][2,])){
            if(is.na(corr_table[[loci[[i]]]][2,q])){
              corr_table[[loci[[i]]]][2,q]<-alignment_positions[[loci[i]]][codon_count]
               codon_count<-codon_count+1
            }
          }
        }
  }

        # create number sequence for alignment (cDNA and gDNA tables only)
        if(source[j]=="cDNA"){ alignment_positions[[loci[i]]]<-as.character(unlist(strsplit(capture.output(cat(DNA_start[[loci[i]]]:length(corr_table[[loci[[i]]]][3,]))), " ")))
        } else if(source[j]=="gDNA"){ if(DNA_start[[loci[i]]] < 0) {alignment_positions[[loci[i]]]<-as.character(unlist(strsplit(capture.output(cat(c(DNA_start[[loci[i]]]:-1,1:length(corr_table[[loci[[i]]]][3,])))), " "))) ## most loci have negative start positions
        } else {alignment_positions[[loci[i]]]<-as.character(unlist(strsplit(capture.output(cat(c(DNA_start[[loci[i]]]:length(corr_table[[loci[[i]]]][3,])))), " ")))  #### HLA-P starts at position 1, but in theory some other locus could start at a different positive position
                            }
                         } # In 3.53.0, only the HLA-P alignment started with 1, so, using the original version above, the first four alignment_positions values are "1"   "0"   "-1"  "1"
                           # pastes number sequence into corr_table accounting for InDels and exon boundaries (cDNA and gDNA tables only)
        if(source[j]=="cDNA"|source[j]=="gDNA"){ # cDNA/gDNA correspondence table block
          cDNAcount<-as.vector(1)
          for(r in 1:length(corr_table[[loci[[i]]]][3,])){
            if(is.na(corr_table[[loci[[i]]]][3,r])){
              corr_table[[loci[[i]]]][3,r]<-alignment_positions[[loci[i]]][cDNAcount]
              cDNAcount<-cDNAcount+1
            }
          }
        } # Close cDNA/gDNA block

        # run if InDels or exon boundaries are present in alignment
        if(any(grepl("INDEL|EXONB", corr_table[[loci[[i]]]][2,]))){ ## if there are no boundaries or indels, this never gets run so N and U cDNA don't get their corr_tables made
          # find which positions have InDels or exon boundaries
          v<-as.numeric(corr_table[[loci[[i]]]][1,][which(grepl("INDEL", corr_table[[loci[[i]]]][2,]))])
          w<-as.numeric(corr_table[[loci[[i]]]][1,][which(grepl("EXONB", corr_table[[loci[[i]]]][2,]))])

          # split cumulative sequences
          vsplit<-split(v, cumsum(c(TRUE, diff(v) != 1L))) # absolute corr_table coordinates (row 1)
          wsplit<-split(w, cumsum(c(TRUE, diff(w) != 1L))) # absolute corr_table coordinates (row 1)

          if(any(grepl("INDEL", corr_table[[loci[[i]]]][2,]))){ # indel block
            # add increments of .1 to any InDels in row 2
            if(source[j]=="AA"|source[j]=="cDNA"){ # AA/cDNA indel block
              for(s in 1:length(vsplit)){
                if(grepl("-", corr_table[[loci[[i]]]][2,][vsplit[[s]][[1]]-1])) {
                  corr_table[[loci[[i]]]][2,][vsplit[[s]]]<-paste(corr_table[[loci[[i]]]][2,][vsplit[[s]][[1]]-1], gsub(0, "", seq(1, length(vsplit[[s]]))/10), sep="")
                    } else {
                  corr_table[[loci[[i]]]][2,][vsplit[[s]]]<-as.numeric(corr_table[[loci[[i]]]][2,][vsplit[[s]][[1]]-1])+as.numeric(paste(0,".",seq(1, length(vsplit[[s]])),sep = ""))
                  }
              }
            } # close AA/cDNA indel block

            # add increments of .1 to any InDels in row 3
            if(source[j]=="cDNA"|source[j]=="gDNA"){ # cDNA/gDNA block
              for(t in 1:length(vsplit)){
                corr_table[[loci[[i]]]][3,][vsplit[[t]]]<-paste(corr_table[[loci[[i]]]][3,][vsplit[[t]][[1]]-1],".",seq(1, length(vsplit[[t]])),sep = "")
              }
            } # close cDNA/gDNAA block
          } # close indel block

          if(any(grepl("EXONB", corr_table[[loci[[i]]]][2,]))){
            # change EXONB to E.(# of last exon)-(# of next exon)
            if(source[j]=="cDNA"){  # make Exon -> E.
              for(u in 1:length(wsplit)){
                corr_table[[loci[[i]]]][2,][wsplit[[u]]]<-paste("E.",u,"-",u+1, sep = "")
                corr_table[[loci[[i]]]][3,][wsplit[[u]]]<-paste("E.",u,"-",u+1, sep = "")
              }
            } # close Exon -> E.

  boundaryIndel <- w[which((w+1) %in% v)] # Identify the position of the ExonBoundary immediately after which an indel occurs
  if(length(boundaryIndel) != 0) {        # These are special procedures for HLA-C and -DQB1 that should be generalizable
  insertion <- FALSE
  toReplace <- toTarget <- NA
            for(ff in 1:length(boundaryIndel)) { # boundary Indel Loop
              if(length(boundaryIndel[ff]) != 0) { # then we have an indel immediately following a boundary
                vsplitSub <- NA
                for(o in 1:length(vsplit)) {
                  if((boundaryIndel[ff]+1) %in% vsplit[[o]]) {vsplitSub <- o}
                  if(!is.na(vsplitSub)) {
                  toReplace <- vsplit[vsplitSub][[as.character(o)]][length(vsplit[vsplitSub][[as.character(o)]])]+1
                  if(length(grep("E",corr_table[[loci[[i]]]][2,toReplace])) != 0) {insertion <- TRUE}
                  toTarget <- vsplit[vsplitSub][[as.character(o)]][1]
                  vsplitSub <- NA
                  }
                }
 #             } # close length(boundaryIndel) loop --- THIS NEEDS TO CLOSE AFTER THE BELOW
        if(insertion) { # Currently pertains only to DQB1 as of 3.54.0
          corr_table[[loci[[i]]]][2,toTarget] <- corr_table[[loci[[i]]]][3,toTarget] <- "Ins"
            } else {
                    corr_table[[loci[[i]]]][2,toTarget] <- corr_table[[loci[[i]]]][2,toReplace]
                    corr_table[[loci[[i]]]][3,toTarget] <- corr_table[[loci[[i]]]][3,toReplace]
                }
              } # close length(boundaryIndel[ff]) loop
            } # close boundaryIndel loop
      } # close boundaryIndel > 0 check

            # change EXONB to corresponding boundary e.g. UTR-E.1
            # create atlas of gDNA boundaries and their locus
            if(source[j]=="gDNA"){ # start of the gDNA construction
              if(loci[[i]] %in% names(fragmentFeatureNames)) { # simplified gDNA atlas builder for genes in the fragemtnFeatureNames object

                gDNA_atlas[[loci[[i]]]]<-matrix(, nrow = 1, ncol = length(wsplit))
                gDNA_atlas[[loci[[i]]]]<- as.numeric(corr_table[[i]][3,w+1])
                newNames <- rep(NA,length(fragmentFeatureNames[names(fragmentFeatureNames) %in% loci[[i]] == TRUE][[1]][[1]])-1)
                    for(f in 1:(length(newNames))) {
                          newNames[f] <- paste(fragmentFeatureNames[names(fragmentFeatureNames) %in% loci[[i]] == TRUE][[1]][[1]][f],"-",fragmentFeatureNames[names(fragmentFeatureNames) %in% loci[[i]] == TRUE][[1]][[1]][f+1], sep="")
                    }
                gDNA_atlas[[loci[[i]]]] <- as.data.frame(t(gDNA_atlas[[loci[[i]]]]))
                colnames(gDNA_atlas[[loci[[i]]]]) <- newNames
                rownames(gDNA_atlas[[loci[[i]]]]) <- "gDNA"

               } else { intron<-as.vector(paste("I.",seq(1,((length(wsplit)-1)/2)-0.5),sep = "")) # } # start of the 'expressed gene' notation section
              exon<-as.vector(paste("E.",seq(1,((length(wsplit)-1)/2)+0.5),sep = ""))
              a<-unlist(strsplit(paste0(exon," ",intron), " "))
              a<-c("U.5",head(a, length(a)-1),"U.3")
              b<-vector()
              for(y in 1:(length(a)-1)){
                b<-c(b,paste(a[[y]],"-",a[[y+1]],sep = ""))}
              gDNA_atlas[[loci[[i]]]]<-matrix(, nrow = 1, ncol = length(wsplit))
              # cover 3 loci
              if(length(b) != length(gDNA_atlas[[loci[[i]]]])) {
                gDNA_atlas[[loci[[i]]]]<-matrix(, nrow = 1, ncol = length(b))
              }
              colnames(gDNA_atlas[[loci[[i]]]])<-b
              rownames(gDNA_atlas[[loci[[i]]]])<-"gDNA"
              for(z in 1:length(wsplit)){
                gDNA_atlas[[loci[[i]]]][1,z]<-corr_table[[loci[[i]]]][3,][wsplit[[z]]+1]
                corr_table[[loci[[i]]]][3,][wsplit[[z]]]<-b[[z]]
              }
             } # End of expressed gene construction
            } # End of the gDNA construction routine

            # create atlas of cDNA boundaries and their locus
            if(source[j]=="cDNA"){

              # if there are no boundaries, make a 1 column atlas
              if(length(exonB[[loci[[i]]]]) == 0) { # if exonB[[loci[[i]]]] has a length of 0, then it has no boundaries, and this alignment represents one feature (HLA-N and -U only as of 3.53.0).

                cDNA_atlas[[loci[[i]]]]<-matrix(, nrow = 1, ncol = 1)
                colnames(cDNA_atlas[[loci[i]]]) <- "No.Exon.Boundaries"
                rownames(cDNA_atlas[[loci[[i]]]])<-c("cDNA")
                cDNA_atlas[[loci[[i]]]][1,1]<-corr_table[[loci[[i]]]][2,1]

              } else {

              cDNA_atlas[[loci[[i]]]]<-matrix(, nrow = 2, ncol = length(wsplit))
              colnames(cDNA_atlas[[loci[[i]]]])<-corr_table[[loci[[i]]]][2,][unlist(wsplit)]
              rownames(cDNA_atlas[[loci[[i]]]])<-c("cDNA","codon")
              for(z in 1:length(wsplit)){
                cDNA_atlas[[loci[[i]]]][1,z]<-corr_table[[loci[[i]]]][3,][wsplit[[z]]+1]
                cDNA_atlas[[loci[[i]]]][2,z]<-corr_table[[loci[[i]]]][2,][wsplit[[z]]+1]
              }
             } # closes if length exon B = 0's else
            } # closes if(source[j] == "cDNA")
          } # closes the if(any(grepl("EXONB" block
        } # closes the if any INDEL|EXONB block

    # Create a cDNA Atlas for N and U (and any other cDNA that doesn't have boundaries)
     if((source[j] == "cDNA") & (length(exonB[[loci[[i]]]]) == 0)) {   # if exonB[[loci[[i]]]] has a length of 0, then it has no boundaries, and this alignment represents one feature (HLA-N and -U only as of 3.53.0).

            cDNA_atlas[[loci[[i]]]]<-matrix(, nrow = 1, ncol = 1)
            colnames(cDNA_atlas[[loci[i]]]) <- "No.Exon.Boundaries"
            rownames(cDNA_atlas[[loci[[i]]]])<-c("cDNA")
            cDNA_atlas[[loci[[i]]]][1,1]<-corr_table[[loci[[i]]]][2,1]

     }
                           # pull the AA atlas from the cDNA atlas
        if(subSource[[j]] == "AA") { #return things back to the way they should be after jumping the AA atlas routine
          source[[j]] <- "AA" ## restore the proper source
          AA_atlas[[loci[[i]]]] <- as.data.frame(cDNA_atlas[[loci[[i]]]])[2,]
          rownames(AA_atlas[[loci[i]]]) <- c("AA")
          atlas[[loci[i]]] <-c(atlas[[loci[i]]], AA = list(AA_atlas[[loci[i]]]),`Version` = alignmentVersion[[loci[[1]]]])
          subSource[[j]] <- 0 ## clear the submitted source
        }

        if(source[j]=="cDNA"){ # With no boundaries, the Atlas is empty
          atlas[[loci[i]]] <-c(atlas[[loci[i]]], cDNA = list(cDNA_atlas[[loci[i]]]),`Version` = alignmentVersion[[loci[[i]]]])
        }else if(source[j]=="gDNA"){
          atlas[[loci[i]]] <-c(atlas[[loci[i]]], gDNA = list(gDNA_atlas[[loci[i]]]),`Version` = alignmentVersion[[loci[[i]]]])
        } # close cDNA & gDNA creation code

        }
      
      if(source[j]=="cDNA" && loci[i] %in% c("R","S","T","V","W")) { ## post-hoc update of cDNA feature names for fragments and pseudogenes 
        
        repLen <- length(colnames(atlas[[loci[i]]]$cDNA))
        featNames <- fragmentFeatureNames[[names(fragmentFeatureNames)[names(fragmentFeatureNames) %in% loci[i]]]]$feature
        featRep <- featNames[which(1:length(featNames)%%2==0)]
        newPairs <- rep(NA,length(featRep)-1)
        
            for(g in 1:length(newPairs)){
              newPairs[g] <- paste(featRep[g],featRep[g+1],sep="-")
            }
        
        colnames(atlas[[loci[i]]]$cDNA) <- newPairs
      }
      
    } # end of the loop if j source values
  } # end of main loop of i length values

  for(k in 1:length(loci)){
    if("AA" %in% names(atlas[[loci[k]]])) { atlas[[loci[k]]]$AA <- as.data.frame(atlas[[loci[k]]]$AA)   }
    if("cDNA" %in% names(atlas[[loci[k]]])) { atlas[[loci[k]]]$cDNA <- as.data.frame(atlas[[loci[k]]]$cDNA)   }
    if("gDNA" %in% names(atlas[[loci[k]]])) { atlas[[loci[k]]]$gDNA <- as.data.frame(atlas[[loci[k]]]$gDNA)   }
  }

  return(atlas)
  } else { return("The 'locus' and/or 'source' parameters were invalid.") } ## end of the LocusValidator check
}
