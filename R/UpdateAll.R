## UpdateAll v02.1.0 21 July 2024

################
##updateAll
#'Update All Package Data Objects Derived from IPD-IMGT/HLA Database Resources
#'
#'@description
#'Applies updateAlleleListHistory(), atlasFull(), buildGazeteer(), extractGeneTypes(), and ffN() to update the alleleListHistory, HLAatlas, HLAgazeteer, IMGTHLAGeneTypes and fragmentFeatureNames data objects.
#'
#'A new alleleListHistory data object should be generated with each IPD-IMGT/HLA Database release. The other data objects will likely only change when new genes are added to the IPD-IMGT/HLA Database.
#'
#'@param updateType A character vector of the names of data objects to be updated. By default, updateAll() builds all five data objects (updateType="all"). Alternatively specific data objects can be updated; e.g., updateType="alleleListHistory" or updateType=c("alleleListHistory","fragmentFeatureNames").
#'@param version A numeric value or character string identifying the version of the ANHIG/IMGTHLA Github repository to build these objects from. By default, updateAll() calls the getLatestVersion() function to identify the most recent IPD-IMGT/HLA Database release.
#'
#'@return No value is returned. The desired versions of the specified data objects are built into the environment that called updateAll().
#'
#'@note Generating a new HLAatlas can take 5 minutes or more to complete.
#'
#'@export
updateAll <-function(updateType="all",version = getLatestVersion()){

    HLTpath <- path.package("HLAtools")
    HLTDpath <- paste(HLTpath,"/data",sep="")
    
  #check and see if a version of an object exists, build new local objects, then load them in turn into the environment

  #IMGTHLAGeneTypes
    if("all" %in% updateType || "IMGTHLAGeneTypes" %in% updateType) {
      rawPage <- readLines("https://hla.alleles.org/genes/index.html",-1,warn = FALSE)
      IHGTversion <- strsplit(strsplit(rawPage[grep("#BeginDate",rawPage,fixed=TRUE)][1],"-->",fixed = TRUE)[[1]][2],"<!--",fixed=TRUE)[[1]][1]
      
      if(IHGTversion == IMGTHLAGeneTypes$version){
        message(paste("IMGTHLAGeneTypes for version",IHGTversion,"is already loaded.",sep=" "))
      } else { 
        
        if(file.exists(paste(HLTDpath,paste(IHGTversion,"IMGTHLAGeneTypes.rda",sep="."),sep="/"))) {
        load(paste(HLTDpath,paste(IHGTversion,"IMGTHLAGeneTypes.rda",sep="."),sep="/"),envir = parent.frame())
        message(paste("IMGTHLAGeneTypes for version",IHGTversion,"has been loaded.",sep=" "))
      } else { 
        
      IMGTHLAGeneTypes <- buildIMGTHLAGeneTypes() # has version, no version argument
        IHGTversion <- IMGTHLAGeneTypes$version
      save(IMGTHLAGeneTypes,file=paste(HLTDpath,paste(IHGTversion,"IMGTHLAGeneTypes.rda",sep="."),sep="/"))
      load(paste(HLTDpath,paste(IHGTversion,"IMGTHLAGeneTypes.rda",sep="."),sep="/"),envir = parent.frame())
      message(paste("IMGTHLAGeneTypes for version",IHGTversion,"has been built and loaded.",sep=" "))
          }
       }
    }
      
    if("all" %in% updateType || "HLAgazeteer" %in% updateType) {
      
      if(version == HLAgazeteer$version){
        message(paste("HLAgazeteer for version",version,"is already loaded.",sep=" "))
      } else {
        if(file.exists(paste(HLTDpath,paste(version,"HLAgazeteer.rda",sep="."),sep="/"))) {
          load(paste(HLTDpath,paste(version,"HLAgazeteer.rda",sep="."),sep="/"),envir = parent.frame())
          message(paste("HLAgazetter for version",version,"has been loaded.",sep=" "))
        } else { 
        
      HLAgazeteer <- buildGazeteer(version) ## has version, takes version
        HGversion <- HLAgazeteer$version
      save(HLAgazeteer,file=paste(HLTDpath,paste(HGversion,"HLAgazeteer.rda",sep="."),sep="/"))
      load(paste(HLTDpath,paste(HGversion,"HLAgazeteer.rda",sep="."),sep="/"),envir = parent.frame())
      message(paste("HLAgazeteer for version",version,"has been built and loaded.",sep=" "))
        }
       }
    }
    
    if("all" %in% updateType || "alleleListHistory" %in% updateType) {
      if(version == alleleListHistory$Version$Version){
        message(paste("alleleListHistory for version",version,"is already loaded.",sep=" "))
      } else {
        if(file.exists(paste(HLTDpath,paste(version,"alleleListHistory.rda",sep="."),sep="/"))) {
          load(paste(HLTDpath,paste(version,"alleleListHistory.rda",sep="."),sep="/"),envir = parent.frame())
          message(paste("alleleListHistory for version",version,"has been loaded",sep=" "))
        } else { 
          
      alleleListHistory <- updateAlleleListHistory() ## has version, takes no versiomn argument
        ALHversion <- alleleListHistory$Version$Version
      save(alleleListHistory,file=paste(HLTDpath,paste(ALHversion,"alleleListHistory.rda",sep="."),sep="/"))
      load(paste(HLTDpath,paste(ALHversion,"alleleListHistory.rda",sep="."),sep="/"),envir = parent.frame())
      message(paste("alleleListHistory for version",version,"has been built and loaded",sep=" "))
        }
      }
    } 

    if("all" %in% updateType || "fragmentFeatureNames" %in% updateType) {
      if(version == fragmentFeatureNames$version){
        message(paste("fragmentFeatureNames for version",version,"is already loaded.",sep=" "))
      } else {      
        if(file.exists(paste(HLTDpath,paste(version,"fragmentfeatureNames.rda",sep="."),sep="/"))) {
          load(paste(HLTDpath,paste(version,"fragmentFeatureNames.rda",sep="."),sep="/"),envir = parent.frame())
          message(paste("fragmentFeatureNames for version",version,"has been loaded.",sep=" "))
        } else {  
        
      fragmentFeatureNames <- ffN(version) ## has version, takes version
        FFNversion <- fragmentFeatureNames$version
      save(fragmentFeatureNames,file=paste(HLTDpath,paste(FFNversion,"fragmentFeatureNames.rda",sep="."),sep="/"))
      load(paste(HLTDpath,paste(FFNversion,"fragmentFeatureNames.rda",sep="."),sep="/"),envir = parent.frame())
      message(paste("fragmentFeatureNames for version",version,"has been built and loaded.",sep=" "))
        }
      }
    }  
     
    if("all" %in% updateType || "HLAatlas" %in% updateType) {
      if(version == HLAatlas$version){
        message(paste("HLAatlas for version",version,"is already loaded.",sep=" "))
      } else {      
        if(file.exists(paste(HLTDpath,paste(version,"HLAatlas.rda",sep="."),sep="/"))) {
          load(paste(HLTDpath,paste(version,"HLAatlas.rda",sep="."),sep="/"),envir = parent.frame())
          message(paste("HLAatlas for version",version,"has been loaded.",sep=" "))
        } else {  
          
      HLAatlas <- atlasFull(version) ## has version, takes version
        HAversion <- HLAatlas$version
      save(HLAatlas,file=paste(HLTDpath,paste(HAversion,"HLAatlas.rda",sep="."),sep="/"))
      load(paste(HLTDpath,paste(HAversion,"HLAatlas.rda",sep="."),sep="/"),envir = parent.frame())
      message(paste("HLAatlas for version",version,"has been built and loaded.",sep=" "))
        }
      }
    }
  ##
}
