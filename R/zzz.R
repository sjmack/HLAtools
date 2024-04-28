
#### Preliminary OnLoad

        onLoad <- function(libname, pkgname) {
  
               pseudo.codon <<- c("DPA2","DPB2","Y") 
              # Prior to release version 3.53.0, the cDNA alignments for the DPA2 and DPB2 pseudogenes included 'AA codon' positions. These lines were removed in release version 3.54.0.
              # Similarly, the HLA-Y pseudogene had an AA codon line prior to version 3.36.0.
              # These exceptions are consumed by atlasMaker(), directing the removal of 'AA codon' lines for these cDNA alignments.
            
          }

        onUnload <- function(libname, pkgname) {
 
              objs <- ls(pos = ".GlobalEnv")

              rm(list = objs[grep("pseudo.codon",objs)], pos = ".GlobalEnv")
          
          }

        .onLoad <-function(libname, pkgname) {
  
              onLoad(libname, pkgname)
 
          }

        .onUnload <- function(libpath,pkgname){
  
              onUnload(libpath,pkgname)
          }
