## HLAtools (development version)

### Version 0.9.14.9000 June 10 to June 12, 2024

- Updated buildAlignments() to build DPA2 and DPB2 cDNA alignments in releases prior to 3.54.0.
- Updated buildAlignments() to build DRB1 gDNA alignments in release 3.48.0.
- Updated buildAlignments() to build HLA-B cDNA alignments in releases 3.44.0 and 3.43.0.
- Updated buildAlignments() to build HLA-N, -S, -U and -Y cDNA alignments in releases prior to 3.36.0.
- Updated buildAlignments() to build HLA-W and -T cDNA alignments in release 3.27.0.
- Updated buildAlignments() to build HFE cDNA alignments in releases prior to 3.28.0.
- Updated buildAlignments() to build DPA AA alignments in release 3.23.0. While the returned object is called "$DPA", the locus in this object is "DPA1".


### Version 0.9.13.9000 June 4, 2024

- Clarified that the oldest available release for GL String Code updating is version 1.05.0.
- Changed the default IMGT/HLA Release version for alignmentFull() to the version of the loaded HLAgazeteer.

### Version 0.9.12.9000 May 28, 2024

- Added motifMatch() and validateMotif() functions. Thanks to Kazu Osoegawa for the suggestion.
- Updated documentation for checkAlignType().
- Corrected documentation for fragmentFeatureNames.

### Version 0.9.11.9000 May 16, 2024

- Updated package title.
- Minor updates to Vignette.
- Changed '*' to '#' for PyPop version references.
- Complete citations for published sources of the HLAgazeteer.
- Limited asterisks in function parameter documentation.
- Added link to PyPop configuration file webpage in pypopHeaders().
- Added locus sanity check in compareSequences().
- Modified getField() so that NA values are not appended to an allele name when a resolution higher than that of the provided allele is requested.

### Version 0.9.10.9000 May 10, 2024 

- Updated buildGazeteer() to account for the DRB2, DRB6, DRB7 and DRB9 gene.
- Updated documentation for HLAgazeteer regarding the DRB2, DRB6, DRB7 and DRB9 genes and the DRB alignment files. 
- Updated documentation for fragmentFeatureNames() regarding the DRB2, DRB6, DRB7 and DRB9 genes and genomic alignments.
- Updated the vignette to reflect these issues.

### Version 0.9.9.9000 May 3, 2024 

- Added BDtoPyPop() function.
- Updates to Vignette, ReadMe and Description.

### Version 0.9.8.9000 April 27, 2024

- Removed reliance on HTexceptions data object.
- Correction to relRisk().

### Version 0.9.7.9000 April 20, 2024

- Standardized use of 'source' in documentation. 
- Changed 'alignType' parameter to 'source' in multiLocusValidation().
- Added checkAlignType() and checkSource() functions for converting between alignmentType (for alignment objects) and source (for alignment files) values.
- Added checkAlignType() to AlignmentFull().
- Suppressed warning messages in alignmentFull() and atlasMaker().
- Restricted compareSequences() to a single 'alignType'.
- Clarified when only a single 'alignType' is allowed for a function.
- Edits to Vignette.

### Version 0.9.6.9000 April 15, 2024

- Added the queryRelease() function for searching the alleleListHistory object.
- Formalized and streamlined the GL String to UNIFORMAT and UNIFORMAT to GL String functions and updated documentation.
- Corrected documentation for customAlign().
- Corrected Markdown formatting in BDstrat.R.
- First draft of complete Vignette.

### Version 0.9.5.9000 April 11, 2024

- Added 'source' selection to multiLocusValidation().
- Modified atlasMaker() to support updated multiLocusValidation().
- Updated packaged data objects (HLAatlas,HLAgazeteer,fragmentFeatureNames and alleleListHistory) to version IPD-IMGT/HLA Database Release version 3.56.0.
- Updated Vignette.

### Version 0.9.4.9000 April 9, 2024

- Added HTexceptions() and the HTexceptions data object, which make defined exceptions for specific use cases available for package functions. 
- Added HTexceptions.rda to the /data folder. 
- Updated atlasMaker() to use HTexceptions for building DPA2, DPB2 and HLA-Y nucleotide atlases.
- Updated atlasMaker() to use multiLocusValidation(), which depends on the HLAgazeteer$version to evaluate locus:version matches. 
- Modified multiLocusValidation() messaging to reflect this dependence on HLAgazeteer$version.
- Updates to Vignette.
- Updates to ReadMe.

### Version 0.9.3.9000 April 7, 2024

- Correctly added GPL3 License to package.
- Added getAlignmentNames() function.
- Added repoVersion() function.
- Added parseAlignmentHead() function.
- Updated buildGazeteer() to use getAlignmentNames().
- Updated atlasMaker(), buildAlignments() and buildGazeteer() to use repoVersion().
- Updates to Vignette.
- Updates to ReadMe.

### Version 0.9.2.9000 April 3, 2024

- Added multiLocusValidation() function.
- Applied multiLocusValidation() in AlignmentFull().
- Closed connections in versionValidation().
- Added documentation to fragmentFeatureNames.R.
- Added content to Vignette.
- Corrected returned list names in BDstrat().
- Updated ReadMe.
- Updated documentation for relRisk().

### Version 0.9.1.9000 April 1, 2024

- Minor update to alignmentFull().
- Updated alignmentFull() documentation to note the absence of DP and DQ gene numbers in older database releases.
- Added sections to Vignette.
- Reformatted News to include bullets.

### Version 0.9.0.9000 March 31, 2024

- Updated BDStrat() to function on multiple alleles at multiple loci. 
- Added sections to Vignette.

### Version 0.8.1.9000 March 24, 2024

- Added the BDstrat() function for stratifying BIGDAWG-formatted datasets by single alleles 
- Added the verifyAllele() function, which leverages AlleleListHistory to determine if and when a full alllele name was valid.
- Updated ReadMe.

### Version 0.8.0.9000 March 21, 2024

- Added the relRisk() function for calculating relative risk using BIGDAWG-formatted datasets.
- Added $align to HLAgazeteer to identify all genes with alignments.
- Added new versions of the HLAatlas, HLAgazeteerm IMGTHLAGeneTypes, alleleListHistory, and fragmentFeatureNames objects to the package.
- Updated ReadMe.
- Added bullet points to News.

### Version 0.7.7.9000 March 19, 2024

- Updated posSort() and and alignmentSearch() to address edge-cases and support searching for intron positions.
- Updated getField() so that expression variant identifiers are optionally appended to truncated versions of full-length expression variants.
- Reestructured descriptions of several functions for clarity. 

### Version 0.7.6.9000 March 19, 2024

- Standardized versioning for the packaged data objects.  
- Added functionality for updating packaged data objects.  

### Version 0.7.5.9000 March 13, 2024

- Expanded HLAgazeteer to define multiple functional and organizational sets of genes.  
- Expanded alignmentFull() to build user-defined sets of alignments.   
- Updated ReadMe and Description to reflect expanded capacities.  

### Version 0.7.4.9000 February 18, 2024

- Consolidated sequence query and custom alignment building functions to support all alignments.  

### Version 0.7.3.9000 February 16, 2024

- Added a nucleotide ('nuc') aligmment to HLAalignments, and changed the previous 'nuc' alignment to 'codon'.  
- Nucleotide sequence query functions support both nucleotide and codon cDNA aligmnents.  

### Version 0.7.2.9000 February 10, 2024

- Added documentation, updated functions for building alignments.  
  
### Version 0.7.1.9000 February 9, 2024

- HLAtools Package contains reference datasets, along side query and analysis tools.  
