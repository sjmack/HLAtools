## HLAtools (development version)

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
  
### Version 0.7.1.9000 Ferbrary 9, 2024

- HLAtools Package contains reference datasets, along side query and analysis tools.  
