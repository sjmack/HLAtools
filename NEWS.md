### HLAtools

### Version 1.1.1

- July 21, 2024
- Added link to ANHIG/IMGTHLA GitHub repo in Description.
- Removed example from BuildIMGTHLAGeneTypes().
- Modified updateAll() to load new objects into the parent.frame() environment, rather than .GlobalEnv.
- Submission to CRAN.

### Version 1.1.0

- July 17, 2024
- Added the GIANT() function for translating vectors and data frames of HLA allele names across IPD-IMGT/HLA Database release versions.
- Updates to ReadMe and Vignette.
- Updates to function examples.
- Submission to CRAN.

### Version 1.0.4

- July 12, 2024
- Modified functions to support and incorporate the new HLA-R gene.
- Corrected fragmentFeaturesNames annotations for HLA-P and -W.
- Updated the alleleListHistory, fragmentFeatureNames, HLAgazeteer and HLAatlas data objects to version 3.57.0, and IMGTHLAGeneTypes to version 09-07-2024, and aded these versions to the package.

### Version 1.0.3

- July 10, 2024
- Corrected typo in queryRelease() documentation.
- Expanded acronyms in DESCRIPTION.
- Added 'value' lines to addCodonLine(), getField(), typeToSource() and updateAll().
- Minor edits to other function descriptions.
- Replaced 'dontrun' with 'donttest' in function documentations, or remove 'dontrun' in some cases.
- Removed note for atlasFull().
- Removed pseudo.codon global variable and zzz.R file.
- Updated GLVhelper() for full release versions 1.*, 2.* and 3.* functionality. 
- Streamlined updateGL() functionality and include informative examples.
- Updated documentation and examples for GLupdate(),  GLV(), GLV2, redec(). 
- Updated GLVhelper() functoinality, documention and messaging.
- Added 'namespace' parameter to GLvalidate(), supporting HLA and KIR GLSC namespaces, although KIR is not supported in the HLAtools package.
- Added BIGDAWG-formatted sHLAdata synthetic dataset for use with BDstrat(), BDtoPyPop() and relRisk().
- Added calls to sHLAdata in  BDstrat(), BDtoPyPop() and relRisk() examples.
- Updated ReadMe and Vignette.

### Version 1.0.2

- June 28, 2024
- Removed GPL license document.
- Adopted canonical URL to BIGDAWG input format in relRisk() and BDStrat().
- Added URL to HLAtools GitHub repository.
- Correct titles for ExpandVersion() and GLtoUN() functions.
- Submission to CRAN.

### Version 1.0.1

- June 27, 2024
- Wrapped individual authors in c() in DESCRIPTION.
- Limited 'cre' role to the package Maintainer.
- Submission to CRAN.

### Version 1.0.0

- June 26, 2024
- Standardized function documentation.
- Remove write.table() calls from atlasMaker().
- Add save.path parameter to BDtoPyPop() and relRisk().
- initial package release.
- Submission to CRAN.

### Version 0.9.16.9000 June 24th to 25th, 2024

- Updated atlasMaker() to build atlases for DRB1 in release 3.48.0.
- Updated atlasMaker() to build atlases for HLA-B in releases 3.44.0 and 3.43.0.
- Updated atlasMaker() to build atlases for HFE in releases 3.27.0 to 3.22.0. 
- Updated atlasMaker() to build atlases for HLA-V in release 3.14.0.
- Updated atlasMaker() to build atlases for DPA2 and DPB2 in releases 3.53.0-3.27.0.
- Updated atlasMaker() to build atlases for HLA-N, -S, -U, and -Y in releases 3.35.0-3.33.0.
- Updated atlasMaker() to build atlases for HLA-Y in releases 3.32.0-3.20.0.
- Updated atlasMaker() to build atlases for HLA-W and -T in release 3.27.0.
- Updated atlasMaker() to build atlases for All loci in release 3.13.0/3.13.1.
- Updated atlasMaker() to build atlases for  HLA-C in release 3.02.0.
- Updated atlasMaker() to build atlases for DPA, DPB, TAP1 and TAP2 in release 3.00.0.
- Stopped atlasMaker() writing correspondence_tables to tempdir.

### Version 0.9.15.9000 June 21, 2024

- Added addCodonLine() function to facilitate fixing nucleotide alignments that are missing "AA codon" lines.

### Version 0.9.14.9000 June 10th to 20th, 2024

- Updated buildAlignments() to build DPA2 and DPB2 cDNA alignments in releases prior to 3.54.0.
- Updated buildAlignments() to build DRB1 gDNA alignments in release 3.48.0.
- Updated buildAlignments() to build HLA-B cDNA alignments in releases 3.44.0 and 3.43.0.
- Updated buildAlignments() to build HLA-N, -S, -U and -Y cDNA alignments in releases prior to 3.36.0.
- Updated buildAlignments() to build HLA-W and -T cDNA alignments in release 3.27.0.
- Updated buildAlignments() to build HFE cDNA alignments in releases prior to 3.28.0.
- Updated buildAlignments() to build HLA-DPA, -DPB, -DQA and -DQB AA alignments in releases prior to 3.24.0. The returned objects are "$DPA", "$DPB", "$DQA" and "$DQB", the loci in these objects are "DPA1", "DPB1", "DQA1" and "DQB1".
- Updated buildAlignments() to build HLA-V in version 3.14.0.
- Updated buildAlignments() to account for version inconsistency (3.13.0 vs 3.13.1) in release 3.13.1 alignments.
- Updated documentation of repoVersion() for accuracy.
- Applied repoVersion() in buildAlignments().
- Modified buildGazeteer() to account for version differences in URL and alignment files for release 3.13.0/3.13.1.
- Updated buildAlignments() to build HLA-C cDNA alignments in release 3.02.0.
- Close connections in HLAgazeteer().
- Updated buildAligmments() to build HLA-DOA AA alignment in release 3.01.0.
- Updated documentation for updateAll().
- Updated repoVersion() for release 3.00.0.
- Updated buildAlignments() to build HLA-DPA and -DPB cDA alignments in release 3.00.0.
- Updated buildAlignments() to build TAP1 and TAP2 cDNA alignments in version 3.00.0.

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
