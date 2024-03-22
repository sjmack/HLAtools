## HLAtools: Functions and Datasets for Human Leukocyte Antigen Informatics

## Version 0.8.0.9000

The Human Leukocyte Antigen (HLA) region is the most polymorphic section of the human genome, with 38,008 allelic variants identified across 46 loci. The key roles played by the class I and class II HLA genes in stem-cell therapy and transplantation, HLA and disease association research, evolutionary biology, and population genetics results in constant discovery of new allele variants. These data are curated and maintained by the [IPD-IMGT/HLA Database](https://www.ebi.ac.uk/ipd/imgt/hla/) and made available on the [ANHIG/IMGTHLA GitHub repository](https://github.com/ANHIG/IMGTHLA) as static text files, which are updated every three months. Standardized use of the data in this key resource can be challenging. To address this, we have developed HLAtools, an R package that automates the consumption of IPD-IMGT/HLA resources, renders them computable, and makes them available alongside tools for data analysis and visualization.   

### Data Resources
The package includes five data objects that foster computation on IPD-IMGT/HLA resources. 

- 'IMGTHLAGeneTypes' describes the [named genes in the HLA region](https://hla.alleles.org/genes/index.html).
- 'HLAgazetteer' defines specific categories of genes supported by the IPD-IMGT/HLA Database. For example:
   - gene fragments (HLAgazeteer$frag : "N" "P" "S" "T" "U" "V" "W" "X" "Z") 
   - non-classical HLA genes (HLAgazeteer$nonclassical : "DMA"  "DMB"  "DOA"  "DOB"  "DPA2" "DPB2" "DQA2" "DQB2" "E"    "F"    "G")
- 'HLAatlas' identifies the boundaries between gene features (exons, introns and untranslated regions) at each gene, pseudogene and gene fragment.
- 'fragmentFeatureNames' identifies and annotates the non-standard features found in some gene fragments, based on the positions of feature boundaries ("|") in the sequence.
- For example, the HLA-U gene fragment includes three features (fragmentFeatureNames\$U\$features : "N.1" "J.1" "S.1") that do not align to other class I gene features. The "N.1" feature is 54 nucleotides of novel (N) sequence, the "J.1" feature is a join (J) of the 3' end of Exon 3, and the 5' end of Intron 3, and the "S.1" feature is a segment (S) of Intron 3. 
- 'alleleListHistory' is a computable index of the [names of all HLA alleles](https://github.com/ANHIG/IMGTHLA/tree/Latest/allelelist) for all IPD-IMGT/HLA Database release versions. 

HLAgazeteer, HLAatlas, and alleleListHistory can be updated with each IPD-IMGT/HLA Database release. 

In addition, the _alignmentFull()_ function builds the 'HLAalignments' object, which includes computable versions of the protein, codon, coding nucleotide and genomic alignments available in the [IMGTHLA GitHub repository](https://github.com/anhig/IMGTHLA), as specified by the user. 'HLAalignments' is not included the package, but can be updated with each IPD-IMGT/HLA Database release.

### Search and Query Functions
The package includes a suite of functions for dissecting and describing similarities and differences between alleles and across loci.

- The compareSequences() function identifies the positions that differ between a pair of alleles at a locus.

```
compareSequences(alignType = "gen", alleles = c("DPA1*01:03:38:01","DPA1*01:03:38:02"))
       allele_name 1544 1723 3318 4149
1 DPA1*01:03:38:01    G    G    C    G
2 DPA1*01:03:38:02    A    A    G    C
```

- The customAlign() function builds customized peptide, codon, nucleotide and genomic alignments for specific alleles at user-specified positions, for alleles at different loci. 

```
customAlign("codon",c("DRB1*01:01","DQB1*02:01","DPB1*01:01"),c(1,2,3,7,8,9,13,14,15))
      Allele   1   2   3   7   8   9  13  14  15
1 DRB1*01:01 GGG GAC ACC TTC TTG TGG TTT GAA TGT
2 DQB1*02:01 AGA GAC TCT TTC GTG TAC GGC ATG TGC
3 DPB1*01:01 AGG GCC ACT TAC GTG TAC CAG GAA TGC
```

### Format Conversion Funtions
alidation, trimming and conversion of HLA allele names across IPD-IMGT/HLA database releases, and conversion between [GL String](https://glstring.org) and [UNIFORMAT](https://hla-net.eu/tools/uniformate/) data formats.  

### Data Analysis Functions
The HLAtools package includes two data-analysis functions. Both use BIGDAWG-formatted genotype datsets as input. 
- relRisk() calculates relative risk (RR) values, confidence interval (CI) values and p-values for non-case-control genotype data. For these analyses, two subject categories are required, but should not be patients and controls; instead, the categories may be, e.g., either of two disease states, where one disease state is coded as 0 and the other is coded as 1 in the second column of the dataset.
- BDstrat() uses BIGDAWG to perform a stratified analyses on a BIGDAWG-formatted case-control dataset. For these analyses, the user specifes an allele at a locus, and BDstrat performs two BIGDAWG analyses; one for all of the case and control subjects that have that allele, and one for all of the case and control subjects that do not.    
