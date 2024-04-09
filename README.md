## HLAtools: Functions and Datasets for Human Leukocyte Antigen Informatics

## Version 0.9.4.9000

The Human Leukocyte Antigen (HLA) region is the most polymorphic section of the human genome, with 38,008 allelic variants identified across 46 loci. The key roles played by the class I and class II HLA genes in stem-cell therapy and transplantation, HLA and disease association research, evolutionary biology, and population genetics results in constant discovery of new allele variants. These data are curated and maintained by the [IPD-IMGT/HLA Database](https://www.ebi.ac.uk/ipd/imgt/hla/) and made available on the [ANHIG/IMGTHLA GitHub repository](https://github.com/ANHIG/IMGTHLA) as static text files, which are updated every three months. Standardized use of the data in this key resource can be challenging. To address this, we have developed HLAtools, an R package that automates the consumption of IPD-IMGT/HLA resources, renders them computable, and makes them available alongside tools for data analysis and visualization.   

### Data Resources
The package includes five data objects that foster computation on IPD-IMGT/HLA resources. 

- 'IMGTHLAGeneTypes' describes the [named genes in the HLA region](https://hla.alleles.org/genes/index.html).
- 'HLAgazetteer' defines specific categories of genes supported by the IPD-IMGT/HLA Database. For example:
   - gene fragments (HLAgazeteer$frag : "N" "P" "S" "T" "U" "V" "W" "X" "Z") 
   - non-classical HLA genes (HLAgazeteer$nonclassical : "DMA"  "DMB"  "DOA"  "DOB"  "DPA2" "DPB2" "DQA2" "DQB2" "E"    "F"    "G")
- 'fragmentFeatureNames' identifies and annotates the non-standard features found in some gene fragments, based on the positions of feature boundaries ("|") in the sequence.
   - For example, the HLA-U gene fragment includes three features (fragmentFeatureNames\$U\$features : "N.1" "J.1" "S.1") that do not align to other class I gene features. The "N.1" feature is 54 nucleotides of novel (N) sequence, the "J.1" feature is a join (J) of the 3' end of Exon 3, and the 5' end of Intron 3, and the "S.1" feature is a segment (S) of Intron 3. 
- 'HLAatlas' identifies the boundaries between gene features (exons, introns and untranslated regions) at each gene, pseudogene and gene fragment.
- 'alleleListHistory' is a computable index of the [names of all HLA alleles](https://github.com/ANHIG/IMGTHLA/tree/Latest/allelelist) for all IPD-IMGT/HLA Database release versions. 

HLAgazeteer, HLAatlas, and alleleListHistory can be updated with each IPD-IMGT/HLA Database release. 

In addition, the _alignmentFull()_ function builds the 'HLAalignments' object, which includes computable versions of the protein, codon, coding nucleotide and genomic alignments available in the [IMGTHLA GitHub repository](https://github.com/anhig/IMGTHLA), as specified by the user. 'HLAalignments' is not included the package, but can be updated with each IPD-IMGT/HLA Database release.

### Search and Query Functions
The package includes a suite of functions for dissecting and describing similarities and differences between alleles and across loci.

- compareSequences() identifies the positions that differ between a pair of alleles at a locus.
```
compareSequences(alignType = "gen", alleles = c("DPA1*01:03:01:04","DPA1*01:03:38:01"))
       allele_name 51 1544 1723 3318 4149
1 DPA1*01:03:01:04  C    A    A    G    C
2 DPA1*01:03:38:01  T    G    G    C    G
```
- customAlign() builds customized peptide, codon, nucleotide and genomic alignments for specific alleles at user-specified positions, for alleles at different loci, and for different positions at each locus. 
```
customAlign("prot",c("DQA1*01:01:01:01","DQB1*05:01:01:01","DPB1*01:01:01:01"),list(1:4,5:8,5:8))
              DQA1 1 2 3 4
1 DQA1*01:01:01:01 E D I V
2             DQB1 5 6 7 8
3 DQB1*05:01:01:01 E D F V
4             DPB1 5 6 7 8
5 DPB1*01:01:01:01 E N Y V

customAlign("codon",c("DRB1*01:01","DQB1*02:01","DPB1*01:01"),c(1,2,3,7,8,9,13,14,15))
      Allele   1   2   3   7   8   9  13  14  15
1 DRB1*01:01 GGG GAC ACC TTC TTG TGG TTT GAA TGT
2 DQB1*02:01 AGA GAC TCT TTC GTG TAC GGC ATG TGC
3 DPB1*01:01 AGG GCC ACT TAC GTG TAC CAG GAA TGC

customAlign("prot",c("DPB1*01:01:01:01","DQA1*01:01:01:01","DQB1*05:01:01:01"),list(19:35,1:4,5:8))
              DPB1 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35
1 DPB1*01:01:01:01  N  G  T  Q  R  F  L  E  R  Y  I  Y  N  R  E  E  Y
2             DQA1  1  2  3  4                                       
3 DQA1*01:01:01:01  E  D  I  V                                       
4             DQB1  5  6  7  8                                       
5 DQB1*05:01:01:01  E  D  F  V                                       
```

- Additional functions include alleleTrim(), which trims HLA allele-names by fields or digits, validateAllele(), which determines if the specified allele-name is present in the 'HLAalignments' object that has been loaded in the R environment, and verifyAllele(), which determines if the specified allele-name is present in the 'AlleleListHistory' object, and optionally identifies the most recent IPD-IMGT/HLA Database release including that allele.

```
alleleTrim(allele = "A*03:01:01", resolution = 2)
[1] "A*03:01"

alleleTrim(allele = "A*030101", resolution = 2,version = 2)
[1] "A*0301"

validateAllele("A*01:01:01:117")
[1] TRUE

validateAllele("A*01:01:01")
A*01:01:01 is not found in version 3.55.0 alignments.
[1] FALSE

verifyAllele("A*01:01:01:01")
[1] TRUE

verifyAllele("A*01:01:01:01",TRUE)
[1] "TRUE"   "3.55.0"

verifyAllele("A*010101",TRUE)
[1] "TRUE"   "2.09.0"

verifyAllele("A*0101",TRUE)
[1] "TRUE"   "1.06.0"
```

### Format Conversion Funtions
The package includes functions that convert [Genotype List (GL) String Codes](https://glstring.org) across IPD/IMGT-HLA Database release versions and nomenclature epochs, and that inter-convert between [GL String](https://glstring.org) and [UNIFORMAT](https://hla-net.eu/tools/uniformate/) formats.  

- GLudpdate() converts HLA allele names in GL String Code objects across IPD/IMGT-HLA Database release versions. 

```
GLupdate("hla#1.07.0#HLA-B*35011", "2.05.0")
[1] "hla#2.05.0#HLA-B*350101"

GLupdate("hla#2.05.0#HLA-B*350101", "3.05.0")
[1] "hla#3.05.0#HLA-B*35:01:01:01"

GLupdate("hla#3.55.0#HLA-B*44:02:01:01","1.05.0")
[1] "hla#1.05.0#HLA-B*4402"

GLupdate("hla#3.25.0#HLA-A*01:01:01:01/HLA-A*01:02+HLA-A*24:02:01:01","3.55.0")
[1] "hla#3.55.0#HLA-A*01:01:01:01/HLA-A*01:02:01:01+HLA-A*24:02:01:01"

GLupdate("hla#3.25.0#HLA-A*01:01:01:01/HLA-A*01:02+HLA-A*24:02:01:01","1.05.0")
[1] "hla#1.05.0#HLA-A*0101/HLA-A*0102+HLA-A*2402101"
```
- GLStoUNI() and UNItoGLS() translate between the GL String and UNIFORMAT grammars
```
GLStoUNI("HLA-A*02:01/HLA-A*02:02+HLA-A*03:01/HLA-A*03:02")
[1] "A*02:01,A*03:01|A*02:01,A*03:02|A*02:02,A*03:01|A*02:02,A*03:02"

UNItoGLS("A*02:01,A*03:01|A*02:01,A*03:02|A*02:02,A*03:01|A*02:02,A*03:02")
[1] "HLA-A*02:01/HLA-A*02:02+HLA-A*03:01/HLA-A*03:02"
```

### Data Analysis Functions
The package includes two data-analysis functions that accept BIGDAWG-formatted genotype datasets as input. 

- relRisk() calculates relative risk (RR) values, confidence interval (CI) values and p-values for non-case-control genotype data. For these analyses, two subject categories are required, but should not be affected/patient and unaffected/control categories; instead, the categories may be, e.g., either of two disease states, where one disease state is coded as 0 and the other is coded as 1 in the second column of the dataset.

- BDstrat() stratifies BIGDAWG-formatted case-control datasets for individual alleles or multiple alleles at multiple loci, and generates two BIGDAWG-formatted datasets; one for the case and control subjects that have those alleles, and one for the case and control subjects that do not.




