## HLAtools: Functions and Datasets for Human Leukocyte Antigen Informatics

## Version 1.0.2

The Human Leukocyte Antigen (HLA) region is the most polymorphic section of the human genome, with 39,886 allelic variants identified across 46 loci. The key roles played by the class I and class II HLA genes in stem-cell therapy and transplantation, HLA and disease association research, evolutionary biology, and population genetics results in constant discovery of new allele variants. These data are curated and maintained by the [IPD-IMGT/HLA Database](https://www.ebi.ac.uk/ipd/imgt/hla/) and made available on the [ANHIG/IMGTHLA GitHub repository](https://github.com/ANHIG/IMGTHLA) as static text files, which are updated every three months. Standardized use of the data in this key resource can be challenging. To address this, we have developed HLAtools, an R package that automates the consumption of IPD-IMGT/HLA resources, renders them computable, and makes them available alongside tools for data analysis, visualization and investigation. This version of the package is compatible with all IPD-IMGT/HLA Database release versions up to release 3.56.0.

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

In addition, the _alignmentFull()_ function builds the 'HLAalignments' object, which includes computable versions of the protein, codon, coding nucleotide and genomic alignments available in the [IMGTHLA GitHub repository](https://github.com/anhig/IMGTHLA), as specified by the user. 'HLAalignments' is not included the package, but can be built for IPD-IMGT/HLA Database releases 3.00.0 to 3.56.0.

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

- motifMatch() searches the HLAalignments object to identify full-length or two-field alleles that share user-specified nucleotide or peptide motifs.
```
motifMatch("A*-21M~2P","prot")
[1] "A*02:774"  "A*11:284"  "A*11:417N" "A*68:216N"

motifMatch("A*196G~301A~3046T","gen",FALSE)
[1] "A*01:09"

motifMatch("A*196G~301A~3046T","gen",TRUE)
[1] "A*01:09:01:01" "A*01:09:01:02"
```

- queryRelease() searches the alleleListHistory object for user-defined allele name variants in a specific IPD-IMGT/HLA release, identifying the number of alleles that match the query term, or a vector of allele names that match the query term. 

```
queryRelease("3.30.0","DRB9",FALSE) 
[1] 1

queryRelease("3.30.0","DRB9",TRUE) 
[1] "DRB9*01:01"

queryRelease("3.31.0","DRB9",FALSE) 
[1] 6

queryRelease("1.05.0","304",TRUE) 
[1] "A*0304"    "A*3304"    "B*1304"    "Cw*03041"  "Cw*03042"  "DQB1*0304" "DRB1*0304" "DRB1*1304" "B*5304"    "A*2304"   
```

- Additional functions include *alleleTrim()*, which trims HLA allele-names by fields or digits, *validateAllele()*, which determines if the specified allele-name is present in the 'HLAalignments' object that has been loaded in the R environment, and *verifyAllele()*, which determines if the specified allele-name is present in the 'AlleleListHistory' object, and optionally identifies the most recent IPD-IMGT/HLA Database release including that allele.

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

### Data Format Conversion Funtions
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

- GLStoUNI() and UNItoGLS() convert between the GL String and UNIFORMAT grammars
```
GLStoUNI("HLA-A*02:01/HLA-A*02:02+HLA-A*03:01/HLA-A*03:02")
[1] "A*02:01,A*03:01|A*02:01,A*03:02|A*02:02,A*03:01|A*02:02,A*03:02"

UNItoGLS("A*02:01,A*03:01|A*02:01,A*03:02|A*02:02,A*03:01|A*02:02,A*03:02")
[1] "HLA-A*02:01/HLA-A*02:02+HLA-A*03:01/HLA-A*03:02"
```

### Data Analysis Functions
The package includes three data-analysis functions that accept BIGDAWG-formatted genotype datasets as input. 

- relRisk() calculates relative risk (RR) values, confidence interval (CI) values and p-values for BIGDAWG-formatted non-case-control genotype datasets. For these analyses, two subject categories are required, but should not be affected/patient and unaffected/control categories; instead, the categories may be, e.g., either of two disease states, where one disease state is coded as 0 and the other is coded as 1 in the second column of the dataset.

```
library(BIGDAWG)
rr <- relRisk(HLA_data[,1:4])

rr$alleles[[1]][c(1,3),]
  Locus     Variant Status_1 Status_0  RelativeRisk        CI.low      CI.high       p.value Significant
1     A 01:01:01:01      176      166  1.0343294132 0.92848272156 1.1522425892 0.54580794275          
3     A    02:05:01      105      142 0.84368316144 0.72728310379 0.9787127917 0.01649632817           *

rr$genotypes[[1]][c(1,52),]
   Locus                 Variant Status_1 Status_0 RelativeRisk        CI.low      CI.high    p.value Significant
1      A 01:01:01:01+01:01:01:01        8        7 1.0693602693 0.66473577000 1.7202796017 0.789560336
52     A       02:05:01+26:01:01        1        7 0.2497492477 0.03990709706 1.5629973447 0.034058218          *
```

- BDstrat() stratifies BIGDAWG-formatted case-control datasets for individual alleles or multiple alleles at multiple loci, and generates two BIGDAWG-formatted datasets; one for the case and control subjects that have those alleles, and one for the case and control subjects that do not.

```
HLA_data.multi.strat <- BDstrat(BIGDAWG::HLA_data,c("DRB1*08:01:03","DRB1*03:01:02","A*26:08"))

HLA_data.multi.strat$`DRB1*08:01:03+DRB1*03:01:02+A*26:08-positive`[1:2,1:6]
  SampleID Disease           A   A.1        DRB1      DRB1.1
2  SCo0002       0 03:01:01:01 68:06    08:01:03 15:01:01:01
3  SCo0003       0       26:08 32:02 07:01:01:01 15:01:01:01

HLA_data.multi.strat$`DRB1*08:01:03+DRB1*03:01:02+A*26:08-negative`[1:2,1:6]
  SampleID Disease           A         A.1     DRB1   DRB1.1
1  SCo0001       0 01:01:01:01 01:01:01:01 01:01:01 01:01:01
4  SCo0004       0 01:01:01:01    32:01:01 01:01:01 11:04:01
```

- BDtoPyPop() converts a BIGDAWG-formatted case-control dataset into two PyPop version 1.\*.\* formatted datasets -- one for all 'case' subjects and one for all 'control' subjects.  

```
HLAdata.PP <- BDtoPyPop(BIGDAWG::HLA_data,"BDHLA",FALSE)

HLAdata.PP$BDHLA.neagtive[1:3,]
     SampleID Disease         A_1         A_2      DRB1_1      DRB1_2      DQB1_1      DQB1_2 DRB3_1 DRB3_2
1     SCo0001       0 01:01:01:01 01:01:01:01    01:01:01    01:01:01 05:03:01:01 05:03:01:01  00:00  00:00
2     SCo0002       0 03:01:01:01       68:06    08:01:03 15:01:01:01    03:02:12 03:01:01:01  00:00  00:00
3     SCo0003       0       26:08       32:02 07:01:01:01 15:01:01:01    03:02:01 03:01:01:01  00:00  00:00

HLAdata.PP$BDHLA.positive[1:3,]
     SampleID Disease         A_1         A_2   DRB1_1      DRB1_2      DQB1_1   DQB1_2      DRB3_1   DRB3_2
1003  SCa0001       1 01:01:01:01 11:01:01:01 01:01:01    04:01:01 05:03:01:01 03:02:01       00:00    00:00
1004  SCa0002       1       32:02    68:01:01 03:01:02    11:01:01    06:02:01 05:02:01 01:01:02:01 03:01:01
1005  SCa0003       1       32:02 11:01:01:02 01:01:01 07:01:01:01 05:03:01:01 02:02:01       00:00    00:00
```

