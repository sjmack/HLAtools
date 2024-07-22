## HLAtools: Functions and Datasets for HLA Informatics

## Version 1.1.1

The "Human Leukocyte Antigen" (HLA) region is the most polymorphic section of the human genome, with 39,886 allelic variants identified across 46 loci. The key roles played by the class I and class II HLA genes in stem-cell therapy and transplantation, HLA and disease association research, evolutionary biology, and population genetics results in constant discovery of new allele variants. These data are curated and maintained by the [ImmunoPolymorphism Database-IMmunoGeneTics/HLA (IPD-IMGT/HLA) Database](https://www.ebi.ac.uk/ipd/imgt/hla/) and made available on the [Anthony Nolan HLA Informatics Group (ANHIG)/IMGTHLA GitHub repository](https://github.com/ANHIG/IMGTHLA) as static text files, which are updated every three months. Standardized use of the data in this key resource can be challenging. To address this, we have developed HLAtools, an R package that automates the consumption of IPD-IMGT/HLA resources, renders them computable, and makes them available alongside tools for data analysis, visualization and investigation. This version of the package is compatible with all IPD-IMGT/HLA Database release versions up to release 3.57.0.

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

In addition, the _alignmentFull()_ function builds the 'HLAalignments' object, which includes computable versions of the protein, codon, coding nucleotide and genomic alignments available in the [IMGTHLA GitHub repository](https://github.com/anhig/IMGTHLA), as specified by the user. 'HLAalignments' is not included the package, but can be built for IPD-IMGT/HLA Database releases 3.00.0 to 3.57.0.

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

### Data Format Conversion and Translation Funtions
The package includes functions that convert [Genotype List (GL) String Codes](https://glstring.org) across IPD/IMGT-HLA Database release versions and nomenclature epochs, inter-convert between [GL String](https://glstring.org) and [UNIFORMAT](https://hla-net.eu/tools/uniformate/) formats, and translate data frames and vectors of HLA allele name data across IPD-IMGT/HLA Database release versions.  

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

- GIANT() translates vectors and data frames of HLA allele names across specified IPD-IMGT/HLA Database release versions. 
```
GIANT(c("A*01:01:01:01","DQA1*01:01:01:01"),"3.56.0","3.00.0")
[1] "A*01:01:01:01" "DQA1*01:01:01"

GIANT(c("A*01:01:01:01","DQA1*01:01:01:01"),"3.56.0","2.20.0")
[1] "A*01010101"  "DQA1*010101"

GIANT(c("A*01:01:01:01","DQA1*01:01:01:01"),"3.56.0","1.09.0")
[1] "A*01011"   "DQA1*0101"

GIANT(sHLAdata,"3.56.0","2.23.0")[1:5,1:14]
   Subject Status        A      A.1      C    C.1      B    B.1   DRB1 DRB1.1   DQA1 DQA1.1   DQB1 DQB1.1
1 UT900-23      0     <NA>     <NA> 010201 020205   1301 180102 160201   0404 050101 050102 030101   0202
2 UT900-24      0 01010101 02010101   0307   0605   1401 390201   0402 160201 040101 050101 060101 030101
3 UT900-25      0     0210   030102   0712 010201   1520   1301 080201 040701   0201 050102 030302   0202
4 UT900-26      0 01010101     0218   0804 120201 350901   4005 090102 080201 030101   0201   0304 030302
5 UT910-01      0   250101 02010101   1507   0307 510103   1401 110101 090102 050102 030101 050301   0304
```

### Data Analysis Functions
The package includes three data-analysis functions that accept Bridging ImmunoGenomic Data-Analysis Workflow Gaps (BIGDAWG)-formatted genotype datasets as input.

To facilitate testing and experimentation with these functions, the 'sHLAdata' data object is included in the package. This BIGDAWG-formatted HLA genotype data object represents completely synthetic IPD-IMGT/HLA release 3.56.0 HLA-A, -C, -B, -DRB1, -DQA1, -DQB1, -DPA1, and -DPB1 genotype data for 24 control subjects and 23 case subjects, and does not represent any true human population.

- relRisk() calculates relative risk (RR) values, confidence interval (CI) values and p-values for BIGDAWG-formatted non-case-control genotype datasets. For these analyses, two subject categories are required, but should not be affected/patient and unaffected/control categories; instead, the categories may be, e.g., either of two disease states, where one disease state is coded as 0 and the other is coded as 1 in the second column of the dataset.

```
rr <- relRisk(sHLAdata)

rr$alleles[[1]][c(1,3),]
  Locus     Variant Status_1 Status_0      RelativeRisk            CI.low          CI.high    p.value Significant
1     A 01:01:01:01        5        7 0.833333333333333 0.411641930809013 1.68701094924813 0.59291160            
2     A 02:01:01:01        8       11 0.830409356725146 0.467284553053116  1.4757168736504 0.50779167            
3     A       02:10        4        5               0.9 0.419645026972693 1.93020278553833 0.77978869            

rr$genotypes[[1]][c(1:3),]
  Locus                 Variant Status_1 Status_0 RelativeRisk     CI.low     CI.high      p.value Significant
1     A 01:01:01:01+02:01:01:01        0        2            0          0         NaN  0.161777420            
2     A       01:01:01:01+02:10        1        0   2.0952380   1.5379358  2.85449005  0.306556273            
3     A       01:01:01:01+02:18        0        4            0          0         NaN  0.042730188           *
```

- BDstrat() stratifies BIGDAWG-formatted case-control datasets for individual alleles or multiple alleles at multiple loci, and generates two BIGDAWG-formatted datasets; one for the case and control subjects that have those alleles, and one for the case and control subjects that do not.

```
HLA_data.multi.strat <- BDstrat(sHLAdata,c("DRB1*16:02:01:01","DRB1*04:07:01:01","A*25:01:01:01"),15)

HLA_data.multi.strat$`DRB1*16:02:01:01+DRB1*04:07:01:01+A*25:01:01:01-positive`[1:2,1:8]
   Subject Status           A         A.1           C      C.1           B         B.1
1 UT900-23      0        <NA>        <NA> 01:02:01:01 02:10:06 13:01:01:01    18:01:02
2 UT900-24      0 01:01:01:01 02:01:01:01 03:07:01:01    06:05 14:01:01:01 39:02:01:01

HLA_data.multi.strat$`DRB1*16:02:01:01+DRB1*04:07:01:01+A*25:01:01:01-negative`[1:2,1:8]
   Subject Status           A   A.1           C         C.1           B         B.1
4 UT900-26      0 01:01:01:01 02:18 08:04:01:01    12:02:01 35:09:01:01 40:05:01:01
6 UT910-02      0       02:10 32:04 18:01:01:01 01:02:01:01    78:02:01 13:01:01:01
```

- BDtoPyPop() converts a BIGDAWG-formatted case-control dataset into two PyPop version 1.\*.\* formatted datasets -- one for all 'case' subjects and one for all 'control' subjects.  

```
HLAdata.PP <- BDtoPyPop(sHLAdata,"sHLA",FALSE)

HLAdata.PP$sHLA.neagtive[1:3,1:10]
   Subject Status         A_1         A_2         C_1         C_2         B_1         B_2      DRB1_1      DRB1_2
1 UT900-23      0        ****        **** 01:02:01:01    02:10:06 13:01:01:01    18:01:02 16:02:01:01 04:04:01:01
2 UT900-24      0 01:01:01:01 02:01:01:01 03:07:01:01       06:05 14:01:01:01 39:02:01:01    04:02:01 16:02:01:01
3 UT900-25      0       02:10    03:01:02       07:12 01:02:01:01       15:20 13:01:01:01 08:02:01:01 04:07:01:01

HLAdata.PP$sHLA.positive[1:3,1:10]
    Subject Status      A_1      A_2         C_1         C_2         B_1         B_2      DRB1_1      DRB1_2
25 UT910-30      1    02:18    02:18 01:02:01:01 03:07:01:01 40:05:01:01 14:01:01:01 11:01:01:01 04:07:01:01
26 UT910-31      1     ****     **** 03:07:01:01        **** 14:01:01:01    78:02:01    04:02:01 08:02:01:01
27 UT910-32      1 03:01:02 03:01:02       06:05 03:07:01:01 39:02:01:01 14:01:01:01 08:02:01:01 08:02:01:01

```

