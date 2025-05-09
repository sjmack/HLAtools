% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/atlasMaker.R
\name{atlasMaker}
\alias{atlasMaker}
\title{Identify the Gene-Feature Boundaries in HLA Region Genes}
\usage{
atlasMaker(loci, source, version = "Latest")
}
\arguments{
\item{loci}{A character string identifying a specific HLA gene (ex. "DRB1"). More than one gene can be specified (e.g., c("DRB1","DQB1")).}

\item{source}{A character string identifying the type of alignment used to build the atlas. The allowed options are "AA", for protein alignments, "cDNA", for nucleotide alignments, and "gDNA", for genomic alignments. More than one alignment type can be specified (e.g., c("AA","gDNA")).}

\item{version}{A character string identifying the release version of the IPD-IMGT/HLA Github repository for which alignments should be generated.}
}
\value{
A list object of 'atlas' dataframes and a 'release version' character string for each locus and alignment type.
}
\description{
The 'AA' atlas identifies the amino acid (AA) positions encoded by codons that span each exon (E) boundary. The 'cDNA' atlas describes the 3' nucleotide position in the cDNA sequence that follows each exon boundary, as well as the codon that spans each exon boundary. The 'gDNA' atlas describes first 3' nucleotide position following each UTR (U), exon or intron (I) boundary. Each feature is followed by its rank (e.g U3 is the 3' UTR). Non-standard gene-features (H, J, N, and S) are detailed in the documentation for the ffN() function.
}
\details{
Builds an 'atlas' of the gene-feature (exon, intron and UTR) boundaries for IPD-IMGT/HLA loci, using IPD-IMGT/HLA Database alignments.
}
\note{
For internal HLAtools use.

Nucleotide atlases for pseudogenes will include a 'codon' row populated with NA values.

Some HLA-DQB1\*05 and \*06 alleles include a 5th exon that is not present in the DQB1\*05:01:01:01 reference allele. In this case, the Exon 4 to Exon 5 boundary is defined as 'Ins' (insertion). For all other alleles there is no E.4-5 boundary.

Boundaries for non-standard hybrid (H), join (J), novel (N) and segment (S) features may be included in gene fragment and pseudogene atlases.
}
