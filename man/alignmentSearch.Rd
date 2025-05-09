% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SearchandBuildAlignments.R
\name{alignmentSearch}
\alias{alignmentSearch}
\title{Search Alignments for Specific Positions in a Specific Allele}
\usage{
alignmentSearch(alignType, allelename, positions, prefix = TRUE, sep = "~")
}
\arguments{
\item{alignType}{The type of alignment being searched. Allowed values are "codon","gen", nuc" and "prot". Only one 'alignType' value is allowed.}

\item{allelename}{A full or 2-field HLA allele name, excluding the "HLA-" prefix.}

\item{positions}{A vector of sequence positions (e.g., c(-17,1,22,130)); in nucleotide and genomic alignments, indel positions are named using decimals. So the first indel between positions 26 and 27 is named 26.1, and the second indel between 26 and 27 is 26.2, etc.}

\item{prefix}{A logical that indicates if the position number(s) should be included in the result.}

\item{sep}{The value that will be used to separate the sequences for each position. No value can be specified (sep=""). If prefix=FALSE and sep="", a sequence block is returned.}
}
\value{
A character string containing the corresponding peptide, codon or nucleotide sequence for each position in 'positions' for 'allelename'. A sequence is not returned if 'allelename' is not not found in the pertinent alignment. A position will be empty if 'allelename' does not have a value at the specified position.
}
\description{
Searches ANHIG/IMGT-HLA alignments and returns protein, codons or nucleotide sequences for a submitted allele name and position(s).
}
\note{
This function requires that the HLAalignments object has been populated with alignments via the alignmentFull() function.

Indel positions must be text-formatted (e.g. "607.12").
}
\examples{
\dontrun{
alignmentSearch("gen","DRB1*15:07:01",11:22)
alignmentSearch("nuc","DRB1*15:07:01",11:22)
alignmentSearch("prot","DRB1*15:07:01",11:22)
alignmentSearch("codon","DRB1*15:07:01",11:22)
alignmentSearch("nuc","DRB1*11:250N",c(605,"607.1","607.12",608))
}

}
