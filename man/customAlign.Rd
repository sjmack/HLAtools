% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SearchandBuildAlignments.R
\name{customAlign}
\alias{customAlign}
\title{Generate a customized peptide, codon or coding or genomic nucleotide sequence alignment.}
\usage{
customAlign(alignType, alleles, positions)
}
\arguments{
\item{alignType}{The type of alignment being searched. Allowed values are "prot", codon", "nuc" and "gen".}

\item{alleles}{A vector of un-prefixed HLA locus names.}

\item{positions}{Either a vector of variant positions, against which all loci will be aligned, or a list of vectors of nucleotide positions, exactly one vector for each allele, against which each corresponding allele will be aligned.}
}
\value{
A data frame of allele names and the corresponding nucleotide sequences for each desired nucleotide position. an error message is returned if input loci is not available in the ANHIG/IMGTHLA Github Repository.
}
\description{
Generates a peptide, codon, coding or genomic nucleotide alignment table for user-specified HLA alleles at user-specified positions.
}
\examples{
\dontrun{
customAlign("codon",c("DRB1*01:01","DQB1*02:01","DPB1*01:01"),c(1,2,3,7,8,9,13,14,15))
customAlign("codon",c("DPB1*01:01:01:01","DQA1*01:01:01:01","DQB1*05:01:01:01"),list(19:35,1:4,6:9))
}
}
