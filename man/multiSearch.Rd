% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SearchandBuildAlignments.R
\name{multiSearch}
\alias{multiSearch}
\title{Search Alignment Sequences at Multiple Positions for a Specified Allele}
\usage{
multiSearch(
  alignType,
  locus,
  allele,
  positions,
  prefix = TRUE,
  sep = "~",
  trimmed = FALSE
)
}
\arguments{
\item{alignType}{The type of alignment being searched. Allowed values are "prot", codon", "nuc" and "gen". Only one 'alignType' value is allowed.}

\item{locus}{A specific locus.}

\item{allele}{An allele name.}

\item{positions}{Nucleotide positions to search.}

\item{prefix}{A logical that indicates if the position number(s) should be included in the result.}

\item{sep}{Defines the separator between position sequences.}

\item{trimmed}{A logical identifying whether 'allele' is a two-field, 'trimmed' allele name (trimmed=TRUE), or a full-length name (trimmed=FALSE). The default value is FALSE.}
}
\value{
codons (possibly with position number) at 'positions' for 'allele'. If 'allele' does not have a codon at 'position', '' will be returned.
}
\description{
Generates a character string of multiple protein, codon or nucleotide or genomic sequences at the specified positions for the specified allele name.
}
\note{
For internal HLAtools use only.
}
