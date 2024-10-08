% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkgDNAStart.R
\name{checkgDNAstart}
\alias{checkgDNAstart}
\title{Identify gDNA Alignments in Which the First Feature Boundary is not Identified as Position +1.}
\usage{
checkgDNAstart(verbose = FALSE)
}
\arguments{
\item{verbose}{A logical indicating if loci with first feature boundary positions that are not 1 should be reported in the console (verbose = TRUE).}
}
\value{
A one-row data frame with one column for each locus with a gDNA alignment.
}
\description{
Checks the position of the first feature boundary (the 5' UTR - Exon 1 boundary in expressed genes) in each gDNA alignment, and identifies those alignments in which that position is not 1.
}
\details{
This function reviews the gDNA atlases in the HLAatlas object and returns a data frame of the first feature boundary position for each locus. For expressed genes and some pseudogenes, this is the position of the start of Exon 1.
}
\note{
This function requires that the HLAalignments object has been populated with alignments via the alignmentFull() function.

For internal HLAtools use.
}
