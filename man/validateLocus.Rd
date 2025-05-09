% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/locusValidator.R
\name{validateLocus}
\alias{validateLocus}
\title{Determine if a Locus Name is in the HLAgazeteer}
\usage{
validateLocus(loci, source)
}
\arguments{
\item{loci}{A character vector of HLA gene names (ex. "DRB1", c("DRB1","DQB1")).}

\item{source}{A character vector of alignment source types. "AA", "cDNA", and "gDNA" are allowed types.}
}
\value{
A logical value. TRUE indicates that all of the names and source types are valid. FALSE indicates that at least one locus name or alignment source type is invalid.
}
\description{
Checks a vector of HLA locus names against the HLAgazeteer to determine if the locus name is valid for a specific type of alignment.
}
\note{
The results of this check should only be considered valid for the version of the HLAgazeteer included in the HLAtools package.
}
\examples{
validateLocus(loci = "DRB1", source = "AA")
validateLocus(loci = c("V"), source = c("cDNA","gDNA"))
validateLocus(loci = c("E","F","G"), source = "gDNA")
}
