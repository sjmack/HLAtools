% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BDtoPyPop.R
\name{formatHead}
\alias{formatHead}
\title{Format PyPop Data Frame Headers}
\usage{
formatHead(colHead)
}
\arguments{
\item{colHead}{A vector of column names. PyPop format requires that paired locus/gene names should end in '_1' and '_2', respectively.}
}
\value{
A vector in which the locus names are suffixed with '_1' and '_2'.
}
\description{
Format the header of a PyPop-formatted data frame.
}
\note{
This function assumes that the first two elements the 'colHead' vector are not locus/gene names.
}
\examples{
formatHead(colHead = colnames(sHLAdata))

}
\references{
Lancaster et al. Front Immunol. 2024 Apr 2;15:1378512. https://pubmed.ncbi.nlm.nih.gov/38629078/

Pappas et al. Hum Immunol. 2016 Mar 77(3):283-287. https://pubmed.ncbi.nlm.nih.gov/26708359/
}
