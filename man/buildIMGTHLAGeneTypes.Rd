% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extractGeneTypes.R
\name{buildIMGTHLAGeneTypes}
\alias{buildIMGTHLAGeneTypes}
\title{Builds an R object describing all of the genes supported by the IPD-IMGT/HLA Database, and identifying those that are pseudogenes or gene fragments.}
\usage{
buildIMGTHLAGeneTypes()
}
\value{
A list objects of two elements -- 'version' and 'GeneTypes'. The 'version' element identifies the date that the source table was generated. The 'GeneTypes' element is a data frame of three columns, identifying each gene supported by the IPD-IMGT/HLA Database, along with its molecular characteristics and its status as either a pseudogene or a gene fragment.
}
\description{
This function scrapes information from 'hla.alleles.org/genes/index.html' and generates a data frame identifying each gene. It requires internet access to function. As such, this function always returns data for the current IPD-IMGT/HLA Database release.
}
\note{
For internal HLAtools use.
}
\examples{
\dontrun{ 
IMGTHLAGeneTypes <- buildIMGTHLAGeneTypes()
}

}
