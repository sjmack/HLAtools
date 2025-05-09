% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{IMGTHLAGeneTypes}
\alias{IMGTHLAGeneTypes}
\title{Molecular characteristics of the Genes Curated by the IPD-IMGT/HLA Database}
\format{
A list object of two elements:
\itemize{
\item(GeneTypes: a data frame of three columns)
\item(version: a character string identifying the date that the source file was written)
}
}
\source{
@source https://hla.alleles.org/pages/genes/genes_list
}
\usage{
data(IMGTHLAGeneTypes)
}
\description{
A data frame of three columns identifying each gene supported by the IPD-IMGT/HLA Database, its molecular characteristics, and its status as a gene fragment or pseudogene.
This object is built by the BuildIMGTHLAGeneTypes() function.
}
\keyword{datasets}
