% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/accessoryFunctions.R
\name{queryRelease}
\alias{queryRelease}
\title{Search Allele Names Across Release Versions}
\usage{
queryRelease(rel, variant = "", all = FALSE)
}
\arguments{
\item{rel}{An IPD-IMGT/HLA Database release version, represented as either a character (e.g., "3.56.0") or a numeric (e.g., 3560) value.}

\item{variant}{A character string. The value of 'variant' can be any part of a locus or allele name (e.g., "DR", "02:01", "DRB1*08:07"). The default ("") specifies all alleles in 'rel'.}

\item{all}{A logical. When 'all' = TRUE, a vector of all instances of 'variant' in 'rel' is returned.  When 'all' = FALSE, the number of instances of 'var' in 'rel' is returned.}
}
\value{
A character vector of all matches to 'variant' in 'rel' or the number of all matches to 'variant' in 'rel'.
}
\description{
Searches specific release versions in the AlleleListHistory object for user-defined allele variants.
}
\examples{
# Identify the number of DRB9 alleles in releases 3.30.0 and 3.31.0.
queryRelease("3.30.0","DRB9",FALSE) 
queryRelease("3.31.0","DRB9",FALSE)

# Identify the total number of alleles in release 3.56.0.
queryRelease(3560)

}
