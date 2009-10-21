\name{genera}
\alias{genera}
\title{Extract genus names}
\description{
  Helper function to extract the genus names from phylogeny tip labels
  in the case when those labels are scientific names using binomial
  nomenclature.
}
\usage{
  genera(phy)
}
\arguments{
  \item{phy}{a \code{phylo4} object (or one that inherits from it)}
}
\details{
  This is a simple function that assumes the tip labels of phylogeny
  \code{phy} are of the form \dQuote{genus_species}. It does no checking
  whatsoever to check whether this is in fact the case.
}
\value{
  A character vector containing the genus portion of the tip labels,
  with element names corresponding to the numeric tip node ID.
}
\author{Jim Regetz (regetz@nceas.ucsb.edu)}
\examples{
  data(weeds)
  genera(weeds$A)
}
