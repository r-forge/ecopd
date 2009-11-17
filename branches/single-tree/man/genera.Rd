\name{genera}
\alias{genera}
\title{Extract genus names}
\description{
  Helper function to extract the genus names from phylogeny tip labels
  in the case when those labels are scientific names using binomial
  nomenclature.
}
\usage{
  genera(phy, comm)
}
\arguments{
  \item{phy}{a \code{phylo4com} object}
  \item{comm}{a vector of community labels; if missing, all communities}
}
\details{
  This is a simple function that assumes the tip labels of phylogeny
  \code{phy} are of the form \dQuote{genus_species}. It does no checking
  whatsoever to check whether this is in fact the case.
}
\value{
  A list of character vectors, one per community, containing the genus
  portion of the tip labels, with element names corresponding to the
  numeric tip node ID.
}
\author{Jim Regetz (regetz@nceas.ucsb.edu)}
\examples{
  data(weeds)
  genera(weeds, "A")
  genera(weeds)
}
