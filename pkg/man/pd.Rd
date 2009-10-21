\name{pd}
\alias{pd}
\title{Calculate phylogenetic diversity}
\description{
  Calculates PD (phylogenetic diversity) for a phylogeny
}
\usage{
  pd(phy)
}
\arguments{
  \item{phy}{a \code{phylo4} object (or one that inherits from it)}
}
\details{
  TODO
}
\value{
  A numeric vector of length 1.
}
\references{
  TODO
}
\author{Jim Regetz (regetz@nceas.ucsb.edu)}
\examples{
  data(weeds)
  pd(weeds$A)
}
