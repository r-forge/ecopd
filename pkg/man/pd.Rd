\name{pd}
\alias{pd}
\title{Calculate phylogenetic diversity}
\description{
  Calculates PD (phylogenetic diversity) for a phylogeny
}
\usage{
  pd(phy, method=c("traditional", "polytomy", "yule"))
}
\arguments{
  \item{phy}{A \code{phylo4d} object.}
  \item{method}{Method for calculating pd.}
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

  minTL(weeds$A) <- getMinTL(weeds$A, genera(weeds$A))
  pd(weeds$A, "polytomy")
  pd(weeds$A, "yule")
}
