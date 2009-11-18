\name{pd}
\docType{methods}
\alias{pd}
\alias{pd-methods}
\alias{pd,phylo4com-method}
\alias{pd,phylo4d-method}
\title{Calculate phylogenetic diversity}
\description{
  Calculates PD (phylogenetic diversity) for a phylogeny
}
\usage{
  \S4method{pd}{phylo4com}(x, method = c("traditional", "polytomy",
    "yule"))
  \S4method{pd}{phylo4d}(x, method = c("traditional", "polytomy",
    "yule"), \dots)
}
\arguments{
  \item{x}{A \code{phylo4com} or \code{phylo4d} object.}
  \item{method}{Method for calculating pd.}
  \item{\dots}{additional arguments to be passed to other methods}
}
\details{
  TODO
}
\value{
  A numeric vector of length equal to the number of communities.
}
\references{
  TODO
}
\author{Jim Regetz (regetz@nceas.ucsb.edu)}
\examples{
  pd(weeds)
#  minTL(weeds$A) <- getMinTL(weeds$A, genera(weeds$A))
#  pd(weeds$A, "polytomy")
#  pd(weeds$A, "yule")
}
