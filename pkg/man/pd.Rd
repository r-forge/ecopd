\name{pd}
\docType{methods}
\alias{pd}
\alias{pd-methods}
\alias{pd,phylo4com-method}
\alias{pd,phylo4d-method}
\alias{pd,phylo4-method}
\title{Calculate phylogenetic diversity}
\description{
  Calculates PD (phylogenetic diversity) for a phylogeny
}
\usage{
  \S4method{pd}{phylo4com}(x, method = c("traditional", "polytomy",
    "yule"))
  \S4method{pd}{phylo4d}(x, method = c("traditional", "polytomy",
    "yule"), \dots)
  \S4method{pd}{phylo4}(x, method = c("traditional"))
}
\arguments{
  \item{x}{A \code{phylo4com}, \code{phylo4d}, or \code{phylo4} object.}
  \item{method}{Method for calculating pd.}
  \item{\dots}{additional arguments to be passed to other methods}
}
\details{
  TODO

  When applied to a \code{phylo4} tree, traditional PD is the only
  option (because the object contains no abundance data in this case).
}
\value{
  For \code{phylo4com} or \code{phylo4d}, a numeric vector of length
  equal to the number of communities. For \code{phylo4}, a single
  numeric value for the tree.
}
\references{
  TODO
}
\author{Jim Regetz (regetz@nceas.ucsb.edu)}
\examples{
  pd(weeds)

  # traditional pd of the master tree
  pd(phylo4(weeds))

#  minTL(weeds$A) <- getMinTL(weeds$A, genera(weeds$A))
#  pd(weeds$A, "polytomy")
#  pd(weeds$A, "yule")

}
