\name{minTL}
\alias{minTL}
\alias{minTL<-}
\title{Get or set minimum tip lengths}
\description{
  Functions to get or set species minimum tip lengths (minTL) for a
  phylogeny.
}
\usage{
  minTL(phy)
  minTL(phy) <- value
}
\arguments{
  \item{phy}{a \code{phylo4d} object}
  \item{value}{Numeric vector of length equal to the number of tips in
    the tree. Order should correspond to the output of nodeId(tree,
    type="tip").}
}
\details{
  TODO
}
\value{
  For \code{minTL}, a named numeric vector of the same length as tips of
  tree. (\code{NULL} is returned if the tree has no assigned minTL
  values.) For \code{minTL<-}, the updated tree.
}
\author{Jim Regetz (regetz@nceas.ucsb.edu)}
\examples{
  data(weeds)
  minTL(weeds$A)
}
