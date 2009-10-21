\name{abundance}
\alias{abundance}
\alias{abundance<-}
\title{Get or set species abundances}
\description{
Functions to get or set the abundances of tip taxa in a phylogenetic
tree.
}
\usage{
  abundance(phy)
  abundance(phy) <- value
}
\arguments{
  \item{phy}{a \code{phylo4d} object}
  \item{value}{Numeric vector of length equal to the number of tips in
    the tree. Order should correspond to the output of nodeId(tree,
    type="tip").}
}
\details{
  Any missing abundances will be represented by \code{NA}, including the
  case when abundance data are missing from the tree altogether. All
  abundance values can be removed by using an assignment value of
  \code{NULL}.

  Note that in the assignment form, subscripting can be used to replace
  particular abundance values; see examples.
}
\value{
  For \code{abundance}, a named vector of the same length as tips of
  tree. For \code{abundance<-}, the updated tree.
}
\author{Jim Regetz (regetz@nceas.ucsb.edu)}
\examples{
  data(weeds)
  abundance(weeds$A)

  ## modify one particular abundance value
  abundance(weeds$A)["Centaurea_alba"] <- 10
  abundance(weeds$A)

  ## remove abundance data altogether
  abundance(weeds$A) <- NULL
  abundance(weeds$A)
}
