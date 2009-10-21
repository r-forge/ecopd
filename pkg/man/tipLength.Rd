\name{tipLength}
\alias{tipLength}
\title{Extract terminal branch lengths}
\description{
  Returns the lengths of all tip edges of a phylogeny.
}
\usage{
  tipLength(phy)
}
\arguments{
  \item{phy}{a \code{phylo4} object (or one that inherits from it)}
}
\value{
  Numeric vector of lengths, with the associated tip node labels as
  names.
}
\author{Jim Regetz (regetz@nceas.ucsb.edu)}
\seealso{
  This function is a simple wrapper of the
  \code{\link[phylobase:phylo4-accessors]{edgeLength}} phylobase method.
}
\examples{
  data(weeds)
  tipLength(weeds$A)
}
