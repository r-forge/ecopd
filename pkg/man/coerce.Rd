\name{coerce}
\docType{methods}
\alias{phylo4,phylo4com-method}
\alias{phylo4d,phylo4com-method}
\title{Extract phylo4com phylogenies}
\description{
  Methods to extract phylo4 and phylo4d trees from phylo4com objects.
}
\usage{
  \S4method{phylo4}{phylo4com}(x, community)
  \S4method{phylo4d}{phylo4com}(x, community)
}
\arguments{
  \item{x}{a \code{phylo4com} object}
  \item{community}{an optional single valued character vector indicating
    which community tree to return; if missing, the full
    (across-communities) tree is returned}
}
\details{
  TODO
}
\value{
  A \code{phylo4} or \code{phylo4d} object.
}
\author{Jim Regetz (regetz@nceas.ucsb.edu)}
\examples{
  phylo4(weeds)
  phylo4(weeds, "A")
  phylo4d(weeds, "A")
}
\keyword{methods}
