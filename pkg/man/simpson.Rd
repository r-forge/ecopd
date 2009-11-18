\name{simpson}
\docType{methods}
\alias{simpson}
\alias{simpson-methods}
\alias{simpson,phylo4com-method}
\alias{simpson,phylo4d-method}
\title{Phylogenetic Simpson's diversity index}
\description{
  Calculates Simpson's diversity index for a set of communities,
  possibly in a way that incorporates phylogenetic information
}
\usage{
  \S4method{simpson}{phylo4com}(x, method = c("phylogenetic",
    "traditional"))
  \S4method{simpson}{phylo4d}(x, method = c("phylogenetic",
    "traditional"))
}
\arguments{
  \item{x}{A \code{phylo4com} or \code{phylo4d} object.}
  \item{method}{character indicating whether to compute the
    \dQuote{phylogenetic} or \dQuote{traditional} Simpson's index}
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
  simpson(weeds)
}
\keyword{methods}
