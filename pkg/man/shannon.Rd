\name{shannon}
\docType{methods}
\alias{shannon}
\alias{shannon-methods}
\alias{shannon,phylo4com-method}
\alias{shannon,phylo4d-method}
\title{Shannon diversity index}
\description{
  Calculates the Shannon diversity index for a set of communities.
}
\usage{
  \S4method{shannon}{phylo4com}(x)
  \S4method{shannon}{phylo4d}(x)
}
\arguments{
  \item{x}{A \code{phylo4com} or \code{phylo4d} object.}
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
  shannon(weeds)
}
\keyword{methods}
