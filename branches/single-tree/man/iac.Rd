\name{iac}
\alias{iac}
\alias{iac-methods}
\alias{iac,phylo4com-method}
\alias{iac,phylo4d-method}
\title{Imbalance of Abundances at the Clade level}
\description{
  Calculates IAC for a set of communities with phylogenetic data.
}
\usage{
  \S4method{iac}{phylo4com}(x, na.rm = TRUE)
  \S4method{iac}{phylo4d}(x, na.rm = TRUE)
}
\arguments{
  \item{x}{A \code{phylo4d} object.}
  \item{na.rm}{logical. Should NA abundances be treated as zero?}
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
  iac(weeds)
}
