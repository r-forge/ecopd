\name{ed}
\alias{ed}
\alias{ed-methods}
\alias{ed,phylo4d-method}
\alias{ed,phylo4com-method}
\alias{aed}
\alias{aed-methods}
\alias{aed,phylo4d-method}
\alias{aed,phylo4com-method}
\alias{eed}
\alias{eed-methods}
\alias{eed,phylo4d-method}
\alias{eed,phylo4com-method}
\alias{hed}
\alias{hed-methods}
\alias{hed,phylo4d-method}
\alias{hed,phylo4com-method}
\alias{eaed}
\alias{eaed-methods}
\alias{eaed,phylo4d-method}
\alias{eaed,phylo4com-method}
\alias{haed}
\alias{haed-methods}
\alias{haed,phylo4d-method}
\alias{haed,phylo4com-method}
\alias{value}
\alias{value-methods}
\alias{value,phylo4d-method}
\alias{value,phylo4com-method}
\title{Abundance-weighted evolutionary distinctiveness}
\description{
  Calculates AED (and related metrics) for a set of communities with
  phylogenetic data.
}
\usage{
  \S4method{aed}{phylo4com}(x, na.rm = TRUE)
  \S4method{aed}{phylo4d}(x, na.rm = TRUE)
}
\arguments{
  \item{x}{A \code{phylo4d} object.}
  \item{na.rm}{logical. Should NA abundances be treated as zero?}
}
\details{
  TODO
}
\value{
  For \code{aed}, \code{ed}, and \code{value}, a list of numeric
  vectors, one per community, with values for each species in the
  community. For \code{eed}, \code{hed}, \code{eaed}, and \code{haed}, a
  numeric vector of length equal to the number of communities.
}
\references{
  TODO
}
\author{Jim Regetz (regetz@nceas.ucsb.edu)}
\examples{
  aed(weeds)
  haed(weeds)
}
