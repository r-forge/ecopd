\name{pae}
\alias{pae}
\alias{pae-methods}
\alias{pae,phylo4com-method}
\alias{pae,phylo4d-method}
\title{Phylogenetic-Abundance Evenness}
\description{
  Calculates PAE for a set of communities with phylogenetic data.
}
\usage{
  \S4method{pae}{phylo4com}(x, na.rm = TRUE)
  \S4method{pae}{phylo4d}(x, na.rm = TRUE)
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
  pae(weeds)
}
