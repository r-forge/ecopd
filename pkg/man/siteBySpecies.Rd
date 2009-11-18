\name{siteBySpecies}
\alias{siteBySpecies}
\title{Generate site-by-species matrix}
\description{
  Returns a site-by-species matrix of abundance or presence/absence
  values from a phylo4com object.
}
\usage{
  siteBySpecies(phy, presence.only=FALSE, na.zero=FALSE, transpose=FALSE)
}
\arguments{
  \item{phy}{a \code{phylo4com} object}
  \item{presence.only}{logical. Should abundance values be replaced
    with 1/0 indicators of presence/absence?}
  \item{na.zero}{logical. should NAs be treated as zero abundances?}
  \item{transpose}{logical. Should matrix be transposed to put 
    species in rows and sites in columns?}
}
\details{
  If \code{na.zero} is false (the default), any \code{NA} values in the
  abundance data will be preserved in the matrix. However, abundance is
  assumed to be 0 for all species absent from the subtree for a given
  site.

  If \code{presence.only} is true, the matrix will contain 1 for species
  with positive abundances in the corresponding tree, and 0 both for
  species with non-positive abundances and for species absent from the
  tree for a given site. As with abundance, if \code{na.zero} is false,
  \code{NA} values in the abundance data will remain as \code{NA} in the
  matrix.

  In the default orientation, species will be in columns and sites in
  rows. To swap rows and columns, use \code{transpose=TRUE} (or manually
  transpose the returned object).
}
\seealso{
  This function is just a convenience wrapper for the
  \code{\link{abundance}} and \code{\link{presence}} functions.
}
\value{
  Numeric matrix 
}
\author{Jim Regetz (regetz@nceas.ucsb.edu)}
\examples{
  data(weeds)
  siteBySpecies(weeds)
  siteBySpecies(weeds, presence.only=TRUE)

  ## change an abundance value to NA
  abundance(weeds, "A", 1) <- NA
  siteBySpecies(weeds, presence.only=TRUE, transpose=TRUE)
}
