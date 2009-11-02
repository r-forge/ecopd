\name{siteBySpecies}
\alias{siteBySpecies}
\title{Generate site-by-species matrix}
\description{
  Returns a site-by-species matrix of abundance values from a phylo4com
  list.
}
\usage{
siteBySpecies(phylo4com, presence.only=FALSE, transpose=FALSE)
}
\arguments{
  \item{phylo4com}{A list of phylo4d objects with abundance data}
  \item{presence.only}{(logical) Should abundance values be replaced
    with 1/0 indicators of presence/absence? Default FALSE.}
  \item{transpose}{(logical) Should matrix be transposed to put 
    species in rows and sites in columns? Default FALSE.}
}
\details{
  Abundance values are taken from the \dQuote{abundance} data column in
  each \code{phylo4d} tree. Any \code{NA} values in the abundance data
  will be preserved in the matrix. However, abundance is assumed to be 0
  for all species absent from the \code{phylo4d} tree for a given site.

  If \code{presence.only=TRUE}, the matrix will contain 1 for species
  with positive abundances in the corresponding tree, and 0 both for
  species with non-positive abundances and for species absent from the
  tree for a given site. Any \code{NA} values in the abundance data will
  remain as \code{NA} in the matrix.

  Species names are taken directly from tip labels in the \code{phylo4d}
  trees, and the total number of species in the matrix will equal the
  number of unique tip labels across all trees. Site names are taken
  directly from the names of elements in the \code{phylo4com} list, and
  the number of sites is simply equal to the number of \code{phylo4d}
  objects. 

  In the default orientation, species will be in columns and sites in
  rows. To swap rows and columns, use \code{transpose=TRUE} (or manually
  transpose the returned object).
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
  abundance(weeds[[1]])[1] <- NA
  siteBySpecies(weeds, presence.only=TRUE, transpose=TRUE)
}
