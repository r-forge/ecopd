\name{richness}
\alias{richness}
\title{Calculate richness for each community}
\description{
  Returns the species richness of each community, i.e. the number of tip
  taxa having positive abundance.
}
\usage{
  richness(phy, comm, na.zero = FALSE)
}
\arguments{
  \item{phy}{a \code{phylo4com} object}
  \item{comm}{a vector of community labels; if missing, all communities}
  \item{na.zero}{logical. should NAs be treated as zero abundances?}
}
\details{
  For each community tree, the upper bound on richness is simply the
  number of tips in the tree. However, the richness value will be lower
  if any species do not have positive abundance. When
  \code{na.zero=TRUE}, species with \code{NA} abundance are considered
  missing, otherwise \code{NA} will propagate through to the result.
}
\value{
  Vector of richness values, with names corresponding to community
  names in the phylo4com object.
}
\author{Jim Regetz (regetz@nceas.ucsb.edu)}
\examples{
  data(weeds)
  richness(weeds)

  abundance(weeds, "A", 4) <- 0
  richness(weeds)

  abundance(weeds, "B", 2:3) <- NA
  richness(weeds, na.zero=TRUE)
}
