\name{richness}
\alias{richness}
\title{Calculate richness for each community}
\description{
  Returns the species richness of each community, i.e. the number of tip
  taxa having positive abundance.
}
\usage{
richness(phylo4com, na.rm = FALSE)
}
\arguments{
  \item{phylo4com}{A list of phylo4d objects with abundance data}
  \item{na.rm}{logical. Should missing abundance values be treated as
    absences?}
}
\details{
  For each community tree, the upper bound on richness is simply the
  number of tips in the tree. However, the richness value will be lower
  if any species do not have positive abundance. When \code{na.rm=TRUE},
  species with \code{NA} abundance are considered missing, otherwise
  \code{NA} will propagate through to the result.
}
\value{
  Vector of richness values, with names corresponding the community
  names in the phylo4com list (if they exist).
}
\author{Jim Regetz}
\examples{
}
