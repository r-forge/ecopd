\name{tipLength}
\alias{tipLength}
\title{Calculate distance to tips from ancestors}
\description{
  Returns the edge length between each tip and its parent, or
  alternatively the total distance to each tip from the root node of the
  phylogeny.
}
\usage{
  tipLength(phy, from=c("parent", "root"))
}
\arguments{
  \item{phy}{a \code{phylo4} object (or one that inherits from it)}
  \item{from}{Should lengths represent distance from parents only (i.e.,
    simply the lengths of terminal branches), or from the root node?}
}
\value{
  Numeric vector of lengths, with the associated tip node labels as
  names.
}
\note{
  When calculating the distance from root, the length of the root edge
  itself is \strong{not} included in the returned value. See examples
  for a simple way to add this additional edge to the result.
}
\author{Jim Regetz (regetz@nceas.ucsb.edu)}
\seealso{
  Methods \code{\link[phylobase:phylo4-accessors]{edgeLength}} and
  \code{\link[phylobase:phylo4-accessors]{sumEdgeLength}} in the
  phylobase package.
}
\examples{
  data(weeds)
  tipLength(weeds$A)

  # distance from the root node to each tip
  tipLength(weeds$A, from="root")

  # as above, but include root edge length
  tipLength(weeds$A, from="root") + edgeLength(weeds$A,
    rootNode(weeds$A))
}
