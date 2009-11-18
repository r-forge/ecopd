\name{pairdist}
\alias{pairdist}
\title{Calculate pairwise distances between nodes}
\description{
  Returns a symmetric matrix containing the pairwise distances between
  all nodes, or possibly the pairwise distances between tips only.
}
\usage{
  pairdist(phy, type=c("all", "tip"), use.labels=FALSE)
}
\arguments{
  \item{phy}{A \code{phylo4} object (or one that inherits from it).}
  \item{type}{Should pairwise distances be returned for "all" nodes
    (default), or only for "tip" nodes?}
  \item{use.labels}{(logical) If FALSE, row and column names are the
    numeric node IDs, otherwise labels will be used (possibly resulting
    in \code{NA} names in the case of missing node labels).}
}
\details{
  TODO
}
\value{
  Numeric matrix of distances, with either node IDs or node labels as
  row and column names.
}
\author{Jim Regetz (regetz@nceas.ucsb.edu)}
\seealso{
  ape package functions \code{cophenetic} and \code{dist.nodes}, on
  which \code{pairdist} is based.
}
\examples{
  # all pairwise distances
  pairdist(weeds)

  # pairwise distances between tips, with labels and dim names
  pairdist(weeds, type="tip", use.labels=TRUE)
}
