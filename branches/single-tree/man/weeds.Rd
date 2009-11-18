\name{weeds}
\alias{weeds}
\alias{weeds.data}
\alias{weeds.tree}
\docType{data}
\title{Weed phylogeny and community data}
\description{
  Phylogeny and abundance data for a small set of weedy plant species
  occurring in three hypothetical plots.
}
\format{
  \describe{
    \item{\code{weeds.data}:}{data frame with one column each for
      plot name, species name, and observed cover}
    \item{\code{weeds.tree}:}{phylo4 tree for all five species included
      in the community dataset}
    \item{\code{weeds}:}{phylo4com object representing the combined
      phylogeny and community data}
  }
}
\details{
  This data set includes five weedy annual plant species commonly found
  across North America. Abundance values represent fabricated data for
  three hypothetical field sites (named A, B, and C). The phylogeny
  component of the dataset is real, and was estimated from four genes:
  rbcL, MatK, ITS1, and 5.8s.
}
\references{
  Cadotte, M., Davies, J., Regetz, J., Kembel, S., Cleland, E., and
  Oakley, T. (In press) Phylogenetic diversity metrics for ecological
  communities: integrating species richness, abundance and evolutionary
  history. \emph{Ecology Letters}.
}
\keyword{datasets}
