\name{phylo4com-class}
\docType{class}
\alias{phylo4com-class}
\title{phylo4com class}
\description{S4 class for community phylogenies.}
\section{Objects from the Class}{
  Objects can be created from trees and community-specific abundance
  data using various forms of the \code{phylo4d} constructor, or using
  \code{new("phylo4com", \dots)} for empty objects.
}
\section{Slots}{
  \describe{
    \item{\code{edge}:}{Matrix of edges}
    \item{\code{edge.label}:}{Character vector of edge (branch) labels}
    \item{\code{edge.length}:}{Numeric vector of edge (branch) lengths}
    \item{\code{label}:}{Character vector of tip (and optionally
      internal node) labels}
    \item{\code{data}:}{data frame for traits of tips and internal nodes}
    \item{\code{subtrees}:}{list of subset trees corresponding to
      individual communities}
    \item{\code{order}:}{character: tree ordering}
    \item{\code{annote}:}{annotation data for tree (currently
      unstructured/unused by methods)}
    \item{\code{metadata}:}{metadata for node/tip data (currently
      used to store some community information)}
  }
}
\section{Extends}{
  Class \code{"phylo4d"} directly;
  Class \code{"phylo4"}, by class \code{"phylo4d"}.
}
\author{Jim Regetz}
\examples{
  showClass("phylo4com")
}
\keyword{classes}
