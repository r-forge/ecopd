\name{phylo4com}
\docType{methods}
\alias{phylo4com}
\alias{phylo4com-methods}
\alias{phylo4com,phylo4,numeric-method}
\alias{phylo4com,phylo4,matrix-method}
\alias{phylo4com,phylo4,data.frame-method}
\alias{phylo4com,phylo4d,missing-method}
\alias{phylo4com,phylo,ANY-method}
\title{Create community phylogenies object}
\description{
  Function to combine species abundance data with phylogenetic
  information for a set of communities.
}
\section{Methods}{
  \describe{
    \item{x = "phylo4", n = "numeric"}{creates community phylogenies
      from a master \code{phylo4} tree, a vector of abundances, a vector
      of corresponding community names/labels, and a vector of
      corresponding taxa names/labels}
    \item{x = "phylo4", n = "matrix"}{creates community phylogenies
      from a master \code{phylo4} tree and a site-by-species abundance
      matrix.}
    \item{x = "phylo4", n = "data.frame"}{creates community phylogenies
      from a master \code{phylo4} tree and a site-by-species abundance
      data.frame}
    \item{x = "phylo", n = "ANY"}{creates a phylobase tree from an
      object of class \code{phylo}}
  }
}
\usage{
  \S4method{phylo4com}{phylo4,numeric}(x, n, communityID, species,
    missing=c("warn", "OK", "fail"))
  \S4method{phylo4com}{phylo4,matrix}(x, n, missing=c("warn", "OK", "fail"))
  \S4method{phylo4com}{phylo4,data.frame}(x, n, missing=c("warn", "OK", "fail"))
  \S4method{phylo4com}{phylo4d,missing}(x, n, cols)
  \S4method{phylo4com}{phylo,ANY}(x, n, \dots, check.node.labels = "keep")
}
\arguments{
  \item{x}{a phylogeny of class \code{phylo4d}, \code{phylo4}, or
    \code{phylo}}
  \item{n}{Vector of abundances of each species in each community, or
    abundance matrix with species as rows and sites as columns.}
  \item{communityID}{Vector of community labels corresponding to
    the abundance values (only used if \code{n} is a vector).}
  \item{species}{Vector of species labels corresponding to
    the abundance values (only used if \code{n} is a vector).}
  \item{cols}{Names of pre-existing phylo4d data columns to use as
    abundance data (only if \code{n} is missing)}
  \item{check.node.labels}{If x is a phylo object, what to do with node
    labels}
  \item{missing}{(not currently implemented) What to do if some
    requested node IDs or names are not in the tree: warn, do nothing,
    or stop with an error.}
  \item{\dots}{additional arguments to be passed to other methods}
}
\details{
  TODO

  If the species vector contains values that do not appear in the tree,
  they are ignored with a warning. If a community has fewer than one
  species matched in the tree, an error will result. In both cases, the
  behavior is simply a result of how the phylobase \code{subset} method
  works; see the help page for this method for more information.
}
\value{
  A \code{phylo4com} objects.
}
\author{Jim Regetz (regetz@nceas.ucsb.edu)}
\seealso{
  \code{\link[phylobase:subset-methods]{subset}} in the phylobase package
}
\examples{
  # Example of phylo4com construction
  data(weeds.data)
  data(weeds.tree)
  phylo4com(weeds.tree, weeds.data$cover, weeds.data$plot,
    weeds.data$taxa)
}
\keyword{methods}
