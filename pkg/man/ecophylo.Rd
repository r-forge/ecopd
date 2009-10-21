\name{ecophylo}
\alias{phylo4com}
\title{Create list of community phylogenies}
\description{
  Function to combine species abundance data with phylogenetic
  information for a set of communities.
}
\usage{
  phylo4com(communityID, species, abundance, tree,
    missing=c("warn", "OK", "fail"))
}
\arguments{
  \item{communityID}{Vector of community identifiers.}
  \item{species}{Vector of species names; must correspond to tip labels
    of the supplied tree.}
  \item{abundance}{Numeric vector of abundances of each species in each
    community.}
  \item{tree}{A \code{phylo4} object (or one that inherits from it).}
  \item{missing}{(not currently implemented) What to do if some
    requested node IDs or names are not in the tree: warn, do nothing,
    or stop with an error.}
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
  A list of \code{phylo4d} objects, one per unique value in the
  \code{communityID} vector. List element names are taken from the
  \code{communityID} vector.
}
\author{Jim Regetz (regetz@nceas.ucsb.edu)}
\seealso{
  \code{\link[phylobase]{subset}} in the phylobase package
}
\examples{
  ##TODO -- contruction example
  # Example of a phylo4com result
  data(weeds)
  print(weeds)
}
