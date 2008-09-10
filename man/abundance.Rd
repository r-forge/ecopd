\name{abundance}
\alias{abundance}
\alias{abundance<-}
\title{ Species abundances in a phylogeny }
\description{
Functions to get or set the names of an object.
}
\usage{
abundance(tree)
abundance(tree) <- value
}
\arguments{
  \item{tree}{ A phylogenetic tree (i.e. an object of class '"phylo"')}
  \item{value}{ Numeric vector of length equal to the number of tips in
tree. Order should correspond to the order of species given by
tree\$tip.label }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
     For 'abundance', 'NULL' or a named vector of the same length as
     tips of tree. ('NULL' is given if the tree has no assigned
     abundance values.) For 'abundance<-', the updated tree.
}
\author{ Jim Regetz }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(tree, value) {

  if (!is.numeric(value)) {
    stop("abundance values must be a numeric vector")
  } else if (length(value)!=length(tree$tip.label)) {
    stop("number of abundance values must equal number of species")
  }

  tree$abundance <- value
  return(tree)

  }
}
