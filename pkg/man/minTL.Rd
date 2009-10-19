\name{minTL}
\alias{minTL<-}
\alias{minTL}
\title{ Species minTLs in a phylogeny }
\description{
Functions to get or set species minimum tip lengths (minTL)
}
\usage{
minTL(tree)
minTL(tree) <- value
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
     For 'minTL', 'NULL' or a named vector of the same length as
     tips of tree. ('NULL' is given if the tree has no assigned
     minTL values.) For 'minTL<-', the updated tree.
}
\author{ Jim Regetz }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
