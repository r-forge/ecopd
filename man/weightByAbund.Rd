\name{weightByAbund}
\alias{weightByAbund}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
weightByAbund(tree, method = c("polytomy", "yule"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{ ~~Describe \code{tree} here~~ }
  \item{method}{ ~~Describe \code{method} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(tree, minLength,
  method=c("polytomy", "yule")) {

###Jim, why the match.arg?
  method <- match.arg(method)

  if (is.null(abundance(tree))) {
    stop("tree contains no abundance information")
  }

  tipRowID <- which(tree$edge[, 2] <= Ntip(tree))
  tipLen <- tree$edge.length[tipRowID]
  tipID <- tree$edge[tipRowID,2]

  n <- abundance(tree)[tipID]
 
  # Test statement:
  if (!identical(names(n), names(minLength))) {
    stop("mismatch between abundance and minTL vectors")
  }

  if (method=="polytomy") {
    tree$edge.length[tipRowID] <- tipLen + (n-1) * minLength
  } else if (method=="yule") {
    C <- 0.57722
    tree$edge.length[tipRowID] <- tipLen - minLength +
      (minLength * (n-1) / (log(n) + C-1))
  }

  return(tree)

  }
}
