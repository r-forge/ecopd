\name{subset.phylo}
\alias{subset.phylo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
subset.phylo(x, tipNames, checkMissing = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object of class \code{"phylo"} }
  \item{tipNames}{ ~~Describe \code{tipNames} here~~ }
  \item{checkMissing}{ ~~Describe \code{checkMissing} here~~ }
  \item{\dots}{ Additional variables }
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
subset(Supertree, c("Acanthaceae", "Zygophyllaceae"))

# Unmatched tipNames will generate an error...
\dontrun{subset(Supertree, c("Acanthaceae", "Zygophyllaceae",
"unknown"))}

# Unless checkMissing is FALSE
subset(Supertree, c("Acanthaceae", "Zygophyllaceae", "unknown"),
checkMissing=FALSE)
}
