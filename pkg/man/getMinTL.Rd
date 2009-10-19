\name{getMinTL}
\alias{getMinTL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getMinTL(tree, genera)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{ ~~Describe \code{tree} here~~ }
  \item{genera}{ ~~Describe \code{families} here~~ }
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
function(tree, families=NULL, supertree, lookup) {

  # TODO: if genus not in FamGen lookup table, need to give warning. but
  # right now the lookup happens before this function gets called

  # if Family is not in the MinTL lookup table give a warning, and later
  # use minTL based on average across all families in the lookup table
  if (!all(families \%in\% row.names(lookup))) {
    warning("one or more taxa missing from lookup table; using mean minTL")
  }

  familiesInSupertree <- lookup[families, "Davies.phylogeny"]
  familiesInSupertree <- as.character(familiesInSupertree)

  # Subset the supertree using families from the user-supplied tree. If
  # any user-supplied taxa cannot be matched to families in the
  # supertree, they are simply ignored
  subsupertree <- subset(supertree, na.omit(familiesInSupertree))
  subsupertree.maxLength <- max(cophenetic(subsupertree))/2
  tree.maxLength <- max(cophenetic(tree))/2

  tableTL <- lookup[familiesInSupertree, "MinBL"]

  # if any TLs are 0 or NA, give warning and use average minTL across
  # all families in the lookup_MinTL table
  if(any(tableTL==0 | is.na(tableTL))) {
    numNA <- sum(tableTL==0 | is.na(tableTL))
    warning("Using meanTL for ", numNA, " tip", if(numNA>1) "s") 
    # calculate average minTL across the entire minTL lookup table,
    # excluding any non-positive or NA values
    meanMinTL <- mean(lookup[["MinBL"]][lookup[["MinBL"]]>0], na.rm=TRUE)
    tableTL[is.na(tableTL)] <- meanMinTL  
    tableTL[tableTL<=0] <- meanMinTL  
  }

  lookupTL <- tableTL * (tree.maxLength / subsupertree.maxLength)
  actualTL <- tipLengths(tree)
  minTL <- ifelse(lookupTL < actualTL, lookupTL, actualTL)

  return(minTL)
 
  }
}
