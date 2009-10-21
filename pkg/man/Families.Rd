\name{Families}
\alias{Families}
\docType{data}
\title{Plant family lookup vector}
\description{
  Built-in vector useful for determining the plant family for a
  particular plant genus.
}
\format{
  A named character vector matching 13,584 genus names with their
  respective families. Vector names are primarily genus names, although
  the 408 family names are also included for self-referential lookups.
}
\source{
  [reference?]
}
\examples{
Families["Liriodendron"]

# Matching on family "Asteraceae" simply returns itself
Families["Asteraceae"]
}
\keyword{datasets}
