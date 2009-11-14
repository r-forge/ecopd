\name{abundance}
\alias{abundance}
\alias{abundance<-}
\title{Get or set species abundances}
\description{
Functions to get or set the abundances of tip taxa in a phylogenetic
tree.
}
\usage{
  abundance(phy, comm, tip, na.zero=TRUE)
  presence(phy, comm, tip, na.zero=TRUE)
  abundance(phy, comm, tip) <- value
}
\arguments{
  \item{phy}{a \code{phylo4com} object}
  \item{comm}{a vector of community labels; if missing, all communities
    (not allowed in the replacement form)}
  \item{tip}{a vector of tip labels or numeric node IDs; if missing, all
    tips}
  \item{na.zero}{logical. should NAs be treated as zero abundances?}
  \item{value}{a numeric vector of abundance values}
}
\details{
  Any missing abundances will be represented by \code{NA}, unless
  \code{na.zero} is true, in which case they will be treated as zero.
  This is also true of the corresponding presence/absence matrix.

  The \code{tip} vector can either be character, in which case tip
  labels are matched, or integer, in which case the numeric node IDs of
  tips in the phylo4d tree will be matched. Users should not assume that
  these integer values will correspond to row positions in the tipData
  data frame, as this is an implementation detail that may change at any
  time.

  In the assignment form, \code{comm} must be a \emph{single} community
  label; it is not possible to update the abundances of multiple
  communities simultaneously. In contrast, \code{tip} can be a single
  value, multiple values, or missing (in which case all tip abundances
  will be replaced). If the \code{value} vector is shorter than the
  number of tips, normal R recycling rules will apply.
}
\value{
  For \code{abundance} and \code{presence}, a data frame.
  For \code{abundance<-}, the updated tree.
}
\author{Jim Regetz (regetz@nceas.ucsb.edu)}
\examples{
  data(weeds)
  abundance(weeds)
  abundance(weeds, "A")

  ## abundance replacement
  abundance(weeds, "A", "Centaurea_alba") <- NA
  abundance(weeds, "B", 2:5) <- c(0, 100, 0)
  abundance(weeds, "C") <- 25

  presence(weeds, na.zero=TRUE)
}
