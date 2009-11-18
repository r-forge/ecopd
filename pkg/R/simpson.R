##
## Simpson's index with and without phylogenetic distances
##

setGeneric("simpson",
  function(x, method=c("phylogenetic", "traditional")) {
    standardGeneric("simpson")
})

setMethod("simpson", signature(x="phylo4d"),
  function(x, method=c("phylogenetic", "traditional")) {
    phyc <- phylo4com(x)
    simpson(phyc, method=method)
})

setMethod("simpson", signature(x="phylo4com"),
  function(x, method=c("phylogenetic", "traditional")) {
    method <- match.arg(method)
    N.relative <- prop.table(t(abundance(x)), 1)
    if (method=="phylogenetic") {
        dmat <- pairdist(x, type="tip")
    } else {
        dmat <- matrix(1, nTips(x), nTips(x))
        diag(dmat) <- 0
    }   
    out <- apply(N.relative, 1, function(n) sum((n %o% n)*dmat))
    return(out) 
})

## earlier version: works on a single phylo4d tree with abundance data
#simpson <- function(phy, method=c("traditional", "phylogenetic")) {
#  method <- match.arg(method)
#  x <- prop.table(abundance(phy))
#  if (method=="phylogenetic") {
#    dmat <- pairdist(phy, type="tip")
#  } else {
#    dmat <- matrix(1, nTips(phy), nTips(phy))
#    diag(dmat) <- 0
#  }   
#  sum((x %o% x) * dmat)
#}
