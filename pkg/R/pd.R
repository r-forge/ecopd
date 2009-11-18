##
## PD methods
##

setGeneric("pd", function(x, ...) {
    standardGeneric("pd")
})

setMethod("pd", signature(x="phylo4"),
  function(x, method=c("traditional")) {
    method <- match.arg(method)
    if (isRooted(x)) {
        nonroot.nodes <- setdiff(nodeId(x), rootNode(x))
        tot.length <- sum(edgeLength(x, nonroot.nodes))
    } else {
        tot.length <- sum(edgeLength(x))
    }
    return(tot.length)
  }
)

setMethod("pd", signature(x="phylo4d"), function(x,
  method=c("traditional", "polytomy", "yule"), ...) {
    phyc <- phylo4com(x, ...)
    pd(phyc, method=method)
})

setMethod("pd", signature(x="phylo4com"), function(x,
  method=c("traditional", "polytomy", "yule")) {

    method <- match.arg(method)

    if (method %in% c("polytomy", "yule")) {
        message("estimating minTL values from Supertree")
        res <- sapply(communities(x), function(community) {
            phyd <- phylo4d(x, community)
            n <- tipData(phyd)[[community]]
            pd(extractTree(weightByAbund(phyd, n, method)))
        })
    } else if (method == "traditional") {
        comms <- x@metadata$comms
        if (is.null(comms)) {
            return(.pd(extractTree(x), method))
        }
        subtrees <- x@subtrees[as.character(comms)]
        res <- sapply(subtrees, pd)[as.character(comms)]
        names(res) <- names(comms)
    }
    return(res)

})

# lookup function for minimum tip length
getMinTL <- function(tree, genera) {

  Supertree <- phylo4(Supertree)

  if (missing(genera)) stop("must supply vector of genera")

  # Families, Supertree, and LookupTL are all system data built into
  # the package

  # Supplied tree must be ultrametric
# Not sure this is necessary?? Removing for now...
#  tree.phy <- suppressWarnings(as(foo, "phylo"))
#  if (!is.ultrametric(tree.phy)) {
#    stop(deparse(substitute(tree)), " is not ultrametric")
#  }

  # TODO: if genus not in Families vector, need to give warning.
  families <- Families[genera]

  # if Family is not in LookupTL give a warning, and later
  # use minTL based on average across all families in the lookup table
  if (!all(families %in% row.names(LookupTL))) {
    warning("one or more taxa missing from lookup table; using mean minTL")
  }

  familiesInSupertree <- LookupTL[families, "supertree"]
  familiesInSupertree <- as.character(familiesInSupertree)

  # Subset the supertree using families from the user-supplied tree. If
  # any user-supplied taxa cannot be matched to families in the
  # supertree, they are simply ignored
  subsupertree <- subset(Supertree, na.omit(familiesInSupertree))
  subsupertree.maxLength <- max(pairdist(subsupertree, type="tip"))/2
  tree.maxLength <- max(tipLength(tree, from="root"))

  tableTL <- LookupTL[familiesInSupertree, "minTL"]

  # if any TLs are 0 or NA, give warning and use average minTL across
  # all families in the LookupTL table
  if(any(tableTL==0 | is.na(tableTL))) {
    numNA <- sum(tableTL==0 | is.na(tableTL))
    warning("Using meanTL for ", numNA, " tip", if(numNA>1) "s") 
    # calculate average minTL across the entire LookupTL table,
    # excluding any non-positive or NA values
    meanMinTL <- mean(LookupTL$minTL[LookupTL$minBL > 0], na.rm=TRUE)
    tableTL[is.na(tableTL)] <- meanMinTL  
    tableTL[tableTL<=0] <- meanMinTL  
  }

  lookupTL <- tableTL * (tree.maxLength / subsupertree.maxLength)
  actualTL <- tipLength(tree, from="parent")
  minTL <- ifelse(lookupTL < actualTL, lookupTL, actualTL)

  return(minTL)
 
}

# function to weight tip length based on abundance
weightByAbund <- function(tree, n, method=c("polytomy", "yule")) {

  method <- match.arg(method)

  if (is.null(minTL(tree))) {
    minLength <- getMinTL(tree, genera(tree))
  } else {
    minLength <- minTL(tree)
  }

  tipLen <- tipLength(tree)

  if (method=="polytomy") {
    newLength <- tipLen + (n-1) * minLength
  } else if (method=="yule") {
    C <- 0.57722
    newLength <- tipLen - minLength +
      (minLength * (n-1) / (log(n) + C-1))
  }

  # Note that this is (unfortunately) tightly coupled to the
  # implementation of phylo4 objects -- an assignment method would be
  # better
  tree@edge.length[getEdge(tree, nodeId(tree, "tip"))] <- newLength

  return(tree)

}
