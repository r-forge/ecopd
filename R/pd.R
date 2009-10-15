# basic pd calculation
pd <- function(tree) {
  # exclude root edge from calculation (if it exists)
  if (isRooted(tree)) {
    nonroot.nodes <- setdiff(nodeId(tree), rootNode(tree))
    tot.length <- sum(edgeLength(tree, nonroot.nodes))
  } else {
    tot.length <- sum(edgeLength(tree))
  }
  return(tot.length)
}

# lookup function for minimum tip length
getMinTL <- function(tree, genera) {

  if (missing(genera)) stop("must supply vector of genera")

  lengthsTipToRoot <- function(tree) {
    # Workaround problem in phylobase -- need to figure out if this is
    # really the right thing to do...
    if (!is.na(tree@root.edge)) {
      warning("omitting root edge")
      is.na(tree@root.edge) <- TRUE
    }
    sapply(tipLabels(tree), function(x) {
      sumEdgeLength(tree, ancestors(tree, x, "ALL"))
    })
  }

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
  subsupertree.maxLength <- max(cophenetic(subsupertree))/2
  tree.maxLength <- max(lengthsTipToRoot(tree))

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
  actualTL <- tipLength(tree)
  minTL <- ifelse(lookupTL < actualTL, lookupTL, actualTL)

  return(minTL)
 
}

# function to weight tip length based on abundance
weightByAbund <- function(tree, method=c("polytomy", "yule")) {

  method <- match.arg(method)

  if (is.null(abundance(tree))) {
    stop("tree contains no abundance information")
  }

  if (is.null(minTL(tree))) {
    stop("tree contains no minTL information")
  }

  n <- abundance(tree)
  minLength <- minTL(tree)

  tipLen <- ancestralEdgeLength(foo, getnodes(foo,
    names(abundance(foo))))
  # Longer form above ensures consistent ordering, but given phylo4
  # rules it might be sufficient to use:
  # tipLen <- ancestralEdgeLength(foo, nodes(foo, "tip"))
 
  # Test statement:
  if (!identical(names(n), names(minLength))) {
    stop("mismatch between abundance and minTL vectors")
  }

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
  tree@edge.length[match(nodes(tree, "tip"), edges(tree)[,2])] <-
    newLength

  return(tree)

}

# calculate PD value for a single community
commPD <- function(tree, method=c("traditional", "polytomy", "yule")) {
  method <- match.arg(method)
  if (method == "traditional") {
    pd(tree)
  } else {
    pd(weightByAbund(tree, method))
  } 
}
