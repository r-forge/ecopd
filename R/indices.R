PAE <- function(tree, abundance) {

  # Calculate PD
  PD <- pd(tree)

  # Calculate lengths of terminal branches
  TL <- tipLength(tree)

  # Match abundances to tips based on terminal labels (won't be
  # necessary if we incorporate abundances into the trait data frame
  # properly)
  #N <- abundance[names(TL)]
# New approach, assumes 'abundance' column in phylo4d data
  N <- tdata(tree)$abundance

  # Calculate PAE
  numer <- PD + sum(TL * (N-1))
  denom <- PD + (mean(N)-1) * sum(TL)
  return(numer/denom)

}

IAC <- function(tree, abundance) {

  # Workaround for problem previously in phylobase -- need to figure out
  # if this has actually been "fixed"
  #if (!is.na(tree@root.edge)) {
  #  warning("omitting root edge")
  #  is.na(tree@root.edge) <- TRUE
  #}

  # Count number of lineages originating at each internal node (i.e.
  # number of splits)
  nodeNames <- as.vector(nodeLabels(tree))
  nSplits <- sapply(nodeNames, function(x) length(children(tree,
    x)))

  # For each tip, take the product of the number of splits across all of
  # its ancestral nodes
  denom <- sapply(phylobase::labels(tree), function(x)
    prod(nSplits[names(ancestors(tree, x))]))

# This next line assumes 'abundance' column in phylo4d data
  abundance <- tdata(tree)$abundance

  # Calculate expected number of individuals under null hypothesis of
  # equal allocation to each lineage at each (node) split 
  expected <- sum(abundance) / denom

  # IAC: summed absolute difference between expected and observed
  # abundances, divide by number of nodes
  sum(abs(expected-abundance)) / nNodes(tree)

}

ED <- function(tree) {
# This function includes its own code for not counting root edge length.
# Maybe this should maybe be done at a higher level?
  nodes <- nodes(tree, which = "all", include.root = TRUE)
  nv <- sapply(nodes, function(node) {
    if (node != rootNode(tree)) {
      S <- length(descendants(tree, node, which = "tip"))
      ans <- as.vector(ancestralEdgeLength(tree, node)/S)
      return(ans)
    } else {
      return(0.0)
    }
  })
  tip.nodes <- nodes(tree, which = "tip")
  EDI <- sapply(tip.nodes, function(n) sum(nv[names(ancestors(tree, n,
    "ALL"))], na.rm = TRUE))
  return(EDI)
}

HED <- function(tree) {
  scaledED <- ED(tree) / pd(tree)
  return(-sum(scaledED * log(scaledED)))
}

EED <- function(tree) {
  HED(tree) / log(nTips(tree))
}

AED <- function(tree) {
 
  # Workaround problem in phylobase -- need to figure out if this is
  # really the right thing to do...
  #if (!is.na(tree@root.edge)) {
  #  warning("omitting root edge")
  #  is.na(tree@root.edge) <- TRUE
  #}

  # Create logical matrix indicating which tips (in columns) are
  # descendants of each interior node (in rows)
  #isDescendant <- sapply(phylobase::labels(tree), function(n) nodeLabels(tree) %in%
  #  names(ancestors(tree, n)))

  # Row corresponding to the root node will have all TRUE values. We
  # want to drop this node from the matrix.
  #isDescendant <- isDescendant[as.logical(rowSums(!isDescendant)),]

  # New way to do the above, using my node-finding function
#  int.nodes <- nodesInternal(tree, include.root=FALSE)
  int.nodes <- nodes(tree, which="int", include.root=FALSE)
  isDescendant <- sapply(phylobase::labels(tree), function(n) int.nodes %in%
    ancestors(tree, n))

  # Create vector of edges on ancestral side of each interior node
  # (having already excluded the root node). Currently assuming that
  # these edges are now in the same order as the nodes (rows) of the
  # above matrix, which is quite possibly wrong
  #intEdgeLengths <- edgeLength(tree)[-seq_len(nTips(tree))]
  #tipEdgeLengths <- edgeLength(tree)[seq_len(nTips(tree))]
  intEdgeLengths <- ancestralEdgeLength(tree, int.nodes)
  tipEdgeLengths <- tipLength(tree)

  # Create matrix containing number of individuals of each species
  # descending from each interior node
# This next line assumes 'abundance' column in phylo4d data
abundance <- tdata(tree)$abundance
  dAbund <- apply(isDescendant, 1, function(x) ifelse(x, abundance[x],
    0))

  # Calculate individual-based AED of each species
  AED <- tipEdgeLengths + colSums(intEdgeLengths * t(apply(dAbund, 2,
    function(x) x/sum(x))))

  AED <- AED/abundance

  return(AED)

}

HAED <- function(tree) {
  # Recast AED in terms of individuals
  AED <- AED(tree)
  abundance <- tdata(tree)$abundance
  scaledIndivAED <- rep(AED, abundance) / pd(tree)
  return(-sum(scaledIndivAED * log(scaledIndivAED)))
}

EAED <- function(tree) {
  HAED(tree) / log(sum(tdata(tree)$abundance))
}

value <- function(tree) {
  aed <- AED(tree)
  aed/sum(aed)
}

