# TODO: This function includes its own code for not counting root edge
# length. Maybe this should maybe be done at a higher level?
ED <- function(tree) {

  # set length of root edge to zero
  edgeLength(tree)[edgeId(tree, "root")] <- 0

  all.nodes <- nodeId(tree, type = "all")
  des <- descendants(tree, all.nodes, type="tips")
  nv <- edgeLength(tree, all.nodes) / sapply(des, length)
  names(nv) <- all.nodes

  tip.nodes <- nodeId(tree, "tip")
  anc <- ancestors(tree, tip.nodes, "ALL")
  EDI <- sapply(anc, function(n) sum(nv[as.character(n)], na.rm=TRUE))
  names(EDI) <- tipLabels(tree)

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

  # get all node IDs, but excluding root node
  nonroot.nodes <- setdiff(nodeId(tree), rootNode(tree))
  # Create logical matrix indicating which tips (in columns) are
  # descendants of each node (in rows), self-inclusive
  isDescendant <- sapply(ancestors(tree, tipLabels(tree), "ALL"),
    function(n) nonroot.nodes %in% n)

  # Create vector of ancestral edge lengths
  edge.length <- edgeLength(tree, nonroot.nodes)

  # Create matrix containing number of individuals of each species
  # descending from each interior node
  abundance <- abundance(tree)
  dAbund <- abundance * t(isDescendant)

  # Calculate individual-based AED of each species
  AED <- colSums(edge.length * t(prop.table(dAbund, margin=2)))
  AED <- AED/abundance

  return(AED)

}

HAED <- function(tree) {
  # Recast AED in terms of individuals
  AED <- AED(tree)
  abundance <- abundance(tree)
  scaledIndivAED <- rep(AED, abundance) / pd(tree)
  return(-sum(scaledIndivAED * log(scaledIndivAED)))
}

EAED <- function(tree) {
  HAED(tree) / log(sum(abundance(tree)))
}

value <- function(tree) {
  aed <- AED(tree)
  aed/sum(aed)
}

