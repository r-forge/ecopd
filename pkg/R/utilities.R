nodes <- function(tree, which=c("all", "internal", "tip"),
  include.root=TRUE) {

  which <- match.arg(which)

  x <- as(tree, "data.frame")
  node.type <- x$node.type
  node <- x$node
  names(node) <- x$label

  tip.nodes <- node[node.type=="tip"]
  if (which=="tip") {
    return(tip.nodes)
  }

  int.nodes <- node[node.type=="internal"]
  if (include.root) {
    root <- node[node.type=="root"]
    int.nodes <- c(root, int.nodes)
  }

  if (which=="internal") {
    return(int.nodes)
  } else if (which=="all") {
    return(c(int.nodes, tip.nodes))
  }

  stop("internal error, contact package maintainer")

}

#
# TODO: Make sure nothing funny is going on here wrt the order of nodes
# as reported
#
ancestralEdgeLength <- function(tree, node=NULL) {
  E <- edges(tree)
  ancestor <- E[, 1]
  if (is.null(node)) {
    node <- E[,2]
  } else if (!all(node %in% E)) {
    stop("one or more nodes not found in tree")
  }
  ## retrieve the ancestor of each node
  idx <- match(node, E[, 2]) # new ordering of the descendants/edges
  ## if (length(ancestor)>0) ancestor <- c(NA, ancestor)
  ancestor <- E[idx, 1]
  ## branch.length <- c(x@root.edge, x@edge.length) # root.edge is not an edge length
  branch.length <- edgeLength(tree)[idx]
  if (is.null(edgeLength(tree))) {
    branch.length <- rep(NA, length(node))
  }
  names(branch.length) <- node
  return(branch.length)
}

# tip length extractor
tipLength <- function(phy) {
  tip.length <- edgeLength(phy, nodeId(phy, type="tip"))
  names(tip.length) <- tipLabels(phy)
  return(tip.length)
}

# abundance extractor
abundance <- function(phy) {
  abund <- tipData(phy)$abundance
  if (is.null(abund)) abund <- rep(NA_real_, nTips(phy))
  names(abund) <- row.names(tipData(phy))
  return(abund)
}

# abundance assignment function
`abundance<-` <- function(phy, value) {
  tipData(phy)$abundance <- value
  return(phy)
}

# minTL extractor
minTL <- function(phy) {
  minTL <- tdata(phy)$minTL
  if (!is.null(minTL)) {
    names(minTL) <- row.names(tdata(phy))
  }
  return(minTL)
}

# minTL assignment function
`minTL<-` <- function(phy, value) {
# MAKE SURE THIS MATCHES WITH THE TAXA ORDER IN THE TREE!
  if (!is.numeric(value)) {
    stop("minTL values must be a numeric vector")
  } else if (length(value)!=nrow(tdata(phy))) {
    stop("number of minTL values must equal number of species")
  }
  phy@tip.data <- replace(phy@tip.data, "minTL", value)
  return(phy)
}

# genera extractor
genera <- function(phy) {
  #From taxa names in tree, remove "_" and species name after
  gsub("_.*$", "", tipLabels(phy))
}

