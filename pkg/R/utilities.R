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
tipLength <- function(tree) {
  tip.length <- edgeLength(tree, nodeId(tree, type="tip"))
  names(tip.length) <- tipLabels(tree)
  return(tip.length)
}

# abundance extractor
abundance <- function(tree) {
  abund <- tipData(tree)$abundance
  if (is.null(abund)) abund <- rep(NA_real_, nTips(tree))
  names(abund) <- row.names(tipData(tree))
  return(abund)
}

# abundance assignment function
`abundance<-` <- function(tree, value) {
  tipData(tree)$abundance <- value
  return(tree)
}

# minTL extractor
minTL <- function(tree) {
  minTL <- tdata(tree)$minTL
  if (!is.null(minTL)) {
    names(minTL) <- row.names(tdata(tree))
  }
  return(minTL)
}

# minTL assignment function
`minTL<-` <- function(tree, value) {
# MAKE SURE THIS MATCHES WITH THE TAXA ORDER IN THE TREE!
  if (!is.numeric(value)) {
    stop("minTL values must be a numeric vector")
  } else if (length(value)!=nrow(tdata(tree))) {
    stop("number of minTL values must equal number of species")
  }
  tree@tip.data <- replace(tree@tip.data, "minTL", value)
  return(tree)
}

# genera extractor
genera <- function(tree) {
  #From taxa names in tree, remove "_" and species name after
  gsub("_.*$", "", tipLabels(tree))
}
