# internal function to compute distances from one node to all others
# NOTE: assumes node IDs are always 1:n (currently true of phylo4 objs)
dijkstra <- function(phy, node) {
  edge <- edges(phy, drop.root=TRUE)
  elen <- edgeLength(phy)
  nodes <- nodeId(phy, "all")
  n <- length(nodes)
  d <- rep(Inf, length=n)
  names(d) <- nodes
  d[node] <- 0
  while (length(nodes)>0) {
    u <- nodes[which.min(d[nodes])[[1]]]
    if (is.infinite(d[u])) break
    nodes <- nodes[-match(u, nodes)]
    anc <- edge[edge[,2]==u,1]
    for (parent in anc[anc %in% nodes]) {
      alt <- d[u] + elen[paste(parent, u, sep="-")]
      if (alt < d[parent]) {
        d[parent] <- alt
      }
    }
    des <- edge[edge[,1]==u,2]
    for (child in des[des %in% nodes]) {
      alt <- d[u] + elen[paste(u, child, sep="-")]
      if (alt < d[child]) {
        d[child] <- alt
      }
    }
  }
  d
}

# tip length extractor
tipLength <- function(phy, from=c("parent", "root")) {
  from <- match.arg(from)
  tips <- nodeId(phy, type="tip")
  if (from=="parent") {
    len <- edgeLength(phy, tips)
  } else if (from=="root") {
    len <- dijkstra(phy, rootNode(phy))
    len <- len[match(tips, names(len))]
  }
  names(len) <- tipLabels(phy)
  return(len)
}

# abundance extractor
abundance <- function(phy, comm) {
  communities <- names(phy@metadata$comms)
  if (missing(comm)) {
    return(tipData(phy)[communities])
  }
  doNotExist <- !comm %in% communities
  if (any(doNotExist)) {
    stop("one or more communities not found in phy: ",
      paste(comm[doNotExist], collapse=", "))
  }
  return(tipData(phy)[comm])
}

# abundance assignment function
`abundance<-` <- function(phy, comm, tip, value) {
  if (!is.atomic(comm) || length(comm)!=1) {
    stop("comm must be a vector of length 1")
  } else if (!comm %in% names(phy@metadata$comms)) {
    stop(paste("community", comm, "not found in phy", sep=" "))
  }
  if (missing(tip)) tip <- tipLabels(phy)
  tipData(phy)[tip, comm] <- value
  return(phy)
}

presence <- function(phy) {
    N <- abundance(phy)
    N[N > 0] <- 1
    N[N <= 0] <- 0
    N[is.na(N)] <- 0
    N
}

# minTL extractor
minTL <- function(phy) {
  minTL <- tipData(phy)$minTL
  if (!is.null(minTL)) {
    names(minTL) <- row.names(tipData(phy))
  }
  return(minTL)
}

# Note: function assumes values are in nodeId(phy, "tip") order
# minTL assignment function
`minTL<-` <- function(phy, value) {
  if (!is.numeric(value)) {
    stop("minTL values must be a numeric vector")
  } else if (length(value)!=nTips(phy)) {
    stop("number of minTL values must equal number of species")
  }
  tipData(phy)[["minTL"]] <- value
  return(phy)
}

# genera extractor
genera <- function(phy) {
  #From taxa names in tree, remove "_" and species name after
  gsub("_.*$", "", tipLabels(phy))
}

## this works as implementation of dist.nodes for phylo4 objects, albeit
## about 1.5x slower than dist.nodes
pairdist <- function(phy, type=c("all", "tip"), use.labels=FALSE) {

  if (hasPoly(phy) || hasRetic(phy)) {
    stop("pairdist can't currently handle polytomies or reticulation")
  }

  type <- match.arg(type)
  edge <- edges(phy, drop.root=TRUE)
  elen <- edgeLength(phy)
  elen <- elen[!names(elen) %in% edgeId(phy, "root")]
  nodes <- nodeId(phy, "all")
  n <- length(nodes)
  d <- matrix(NA_real_, nrow=n, ncol=n)
  diag(d) <- 0
  d[edge] <- d[edge[,2:1]] <- elen
  ntips <- nTips(phy)
  tips <- nodeId(phy, "tip")
  tip.parents <- sapply(tips, function(n) edge[edge[,2]==n, 1])
  nodes.todo <- tip.parents[duplicated(tip.parents)]
  done <- tips
  root <- if (isRooted(phy)) rootNode(phy) else ntips+1
  descendants <- matrix(edge[order(edge[,1]),2], nrow=2)
  ancestors <- edge[match(nodes, edge[,2]), 1]

  while(length(nodes.todo)>0) {
    nod <- nodes.todo[1]
    if (nod==root && length(nodes.todo)>1) {
      nod <- nodes.todo[2]
    }
    des1 <- descendants[1, nod-ntips]
    des2 <- descendants[2, nod-ntips]
    if (des1>ntips) des1 <- which(!is.na(d[des1, ]))
    if (des2>ntips) des2 <- which(!is.na(d[des2, ]))
    for (y in des1) {
      d[y, des2] <- d[des2, y] <- d[nod, y] + d[nod, des2]
    }
    if (nod!=root) {
      anc <- ancestors[nod]
      l <- elen[paste(anc, nod, sep="-")]
      d[des2, anc] <- d[anc, des2] <- d[nod, des2] + l
      d[des1, anc] <- d[anc, des1] <- d[nod, des1] + l
      done <- c(done, nod)
      if (all(descendants[, anc-ntips] %in% done)) {
        nodes.todo <- c(nodes.todo, anc)
      }
    }
    nodes.todo <- nodes.todo[nod!=nodes.todo]
  }

  if (type=="tip") {
    d <- d[1:ntips, 1:ntips]
  }

  if (use.labels) {
    rownames(d) <- colnames(d) <- unname(labels(phy, type=type))
  } else {
    rownames(d) <- colnames(d) <- nodeId(phy, type=type)
  }

  return(d)
}
