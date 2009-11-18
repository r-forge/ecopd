##
## ED and related methods
##

setGeneric("ed", function(x, na.rm=TRUE) {
    standardGeneric("ed")
})

setMethod("ed", signature(x="phylo4d"), function(x, na.rm=TRUE) {
    phyc <- phylo4com(x)
    ed(phyc, na.rm=na.rm)
})

# TODO: This function includes its own code for not counting root edge
# length. Maybe this should maybe be done at a higher level?
setMethod("ed", signature(x="phylo4com"), function(x, na.rm=TRUE) {

    comms <- x@metadata$comms
    if (is.null(comms)) {
        stop("no community data specified in phylo4com object")
    }
    subtrees <- x@subtrees[unique(as.character(comms))]
    .edi <- function(tree) {
        # set length of root edge to zero
        edgeLength(tree)[edgeId(tree, "root")] <- 0

        all.nodes <- nodeId(tree, type = "all")
        des <- descendants(tree, all.nodes, type="tips")
        nv <- edgeLength(tree, all.nodes) / sapply(des, length)
        names(nv) <- all.nodes

        tip.nodes <- nodeId(tree, "tip")
        anc <- ancestors(tree, tip.nodes, "ALL")

        res <- sapply(anc, function(n) sum(nv[as.character(n)], na.rm=TRUE))
        names(res) <- tipLabels(tree)
        res
    }
    res <- lapply(subtrees, .edi)[as.character(comms)]
    names(res) <- names(comms)
    return(res)

})

setGeneric("hed", function(x, na.rm=TRUE) {
    standardGeneric("hed")
})

setMethod("hed", signature(x="phylo4d"), function(x, na.rm=TRUE) {
    phyc <- phylo4com(x)
    hed(phyc, na.rm=na.rm)
})

setMethod("hed", signature(x="phylo4com"), function(x, na.rm=TRUE) {
    ED <- ed(x)
    PD <- pd(x)
    res <- sapply(seq_len(length(PD)), function(i) {
        scaledED <- ED[[i]] / PD[[i]]
        -sum(scaledED * log(scaledED))
    })
    names(res) <- names(PD)
    return(res)
})

setGeneric("eed", function(x, na.rm=TRUE) {
    standardGeneric("eed")
})

setMethod("eed", signature(x="phylo4d"), function(x, na.rm=TRUE) {
    phyc <- phylo4com(x)
    eed(phyc, na.rm=na.rm)
})

setMethod("eed", signature(x="phylo4com"), function(x, na.rm=TRUE) {
    comms <- x@metadata$comms
    if (is.null(comms)) {
        stop("no community data specified in phylo4com object")
    }
    subtrees <- x@subtrees[unique(as.character(comms))]
    hed(x) / log(sapply(subtrees, nTips)[as.character(comms)])
})

setGeneric("aed", function(x, na.rm=TRUE) {
    standardGeneric("aed")
})

setMethod("aed", signature(x="phylo4d"), function(x, na.rm=TRUE) {
    phyc <- phylo4com(x)
    aed(phyc, na.rm=na.rm)
})

# TODO: This function includes its own code for not counting root edge
# length. Maybe this should maybe be done at a higher level?
setMethod("aed", signature(x="phylo4com"),
  function(x, na.rm=TRUE) {

    comms <- x@metadata$comms
    if (is.null(comms)) {
        stop("no community data specified in phylo4com object")
    }
    subtrees <- x@subtrees[unique(as.character(comms))]
    .isD <- function(tree) {
        # get all node IDs, but excluding root node
        nonroot.nodes <- setdiff(nodeId(tree), rootNode(tree))
        # Create logical matrix indicating which tips (in columns) are
        # descendants of each node (in rows), self-inclusive
        t(sapply(ancestors(tree, tipLabels(tree), "ALL"),
            function(n) nonroot.nodes %in% n))
    }

    .elen <- function(tree) {
        nonroot.nodes <- setdiff(nodeId(tree), rootNode(tree))
        edgeLength(tree, nonroot.nodes)
    }
    isDescendant <- lapply(subtrees, .isD)[as.character(comms)]
    edge.length <- lapply(subtrees, .elen)[as.character(comms)]
    N <- abundance(x)

    res <- lapply(seq_along(N), function(i) {
        spp <- row.names(isDescendant[[i]])        
        dAbund <- N[spp, i] * isDescendant[[i]]

        # Calculate individual-based AED of each species
        AED <- colSums(edge.length[[i]] * t(prop.table(dAbund, margin=2)))
        AED/N[spp, i]
    })
    names(res) <- names(comms)
    return(res)

})

#    .isD <- function(tree, template) {
#        # get all node IDs, but excluding root node
#        nonroot.nodes <- setdiff(nodeId(tree), rootNode(tree))
#        # Create logical matrix indicating which tips (in columns) are
#        # descendants of each node (in rows), self-inclusive
#        res <- t(sapply(ancestors(tree, tipLabels(tree), "ALL"),
#          function(n) nonroot.nodes %in% n))
#        template[match(rownames(res), rownames(template)), ] <- res
#        template
#    }
#
#    .el <- function(tree) {
#        nonroot.nodes <- setdiff(nodeId(tree), rootNode(tree))
#        edgeLength(tree, nonroot.nodes)
#    }
#    tmp <- setNames(rep(NA, nTips(x)), tipLabels(x))
#    tmp <- matrix(NA, nrow=nTips(x), ncol=nNodes(x)+nTips(x)-1)
#    rownames(tmp) <- tipLabels(x)
#    isDescendant <- lapply(subtrees, .isD, tmp)
#    isDescendant <- array(unlist(isDescendant), dim=c(nTips(x),
#        nTips(x)+nNodes(x)-1, length(subtrees)))
#
#    # Create vector of ancestral edge lengths
#    edge.length <- sapply(subtrees, .elen)
#
#    # Create matrix containing number of individuals of each species
#    # descending from each interior node
#    N <- as.matrix(abundance(x))
#    dAbund <- sweep(isDescendant, c(1,3), N, "*")
#
#    # Calculate individual-based AED of each species
#    pt <- prop.table(dAbund, margin=c(2,3))
#    AED <- sweep(pt, c(2,3), edge.length, "*")
#    AED <- apply(AED, 3, rowSums) / N
#
#    return(AED)

setGeneric("haed", function(x, na.rm=TRUE) {
    standardGeneric("haed")
})

setMethod("haed", signature(x="phylo4d"), function(x, na.rm=TRUE) {
    phyc <- phylo4com(x)
    haed(phyc, na.rm=na.rm)
})

setMethod("haed", signature(x="phylo4com"), function(x, na.rm=TRUE) {
    # Recast AED in terms of individuals
    AED <- aed(x)
    PD <- pd(x)
    N <- abundance(x)
    scaled.AED <- lapply(seq_along(N), function(i) {
        spp <- names(AED[[i]])        
        rep(unname(AED[[i]]), N[spp, i]) / PD[[i]]
    })
    res <- sapply(scaled.AED, function(x) -sum(x * log(x)))
    names(res) <- names(AED)
    return(res)
})

setGeneric("eaed", function(x, na.rm=TRUE) {
    standardGeneric("eaed")
})

setMethod("eaed", signature(x="phylo4d"), function(x, na.rm=TRUE) {
    phyc <- phylo4com(x)
    eaed(phyc, na.rm=na.rm)
})

setMethod("eaed", signature(x="phylo4com"), function(x, na.rm=TRUE) {
  haed(x) / log(colSums(abundance(x)))
})

setGeneric("value",
  function(x, na.rm=TRUE) {
    standardGeneric("value")
})

setMethod("value", signature(x="phylo4d"),
  function(x, na.rm=TRUE) {
    phyc <- phylo4com(x)
    value(phyc, na.rm=na.rm)
})

setMethod("value", signature(x="phylo4com"),
  function(x, na.rm=TRUE) {
    AED <- aed(x)
    lapply(AED, function(x) x/sum(x))
})

