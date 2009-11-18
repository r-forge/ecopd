##
## IAC methods
##

setGeneric("iac", function(x, na.rm=TRUE) {
    standardGeneric("iac")
})

setMethod("iac", signature(x="phylo4d"), function(x, na.rm=TRUE) {
    phyc <- phylo4com(x)
    iac(phyc, na.rm=na.rm)
})

setMethod("iac", signature(x="phylo4com"), function(x, na.rm=TRUE) {

    comms <- x@metadata$comms
    if (is.null(comms)) {
        stop("no community data specified in phylo4com object")
    }
    subtrees <- x@subtrees[unique(as.character(comms))]
    .denom <-  function(tree, template) {
        # Count number of lineages originating at each internal node
        # (i.e. number of splits)
        int.nodes <- nodeId(tree, "internal")
        nSplits <- sapply(int.nodes, function(x) length(children(tree,
            x)))
        names(nSplits) <- int.nodes
        # For each tip, take the product of the number of splits across
        # all of its ancestral nodes
        res <- sapply(ancestors(tree, tipLabels(tree)), function(x)
            prod(nSplits[as.character(x)]))
        template[match(names(res), names(template))] <- res
        template
    }

    # now for each subtree...
    tmp <- setNames(rep(NA, nTips(x)), tipLabels(x))
    denom <- lapply(subtrees, .denom, tmp)
    denom <- do.call("cbind", denom[as.character(comms)])
    nnodes <- sapply(x@subtrees, nNodes)[as.character(comms)]

    # Calculate expected number of individuals under null hypothesis
    # of equal allocation to each lineage at each (node) split 
    N <- abundance(x)
    expected <- t(colSums(N, na.rm=na.rm) / t(denom))

    # IAC: summed absolute difference between expected and observed
    # abundances, divided by number of nodes
    colSums(abs(expected - N), na.rm=na.rm) / nnodes

})

