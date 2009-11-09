##
## PAE methods
##

setGeneric("pae", function(x, na.rm=TRUE) {
    standardGeneric("pae")
})

setMethod("pae", signature(x="phylo4d"), function(x, na.rm=TRUE) {
    phyc <- phylo4com(x)
    pae(phyc, na.rm=na.rm)
})

setMethod("pae", signature(x="phylo4com"), function(x, na.rm=TRUE) {

    N <- abundance(x)
    comms <- x@metadata$comms
    if (is.null(comms)) {
        stop("no community data specified in phylo4com object")
    }
    subtrees <- x@subtrees[unique(as.character(comms))]

    # now for each subtree...
    #PD <- sapply(subtrees, pd)[as.character(comms)]
    PD <- pd(x)
    tmp <- setNames(rep(0, nTips(x)), tipLabels(x))
    TL <- lapply(subtrees, function(tree) {
        res <- tipLength(tree)
        tmp[match(names(res), names(tmp))] <- res
        tmp
    })
    TL <- do.call("cbind", TL[as.character(comms)])

    numer <- PD + colSums(TL * (N - 1))
    denom <- PD + (colSums(N, na.rm = na.rm) / colSums(presence(x),
        na.rm=na.rm) - 1) * colSums(TL)
    res <- numer/denom
    names(res) <- names(comms)
    return(res)

})
