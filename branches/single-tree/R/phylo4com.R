setClass("phylo4com",
    representation(subtrees="list"),
    prototype = list(subtrees = list()),
    validity = checkPhylo4,
    contains="phylo4d")

## create single tree with all data in it, plus a list of subtrees (sans
## data) for each unique species composition
setGeneric("phylo4com", function(x, n, ...) {
    standardGeneric("phylo4com")
})

setMethod("phylo4com", c("phylo", "ANY"),
  function(x, n, ..., check.node.labels="keep") {
    x <- phylo4(x, check.node.labels)
    phylo4com(x, n, ...)
})

setMethod("phylo4com", c("phylo4", "numeric"),
  function(x, n, communityID, species,
    missing=c("warn", "OK", "fail")) {

    ## create site-by-species abundance matrix
    comm <- factor(communityID, levels=unique(communityID))
    taxa <- factor(species)
    dat <- matrix(0, nrow=nlevels(taxa), ncol=nlevels(comm),
        dimnames=list(levels(taxa), unique(comm)))
    dat[cbind(taxa, comm)] <- n

    ## hand off to the phylo4com matrix method
    phylo4com(x, dat, missing)

})

setMethod("phylo4com", c("phylo4", "matrix"),
  function(x, n, missing=c("warn", "OK", "fail")) {

    taxa <- rownames(n)
    if (any(duplicated(taxa))) {
        stop("duplicated taxa are not permitted")
    }
    comm <- colnames(n)
    if (any(duplicated(comm))) {
        stop("duplicated community IDs are not permitted")
    }

    phy <- subset(x, tips.include=taxa)
    phyd <- addData(phy, tip.data=n, extra.data="OK")
    phylo4com(phyd, cols=comm)

})

setMethod("phylo4com", c("phylo4", "data.frame"),
  function(x, n, missing=c("warn", "OK", "fail")) {
    n <- as.matrix(n)
    phylo4com(x, n, missing)
})

setMethod("phylo4com", c("phylo4d", "missing"),
  function(x, n, cols) {

    if (missing(cols)) {
        cols <- names(tipData(x))
    }
    if (is.null(cols)) {
        res <- as(x, "phylo4com")
        res@metadata$comms <- NULL
        return(res)
    }
    x@metadata$comms <- setNames(rep(NA, length(cols)),
        make.names(cols))

    # create trees for each unique community wrt composition only
    P <- presence(x)
    phy <- extractTree(x)
    ids <- as.character(as.numeric(factor(sapply(P, paste,
        collapse=""))))
    subtrees <- lapply(P[!duplicated(ids)],
        function(n) subset(phy, rownames(P)[n %in% 1]))
    names(subtrees) <- ids[!duplicated(ids)]

    res <- as(x, "phylo4com")
    res@subtrees <- subtrees
    res@metadata$comms[] <- ids
    
    return(res)
})


## older approach: create a single phylo4d object with all data in it
#phylo4com <- function(communityID, species, abundance, tree,
#    missing=c("warn", "OK", "fail")) {
#
#    comm <- factor(communityID, levels=unique(communityID))
#    taxa <- factor(species)
#    dat <- matrix(0, nrow=nlevels(taxa), ncol=nlevels(comm),
#      dimnames=list(levels(taxa), unique(comm)))
#    dat[cbind(taxa, comm)] <- abundance
#
#    subtree <- subset(tree, tips.include=levels(taxa))
#
#    res <- addData(subtree, tip.data=dat, extra.data="OK")
#    res@metadata$comms <- make.names(unique(comm))
#    return(res)
#
#}

## similar to above, but cast function from reshape package
#phylo4com <- function(communityID, species, abundance, tree,
#    missing=c("warn", "OK", "fail")) {
#
#    raw <- data.frame(.taxa=factor(species), .comm=factor(communityID),
#      value=abundance, stringsAsFactors=FALSE)
#    dat <- cast(.taxa ~ .comm, data=raw)
#    row.names(dat) <- dat[[1]]
#    dat <- dat[-1]
#
#    subtree <- subset(tree, tips.include=row.names(dat))
#
#    addData(subtree, tip.data=dat,
#      metadata=list(comms=make.names(colnames(dat))), extra.data="OK")
#
#}
