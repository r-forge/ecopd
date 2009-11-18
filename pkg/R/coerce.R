setMethod("phylo4d", c("phylo4com"),
  function(x, community) {
    if (missing(community)) {
        return(as(x, "phylo4d"))
    }
    if (length(community)!=1) {
        stop("a single community label must be provided")
    }

    ## get the tree
    tree <- phylo4(x, community)

    ## get abundance values, and drop species with zero abundance
    N <- abundance(x, community)
    N <- N[N[[community]]!=0, , drop=FALSE]

    ## combine and return phylo4d object
    phylo4d(tree, tip.data=N)
})

setMethod("phylo4", c("phylo4com"),
  function(x, community) {
    if (missing(community)) {
        return(as(x, "phylo4"))
    }
    if (length(community)!=1) {
        stop("a single community label must be provided")
    }
    getSubtrees(x, community)[[1]]
})

## internal helper function for extracting a list of the subtrees
## associated with each communities
getSubtrees <- function(x, community) {
    communities <- x@metadata$comms
    if (missing(community)) {
        community <- names(communities)
    }
    doNotExist <- !community %in% names(communities)
    if (any(doNotExist)) {
        stop("one or more communities not found in x: ",
            paste(community[doNotExist], collapse=", "))
    }
    res <- x@subtrees[communities[names(communities) %in% community]]
    names(res) <- community
    return(res)
}
    
