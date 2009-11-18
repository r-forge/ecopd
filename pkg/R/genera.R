setGeneric("genera", function(x, ...) {
    standardGeneric("genera")
})

setMethod("genera", c("phylo4"),
  function(x) {
    # From taxa names in tree, remove "_" and species name after
    gsub("_.*$", "", tipLabels(x))
})

setMethod("genera", c("phylo4com"),
  function(x, community) {
    if (missing(community)) {
        community <- communities(x)
    }
    doNotExist <- !community %in% communities(x)
    if (any(doNotExist)) {
        stop("one or more communities not found in x: ",
            paste(community[doNotExist], collapse=", "))
    }
    g <- lapply(getSubtrees(x, community), genera)
    names(g) <- community
    return(g)
})
