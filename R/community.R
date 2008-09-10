# Apply arbitrary function to each community subtree
doByCommunity <- function(ecophylo, FUN, simplify=TRUE, ...) {
  communities <- ecophylo$data
  tree <- ecophylo$tree
  answer <- lapply(communities, function(community) {
    sub.tree <- subset(tree, names(community), checkMissing=FALSE)
    if (!is.null(minTL(tree))) {
      minTL(sub.tree) <- minTL(tree)[sub.tree$tip.label]
    }
    abundance(sub.tree) <- community[sub.tree$tip.label]
    FUN(sub.tree, ...)
  })
  if (simplify && length(answer) && length(common.len <-
    unique(unlist(lapply(answer, length)))) == 1) {
    if (common.len == 1) {
      answer <- unlist(answer, recursive = FALSE)
    } else if (common.len > 1) {
      answer <- do.call("rbind", answer)
    }
  }
  return(answer)
}

# Generate site-by-species matrix
siteBySpecies <- function(ecophylo) {
  sapply(ecophylo$data, function(x) {
    spp <- unique(unlist(lapply(ecophylo$data, names)))
    vec <- x[spp]
    vec[is.na(vec)] <- 0
    names(vec) <- spp
    vec
  })
}

# Calculate species richness
richness <- function(ecophylo, onlyInTree=FALSE) {
  if (onlyInTree) {
    sapply(ecophylo$data, function(species) sum(names(species) %in%
      ecophylo$tree$tip.label))
  } else {
    sapply(ecophylo$data, length)
  }
}
