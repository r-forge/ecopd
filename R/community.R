# Generate site-by-species matrix
siteBySpecies <- function(phylo4com) {
  spp <- unique(unlist(lapply(phylo4com, phylobase::labels)))
  sapply(phylo4com, function(x) {
    vec <- tdata(x)[spp, "abundance"]
    vec[is.na(vec)] <- 0
    names(vec) <- spp
    vec
  })
}

# Calculate species richness
richness <- function(phylo4com) {
  sapply(phylo4com, function(x) nrow(tdata(x)))
}

