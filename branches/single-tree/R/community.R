# Generate site-by-species matrix
siteBySpecies <- function(phylo4com, presence.only=FALSE,
  transpose=FALSE) {

  ## create and populate the matrix
  spp <- unique(unlist(lapply(phylo4com, tipLabels)))
  mat <- sapply(phylo4com, function(x) {
    vec <- tipData(x)[spp, "abundance"]
    vec[!spp %in% tipLabels(x)] <- 0
    vec
  })
  rownames(mat) <- spp

  if (presence.only) {
    mat[mat>0] <- 1
    mat[mat<=0] <- 0
  }

  if (!transpose) mat <- t(mat)

  return(mat)

}

# Calculate species richness
richness <- function(phylo4com, na.rm=FALSE) {
  sapply(phylo4com, function(x) sum(abundance(x)>0, na.rm=na.rm))
}
