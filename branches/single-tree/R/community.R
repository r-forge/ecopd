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

# simpson's index with and without phylogenetic distances using
# commmunity matrix 'x', species names as colnames and community names
# as rownames. phy=FALSE returns traditional simpson's index
simp.phy <- function(x, tr, phy=TRUE) {

  x <- as.matrix(x)
  x <- x[, order(colnames(x)), drop=FALSE]
  x <- prop.table(x, 1)
    
  if (phy==TRUE){
    phy.mat <- cophenetic(tr)
    phy.mat <- phy.mat[order(rownames(phy.mat)), order(colnames(phy.mat))]
    out <- apply(x, 1, function(x) sum((x %o% x)*phy.mat))
  } else {
    bin.mat <- matrix(1, dim(x)[2], dim(x)[2])
    diag(bin.mat) <- 0
    out <- apply(x, 1, function(x) sum((x %o% x)*bin.mat))
  }   
        
  return(out) 
}

