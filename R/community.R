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

# simpson's index with and without phylogenetic distances using
# commmunity matrix 'x', species names as colnames and community names
# as rownames. phy=FALSE returns traditional simpson's index
simp.phy <- function(x, tr, phy=TRUE) {

  x <- as.matrix(x)
  x <- x[, order(colnames(x))]
  total <- apply(x, 1, sum)
  x <- sweep(x, 1, total, "/")
    
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

