# ecophylo object constructor
ecophylo <- function (communityID, species, abundance, tree) {

  if (all(species %in% tree$tip.label)) {
    missing <- NULL
  } else {
    missing <- setdiff(species, tree$tip.label)
  }

  communities <- split(data.frame(species, abundance), 
    factor(communityID, unique(as.character(communityID))))
  
  spAbundances <- sapply(communities, function(community) {
    vec <- community$abundance
    names(vec) <- community$species
    return(vec)
    }, SIMPLIFY=FALSE)

  obj <- list(data=spAbundances, tree=tree, missing=missing)
  class(obj) <- "ecophylo"

  return(obj)

}

# ecophylo print method
print.ecophylo <- function(x, ...) {
  cat("\n--- Phylogeny ---\n")
  print(x$tree)
  cat("\n--- Ecological data ---\n")
  ncom <- length(x$data)
  nspp <- length(unique(unlist(sapply(x$data, names))))
  cat(paste(nspp, "total species across", ncom, "communities\n"))
  cat("\n")
  if (!is.null(x$missing)) {
    cat("Species in data but not in tree:\n")
    writeLines(strwrap(paste(x$missing, collapse = ", "), indent=4,
    exdent=8))
  }
  cat("\n")
  return(x)
}
