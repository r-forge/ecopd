phylo4com <- function(communityID, species, abundance, tree) {

  # need to fix this. my subset.phylo function automatically dropped any
  # species not found in tree, but phylobase subset returns an error.
  # thus, need to first remove any such species (with warning?) before
  # calling the subset method.
  if (all(species %in% tree@tip.label)) {
    missing <- NULL
  } else {
    missing <- setdiff(species, tree@tip.label)
    error("one or more species not found in tree")
  }

  communities <- split(data.frame(species, abundance), 
    factor(communityID, unique(as.character(communityID))))
  
  subtrees <- lapply(communities, function(community) {
    cdata <- with(community, data.frame(abundance, row.names=species))
    subtree <- phylobase::subset(tree, tips.include=community$species)
    tdata(subtree) <- cdata
    return(subtree)
  })

  return(subtrees)

  # alternative conceptualization...
  #result <- list(tree=tree, community=subtrees)
  #class(result) <- "phylo4com"
  #return(result)

}
