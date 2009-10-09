phylo4com <- function(communityID, species, abundance, tree,
  missing=c("warn", "OK", "fail")) {

  ## species should only ever match against tip *labels* (not node IDs)
  species <- as.character(species)

  ## identify and exclude any species missing from tree
#  tips <- getNode(tree, unique(species), type="tip", missing=missing)
#  isInTree <- species %in% names(tips[!is.na(tips)])
#  species <- species[isInTree]
#  abundance <- abundance[isInTree]
#  communityID <- communityID[isInTree]

  communities <- split(data.frame(species, abundance,
    stringsAsFactors=FALSE), factor(communityID,
    unique(as.character(communityID))))
  
  subtrees <- lapply(communities, function(community) {
    subtree <- subset(tree, tips.include=community$species)
    cdata <- with(community, data.frame(abundance, row.names=species))
    tipData(subtree, extra.data="OK") <- cdata
    subtree
  })

  return(subtrees)

  # alternative conceptualization...
  #result <- list(tree=tree, community=subtrees)
  #class(result) <- "phylo4com"
  #return(result)

}
