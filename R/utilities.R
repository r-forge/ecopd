# abundance extractor
abundance <- function(tree) {
  abund <- tree$abundance
  if (!is.null(abund)) {
    names(abund) <- tree$tip.label
  }
  return(abund)
}

# abundance assignment function
`abundance<-` <- function(tree, value) {
  if (!is.numeric(value)) {
    stop("abundance values must be a numeric vector")
  } else if (length(value)!=length(tree$tip.label)) {
    stop("number of abundance values must equal number of species")
  }
  tree$abundance <- value
  return(tree)
}

# minTL extractor
minTL <- function(tree) {
  minTL <- tree$minTL
  if (!is.null(minTL)) {
    names(minTL) <- tree$tip.label
  }
  return(minTL)
}

# minTL assignment function
`minTL<-` <- function(tree, value) {
# MAKE SURE THIS MATCHES WITH THE TAXA ORDER IN THE TREE!
  if (!is.numeric(value)) {
    stop("minTL values must be a numeric vector")
  } else if (length(value)!=length(tree$tip.label)) {
    stop("number of minTL values must equal number of species")
  }
  tree$minTL <- value
  return(tree)
}

# tip length extractor
tipLengths <- function(tree) {
  tipRowID <- which(tree$edge[, 2] <= Ntip(tree))
  tipLen <- tree$edge.length[tipRowID]
  names(tipLen) <- tree$tip.label[tree$edge[tipRowID,2]]
  return(tipLen)
}

# genera extractor
genera <- function(tree) {
  #From taxa names in tree, remove "_" and species name after
  gsub("_.*$", "", tree$tip.label)
}

# phylo subset function
subset.phylo <- function(x, tipNames, checkMissing=TRUE, ...) {

  # Get all tip names from tree
  allTipNames <- x$tip.label

  # If desired, ensure supplied tip names exist in the tree
  if (checkMissing) {
    missingNames <- setdiff(tipNames, allTipNames)
    if (length(missingNames)>0) {
      stop("one or more tip names not found in main tree:\n",
        paste(missingNames, collapse=", "))
    }
  }

  foundInTree <- allTipNames %in% unique(tipNames)

  # Make sure the subset will have at least two tips
  if (sum(foundInTree) < 2) {
  	stop("subset tree has fewer than 2 species")
  }

  # Drop tips that are not among the supplied tip nmaes
  dropme <- allTipNames[!foundInTree]
  if (length(dropme)>0) {
    sub.tree <- drop.tip(x, dropme)
  } else {
    sub.tree <- x
  }

  return(sub.tree)

}

