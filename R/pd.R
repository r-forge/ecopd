# basic pd calculation
pd <- function(tree) {
  sum(tree$edge.length)
}

# lookup function for minimum tip length
getMinTL <- function(tree, genera) {

  # Families, Supertree, and LookupTL are all system data built into
  # the package

  # Supplied tree must be ultrametric
  if (!is.ultrametric(tree)) {
    stop(deparse(substitute(tree)), " is not ultrametric")
  }

  # TODO: if genus not in Families vector, need to give warning.
  families <- Families[genera]

  # if Family is not in LookupTL give a warning, and later
  # use minTL based on average across all families in the lookup table
  if (!all(families %in% row.names(LookupTL))) {
    warning("one or more taxa missing from lookup table; using mean minTL")
  }

  familiesInSupertree <- LookupTL[families, "supertree"]
  familiesInSupertree <- as.character(familiesInSupertree)

  # Subset the supertree using families from the user-supplied tree. If
  # any user-supplied taxa cannot be matched to families in the
  # supertree, they are simply ignored
  subsupertree <- subset(Supertree, na.omit(familiesInSupertree))
  subsupertree.maxLength <- max(cophenetic(subsupertree))/2
  tree.maxLength <- max(cophenetic(tree))/2

  tableTL <- LookupTL[familiesInSupertree, "minTL"]

  # if any TLs are 0 or NA, give warning and use average minTL across
  # all families in the LookupTL table
  if(any(tableTL==0 | is.na(tableTL))) {
    numNA <- sum(tableTL==0 | is.na(tableTL))
    warning("Using meanTL for ", numNA, " tip", if(numNA>1) "s") 
    # calculate average minTL across the entire LookupTL table,
    # excluding any non-positive or NA values
    meanMinTL <- mean(LookupTL$minTL[LookupTL$minBL > 0], na.rm=TRUE)
    tableTL[is.na(tableTL)] <- meanMinTL  
    tableTL[tableTL<=0] <- meanMinTL  
  }

  lookupTL <- tableTL * (tree.maxLength / subsupertree.maxLength)
  actualTL <- tipLengths(tree)
  minTL <- ifelse(lookupTL < actualTL, lookupTL, actualTL)

  return(minTL)
 
}

# function to weight tip length based on abundance
weightByAbund <- function(tree, method=c("polytomy", "yule")) {

  method <- match.arg(method)

  if (is.null(abundance(tree))) {
    stop("tree contains no abundance information")
  }

  if (is.null(minTL(tree))) {
    stop("tree contains no minTL information")
  }

  tipRowID <- which(tree$edge[, 2] <= Ntip(tree))
  tipLen <- tree$edge.length[tipRowID]
  tipID <- tree$edge[tipRowID,2]

  n <- abundance(tree)[tipID]
  minLength <- minTL(tree)[tipID]
 
  # Test statement:
  if (!identical(names(n), names(minLength))) {
    stop("mismatch between abundance and minTL vectors")
  }

  if (method=="polytomy") {
    tree$edge.length[tipRowID] <- tipLen + (n-1) * minLength
  } else if (method=="yule") {
    C <- 0.57722
    tree$edge.length[tipRowID] <- tipLen - minLength +
      (minLength * (n-1) / (log(n) + C-1))
  }

  return(tree)

}

# calculate PD value for a single community
commPD <- function(tree, method=c("traditional", "polytomy", "yule")) {
  method <- match.arg(method)
  if (method == "traditional") {
    pd(tree)
  } else {
    pd(weightByAbund(tree, method))
  } 
}

ecoPD <- function (ecophylo, method=c("traditional", "polytomy",
  "yule"), use.supertree=TRUE, nsims=0) {

  method <- match.arg(method)

  if (use.supertree) {
    #message("Looking up minTL values from supertree")
    minTL(ecophylo$tree) <- getMinTL(ecophylo$tree, genera(ecophylo$tree))
  } else {
    if (is.null(minTL(ecophylo$tree))) {
      stop("ecophylo tree must have minTL values if use.supertree=FALSE")
    }
  }

  # define pd function to apply to each community
  pdfunc <- function(sub.tree, method, nsims) {
      pd <- commPD(sub.tree, method=method)
      if (nsims>0) {
        pd.sims <- replicate(nsims, {
          abundance(sub.tree) <- sample(abundance(sub.tree))
          commPD(sub.tree, method=method)
        })
        pd.quants <- quantile(pd.sims, c(0.5, 0.025, 0.975))
        names(pd.quants) <- c("median", "0.025", "0.975")
        pd <- c(pd.obs=pd, mean = mean(pd.sims), pd.quants)
      }
      pd
    }

  answer <- doByCommunity(ecophylo, pdfunc, simplify=TRUE,
    method=method, nsims=nsims)
  attr(answer, "method") <- method
  return(answer)

}

# calculate PD values for a whole set of communities
ecoPD.old <- function (ecophylo, method, use.supertree=TRUE, nsims=0) {

  communities <- ecophylo$data
  tree <- ecophylo$tree
  if (use.supertree) {
    #message("Looking up minTL values from supertree")
    minTL(tree) <- getMinTL(tree, genera(tree))
  } else {
    if (is.null(minTL(tree))) {
      stop("ecophylo tree must have minTL values if use.supertree=FALSE")
    }
  }

  PD <- lapply(communities, function(community) {
    sub.tree <- subset(tree, names(community), checkMissing=FALSE)
    minTL(sub.tree) <- minTL(tree)[sub.tree$tip.label]
    abundance(sub.tree) <- community[sub.tree$tip.label]
    pd <- c(pd=commPD(sub.tree, method=method))
    if (nsims>0) {
      pd.sims <- replicate(nsims, {
        abundance(sub.tree) <- sample(abundance(sub.tree))
        commPD(sub.tree, method=method)
      })
      pd.quants <- quantile(pd.sims, c(0.5, 0.025, 0.975))
      names(pd.quants) <- c("median", "0.025", "0.975")
      pd <- c(pd, pd.mean = mean(pd.sims), pd.quants)
    }
    pd
  })

  return(data.frame(do.call("rbind", PD)))

}
