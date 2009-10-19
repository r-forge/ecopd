#
# --- Create test trees (ape::phylo objects) ---
#

weed.tree <- read.tree(text="(((Taraxacum_officinale:0.2601451047,(Silybum_marianum:0.1088645027,Centaurea_alba:0.1088645027)0.997000:0.151280602)0.547000:0.5232373813,Torilis_arvensis:0.783382486)0.878000:0.216617514,Trifolium_repens:1)1.000000:0.04059;")

weed.tree.notultra <- read.tree(text="(((Taraxacum_officinale:0.058133000000000004,(Silybum_marianum:0.022662,Centaurea_alba:0.02092)0.997000:0.026418999999999998)0.547000:0.075571,Torilis_arvensis:0.105052)0.878000:0.037612,Trifolium_repens:0.25731000000000004)1.000000:0.04059;")

#
# --- Test basic tree functions ---
#
 
test.phylo <- function() {
  checkTrue(class(weed.tree)=="phylo")
  checkTrue(class(weed.tree.notultra)=="phylo")
  checkTrue(is.ultrametric(weed.tree))
  checkTrue(!is.ultrametric(weed.tree.notultra))
}

test.pd <- function() {
  DEACTIVATED("no longer calculating pd for ape phylo objects")
  checkEqualsNumeric(pd(weed.tree), 3.1523920934)
  checkEqualsNumeric(commPD(weed.tree), 3.1523920934)
  checkEqualsNumeric(commPD(weed.tree, "trad"), 3.1523920934)
  checkException(commPD(weed.tree, "poly"))
  checkException(commPD(weed.tree, "yule"))

  checkEqualsNumeric(pd(weed.tree.notultra), 0.603679)
  checkEqualsNumeric(commPD(weed.tree.notultra), 0.603679)
  checkEqualsNumeric(commPD(weed.tree.notultra, "trad"), 0.603679)
  checkException(commPD(weed.tree.notultra, "poly"))
  checkException(commPD(weed.tree.notultra, "yule"))
}

test.genera <- function() {
  DEACTIVATED("no longer extracting genera from ape phylo objects")
  ans <- c("Taraxacum", "Silybum", "Centaurea", "Torilis", "Trifolium")
  checkEquals(ans, genera(weed.tree))
}

test.tipLengths <- function() {
  DEACTIVATED("no longer calculating tiplengths for ape phylo objects")
  ans <- structure(c(0.2601451047, 0.1088645027, 0.1088645027,
    0.783382486, 1), .Names = c("Taraxacum_officinale",
    "Silybum_marianum", "Centaurea_alba", "Torilis_arvensis",
    "Trifolium_repens"))
  checkEquals(ans, tipLengths(weed.tree))
}

test.subset.phylo <- function() {
  DEACTIVATED("no longer using subset for ape phylo objects")
  sub2 <- read.tree(text=
    "(Taraxacum_officinale:0.2601451047,Silybum_marianum:0.2601451047)0.547000;")
  sub3 <- read.tree(text=
    "((Taraxacum_officinale:0.2601451047,Silybum_marianum:0.2601451047)0.547000:0.7398548953,Trifolium_repens:1)1.000000;")
  checkEquals(sub2, subset(weed.tree, c("Taraxacum_officinale",
    "Silybum_marianum")))
  checkEquals(sub3, subset(weed.tree, c("Taraxacum_officinale",
    "Silybum_marianum", "Trifolium_repens")))
  checkEquals(weed.tree, subset(weed.tree, weed.tree$tip.label))
  checkException(subset(weed.tree, "Taraxacum_officinale"))
  checkException(subset(weed.tree, c("Taraxacum_officinale", "xxx")))
  checkException(subset(weed.tree, "xxx"))
}

