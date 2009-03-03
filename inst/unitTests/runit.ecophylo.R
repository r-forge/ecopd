#
# --- Create test data ---
#

weed.com <- read.table(textConnection(
  "plot   cover   taxa
  A   5   Centaurea_alba
  A   5   Silybum_marianum
  A   10  Taraxacum_officinale
  A   25  Torilis_arvensis
  A   55  Trifolium_repens
  B   35  Centaurea_alba
  B   35  Silybum_marianum
  B   15  Taraxacum_officinale
  B   10  Torilis_arvensis
  B   5   Trifolium_repens
  C   20  Centaurea_alba
  C   20  Silybum_marianum
  C   20  Taraxacum_officinale
  C   20  Torilis_arvensis
  C   20  Trifolium_repens"
  ), header=TRUE)
closeAllConnections()

weed.com2 <- weed.com[c(1:5, 7:10, 13:15),]

weed.tree.notultra <- read.tree(text="(((Taraxacum_officinale:0.058133000000000004,(Silybum_marianum:0.022662,Centaurea_alba:0.02092)0.997000:0.026418999999999998)0.547000:0.075571,Torilis_arvensis:0.105052)0.878000:0.037612,Trifolium_repens:0.25731000000000004)1.000000:0.04059;")

weed.tree <- read.tree(text="(((Taraxacum_officinale:0.2601451047,(Silybum_marianum:0.1088645027,Centaurea_alba:0.1088645027)0.997000:0.151280602)0.547000:0.5232373813,Torilis_arvensis:0.783382486)0.878000:0.216617514,Trifolium_repens:1)1.000000:0.04059;")

#
# --- Test community-level ecoPD functions ---
#

test.phylo4com <- function() {
  # Test for error when any supplied taxa are missing from supplied tree
  BADTAXA <- weed.com$taxa
  levels(BADTAXA) <- paste("Unknown", 1:5, sep="")
  checkException(phylo4com(weed.com$plot, BADTAXA, weed.com$cover,
    phylo4d(weed.tree)))
}
 
test.pd.trad <- function() {
  DEACTIVATED("no longer calculating pd for ape-based ecophylo objects")

  # Create ecophylo objects out of weed community and phylogeny data
  ecophydata1 <- ecophylo(weed.com$plot, weed.com$taxa, weed.com$cover,
    weed.tree)
  ecophydata2 <- ecophylo(weed.com2$plot, weed.com2$taxa,
    weed.com2$cover, weed.tree)
  ecophydata.notultra <- ecophylo(weed.com$plot, weed.com$taxa,
    weed.com$cover, weed.tree.notultra)

  ans <- c(3.1523920934, 3.1523920934, 3.1523920934)
  checkEqualsNumeric(ans, ecoPD(ecophydata1, method="trad"))

  ans <- c(3.1523920934, 3.0435275907, 2.783382486)
  checkEqualsNumeric(ans, ecoPD(ecophydata2, method="trad"))

  checkException(ecoPD(ecophydata.notultra, method="trad"))
}

test.pd.poly <- function() {
  DEACTIVATED("no longer calculating pd for ape-based ecophylo objects")

  # Create ecophylo objects out of weed community and phylogeny data
  ecophydata1 <- ecophylo(weed.com$plot, weed.com$taxa, weed.com$cover,
    weed.tree)
  ecophydata2 <- ecophylo(weed.com2$plot, weed.com2$taxa,
    weed.com2$cover, weed.tree)
  ecophydata.notultra <- ecophylo(weed.com$plot, weed.com$taxa,
    weed.com$cover, weed.tree.notultra)

  ans <- c(8.26633721124392, 7.05995380972657, 7.54247469939537)
  checkEqualsNumeric(ans, ecoPD(ecophydata1, method="poly"))

  ans <- c(8.26633721124392, 5.64834888836452, 5.7174610946672)
  checkEqualsNumeric(ans, ecoPD(ecophydata2, method="poly"))

  checkException(ecoPD(ecophydata.notultra, method="poly"))
}

test.pd.yule <- function() {
  DEACTIVATED("no longer calculating pd for ape-based ecophylo objects")

  # Create ecophylo objects out of weed community and phylogeny data
  ecophydata1 <- ecophylo(weed.com$plot, weed.com$taxa, weed.com$cover,
    weed.tree)
  ecophydata2 <- ecophylo(weed.com2$plot, weed.com2$taxa,
    weed.com2$cover, weed.tree)
  ecophydata.notultra <- ecophylo(weed.com$plot, weed.com$taxa,
    weed.com$cover, weed.tree.notultra)

  ans <- c(4.72177974392577, 4.46454888807077, 4.6275784250289)
  checkEqualsNumeric(ans, ecoPD(ecophydata1, method="yule"))

  ans <- c(4.72177974392577, 3.97813049106849, 3.76931214789738)
  checkEqualsNumeric(ans, ecoPD(ecophydata2, method="yule"))

  checkException(ecoPD(ecophydata.notultra, method="yule"))
}

test.siteBySpecies <- function() {
  DEACTIVATED("no longer calculating siteBySpecies for ape-based ecophylo objects")

  # Create ecophylo objects out of weed community and phylogeny data
  ecophydata1 <- ecophylo(weed.com$plot, weed.com$taxa, weed.com$cover,
    weed.tree)

  ans <- structure(c(5, 5, 10, 25, 55, 35, 35, 15, 10, 5, 20, 20, 20,
    20, 20), .Dim = c(5L, 3L), .Dimnames = list(c("Centaurea_alba",
    "Silybum_marianum", "Taraxacum_officinale", "Torilis_arvensis",
    "Trifolium_repens"), c("A", "B", "C")))
  checkEquals(ans, siteBySpecies(ecophydata1))
}

test.richness <- function() {
  DEACTIVATED("no longer calculating richness for ape-based ecophylo objects")

  # Create ecophylo objects out of weed community and phylogeny data
  ecophydata1 <- ecophylo(weed.com$plot, weed.com$taxa, weed.com$cover,
    weed.tree)
  ecophydata2 <- ecophylo(weed.com2$plot, weed.com2$taxa,
    weed.com2$cover, weed.tree)

  ans <- structure(c(5L, 5L, 5L), .Names = c("A", "B", "C"))
  checkEquals(richness(ecophydata1), ans)
  ans <- structure(c(5L, 4L, 3L), .Names = c("A", "B", "C"))
  checkEquals(richness(ecophydata2), ans)
}
