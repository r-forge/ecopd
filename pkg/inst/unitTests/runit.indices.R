#
# --- Create test trees (ape::phylo objects) ---
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

weed.tree <- read.tree(text="(((Taraxacum_officinale:0.2601451047,(Silybum_marianum:0.1088645027,Centaurea_alba:0.1088645027)0.997000:0.151280602)0.547000:0.5232373813,Torilis_arvensis:0.783382486)0.878000:0.216617514,Trifolium_repens:1)1.000000:0.04059;")

weed.tree.notultra <- read.tree(text="(((Taraxacum_officinale:0.058133000000000004,(Silybum_marianum:0.022662,Centaurea_alba:0.02092)0.997000:0.026418999999999998)0.547000:0.075571,Torilis_arvensis:0.105052)0.878000:0.037612,Trifolium_repens:0.25731000000000004)1.000000:0.04059;")

comtrees <- phylo4com(weed.com$plot, weed.com$taxa, weed.com$cover,
  phylo4d(weed.tree))

comtrees.notultra <- phylo4com(weed.com$plot, weed.com$taxa,
  weed.com$cover, phylo4d(weed.tree.notultra))

#
# --- Test utility functions ---
#

test.phylobase <- function() {
  checkTrue(unique(sapply(comtrees, class))=="phylo4d")
  checkTrue(unique(sapply(comtrees.notultra, class))=="phylo4d")
}

test.abundance <- function() {
  checkEqualsNumeric(abundance(comtrees$A), c(10, 5, 5, 25, 55))
  checkEqualsNumeric(abundance(comtrees$B), c(15, 35, 35, 10, 5))
  checkEqualsNumeric(abundance(comtrees$C), c(20, 20, 20, 20, 20))
}

test.genera <- function() {
  target <- structure(c("Taraxacum", "Silybum", "Centaurea", "Torilis",
    "Trifolium"), .Names = c("1", "2", "3", "4", "5"))
  checkIdentical(genera(comtrees$A), target)
  checkIdentical(genera(comtrees$B), target)
  checkIdentical(genera(comtrees$C), target)
}

test.siteBySpecies <- function() {
  target <- structure(c(10, 5, 5, 25, 55, 15, 35, 35, 10, 5, 20, 20, 20,
    20, 20), .Dim = c(5L, 3L), .Dimnames =
    list(c("Taraxacum_officinale", "Silybum_marianum", "Centaurea_alba",
    "Torilis_arvensis", "Trifolium_repens"), c("A", "B", "C")))
  checkIdentical(siteBySpecies(comtrees), t(target))
  checkIdentical(siteBySpecies(comtrees, transpose=TRUE), target)
  # presence/absence only
  target[] <- 1
  checkIdentical(siteBySpecies(comtrees, presence=TRUE, transpose=TRUE),
    target)

}

test.richness <- function() {
  target <- structure(c(5L, 5L, 5L), .Names = c("A", "B", "C"))
  checkIdentical(richness(comtrees), target)
}


#
# --- Test ecoPD indices ---
#

test.pd <- function() {
  # Note: pd would be 3.1523920934 if root edge is included; MC confirms
  # that we don't want the root edge included
  checkEqualsNumeric(pd(comtrees$A), 3.1523920934)
  checkEqualsNumeric(pd(comtrees$B), 3.1523920934)
  checkEqualsNumeric(pd(comtrees$C), 3.1523920934)

  checkEqualsNumeric(pd(comtrees.notultra$A), 0.603679)
  checkEqualsNumeric(pd(comtrees.notultra$B), 0.603679)
  checkEqualsNumeric(pd(comtrees.notultra$C), 0.603679)
}

test.PAE <- function() {
  # was 1.716026353 when root edge was included
  checkEquals(PAE(comtrees$A), 1.716656575899907)
  # was 0.5478761669 when root edge was included
  checkEquals(PAE(comtrees$B), 0.547478222537882)
  # was 1.0 when root edge was included
  checkEquals(PAE(comtrees$C), 1.0)

  checkEquals(PAE(comtrees.notultra$A), 1.880577959657120)
  checkEquals(PAE(comtrees.notultra$B), 0.517350975072873)
  checkEquals(PAE(comtrees.notultra$C), 1.000000000000000)
}

test.IAC <- function() {
  checkEquals(IAC(comtrees$A), 2.5)
  checkEquals(IAC(comtrees$B), 30.0)
  checkEquals(IAC(comtrees$C), 17.5)
  checkEquals(IAC(comtrees.notultra$A), 2.5)
  checkEquals(IAC(comtrees.notultra$B), 30.0)
  checkEquals(IAC(comtrees.notultra$C), 17.5)
}

test.ED <- function() {

  checkEqualsNumeric(ED(comtrees$A), c(0.488711943633333,
    0.413071642633333, 0.413071642633333, 0.8375368645, 1))
  checkEqualsNumeric(ED(comtrees$B), c(0.488711943633333,
    0.413071642633333, 0.413071642633333, 0.8375368645, 1))
  checkEqualsNumeric(ED(comtrees$C), c(0.488711943633333,
    0.413071642633333, 0.413071642633333, 0.8375368645, 1))

  checkEquals(EED(comtrees$A), 0.955592201964205)
  checkEquals(EED(comtrees$B), 0.955592201964205)
  checkEquals(EED(comtrees$C), 0.955592201964205)
  checkEquals(HED(comtrees$A), 1.53796631866758)
  checkEquals(HED(comtrees$B), 1.53796631866758)
  checkEquals(HED(comtrees$C), 1.53796631866758)

  checkEquals(EED(comtrees.notultra$A), 0.910003541733306)
  checkEquals(EED(comtrees.notultra$B), 0.910003541733306)
  checkEquals(EED(comtrees.notultra$C), 0.910003541733306)
  checkEquals(HED(comtrees.notultra$A), 1.46459420051489)
  checkEquals(HED(comtrees.notultra$B), 1.46459420051489)
  checkEquals(HED(comtrees.notultra$C), 1.46459420051489)
}

test.AED <- function() {
  checkEquals(EAED(comtrees$A), 0.970281697458673)
  checkEquals(EAED(comtrees$B), 0.880992381357733)
  checkEquals(EAED(comtrees$C), 0.984480140607184)
  checkEquals(HAED(comtrees$A), 4.46831234514660)
  checkEquals(HAED(comtrees$B), 4.05711984871128)
  checkEquals(HAED(comtrees$C), 4.53369859222157)

  checkEquals(EAED(comtrees.notultra$A), 0.983944053257198)
  checkEquals(EAED(comtrees.notultra$B), 0.838740167266669)
  checkEquals(EAED(comtrees.notultra$C), 0.968547587587551)
  checkEquals(HAED(comtrees.notultra$A), 4.53122981874033)
  checkEquals(HAED(comtrees.notultra$B), 3.86254121208713)
  checkEquals(HAED(comtrees.notultra$C), 4.46032647406888)
}

test.value <- function() {
  checkEqualsNumeric(value(comtrees$A), c(0.230660009755555,
    0.274721498230706, 0.274721498230706, 0.146308454598398,
    0.0735885391846351))
  checkEqualsNumeric(value(comtrees$B), c(0.0772258038456494,
    0.0410634449172336, 0.0410634449172336, 0.241508256572558,
    0.599139049747326))
  checkEqualsNumeric(value(comtrees$C), c(0.155028920627140,
    0.131034348010884, 0.131034348010884, 0.265682960648679, 
    0.317219422702413))
}

