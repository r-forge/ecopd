#
# --- Load and create test data ---
#

data(weeds)

weeds.tree.nu <- read.tree(text="(((Taraxacum_officinale:0.058133000000000004,(Silybum_marianum:0.022662,Centaurea_alba:0.02092)0.997000:0.026418999999999998)0.547000:0.075571,Torilis_arvensis:0.105052)0.878000:0.037612,Trifolium_repens:0.25731000000000004)1.000000:0.04059;")

pcom <- phylo4com(weeds.tree, weeds.data$cover, weeds.data$plot,
  weeds.data$taxa)

pcom.nu <- phylo4com(weeds.tree.nu, weeds.data$cover,
  weeds.data$plot, weeds.data$taxa)

#
# --- Test utility functions ---
#

test.phylo4com <- function() {
  checkTrue(is(pcom, "phylo4com"))
  checkTrue(is(pcom.nu, "phylo4com"))
}

test.abundance <- function() {
  N <- data.frame(A=c(10, 5, 5, 25, 55), B=c(15, 35, 35, 10, 5),
    C=c(20, 20, 20, 20, 20), row.names=tipLabels(pcom))
  checkIdentical(abundance(pcom), N)
  checkIdentical(abundance(pcom, "B"), N["B"])
}

test.presence <- function() {
  P <- data.frame(A=c(1, 1, 1, 1, 1), B=c(1, 1, 1, 1, 1),
    C=c(1, 1, 1, 1, 1), row.names=tipLabels(pcom))
  checkIdentical(presence(pcom), P)
  checkIdentical(presence(pcom, "B"), P["B"])
}

test.siteBySpecies <- function() {
  checkIdentical(siteBySpecies(pcom), t(as.matrix(abundance(pcom))))
  checkIdentical(siteBySpecies(pcom, transpose=TRUE),
    as.matrix(abundance(pcom)))
  checkIdentical(siteBySpecies(pcom, presence.only=TRUE),
    t(as.matrix(presence(pcom))))
}

test.richness <- function() {
  R <- structure(c(5, 5, 5), .Names = c("A", "B", "C"))
  checkIdentical(richness(pcom), R)
}

test.genera <- function() {
  G <- structure(c("Taraxacum", "Silybum", "Centaurea", "Torilis",
    "Trifolium"), .Names = c("1", "2", "3", "4", "5"))
  checkIdentical(genera(pcom), list(A=G, B=G, C=G))
  checkIdentical(genera(pcom, "B"), list(B=G))
}


#
# --- Test ecoPD indices ---
#

test.pd <- function() {
  # Note: pd would be higher if root edge is included; MC confirms
  # that we don't want the root edge included
  PD <- setNames(rep(3.1523920934, 3), c("A", "B", "C"))
  checkEqualsNumeric(pd(pcom), PD)

  PD.nu <- setNames(rep(0.603679, 3), c("A", "B", "C"))
  checkEqualsNumeric(pd(pcom.nu), PD.nu)
}

test.pae <- function() {

  # was c(1.716026353, 0.5478761669, 1.0) when root edge was included
  PAE <- setNames(c(1.716656575899907, 0.547478222537882, 1.0),
    c("A", "B", "C"))
  checkEquals(pae(pcom), PAE)

  PAE.nu <- setNames(c(1.880577959657120, 0.517350975072873, 1.0),
    c("A", "B", "C"))
  checkEquals(pae(pcom.nu), PAE.nu)
}

test.iac <- function() {
  IAC <- setNames(c(2.5, 30.0, 17.5), c("A", "B", "C"))
  checkEquals(iac(pcom), IAC)

  IAC.nu <- setNames(c(2.5, 30.0, 17.5), c("A", "B", "C"))
  checkEquals(iac(pcom.nu), IAC)
}

test.ed <- function() {
  ED <- list(
    A = setNames(c(0.488711943633333, 0.413071642633333,
      0.413071642633333, 0.8375368645, 1), tipLabels(weeds)),
    B = setNames(c(0.488711943633333, 0.413071642633333,
      0.413071642633333, 0.8375368645, 1), tipLabels(weeds)),
    C = setNames(c(0.488711943633333, 0.413071642633333,
      0.413071642633333, 0.8375368645, 1), tipLabels(weeds)))
  checkEquals(ed(pcom), ED)

  EED <- setNames(rep(0.955592201964205, 3), c("A", "B", "C"))
  checkEquals(eed(pcom), EED)
  HED <- setNames(rep(1.53796631866758, 3), c("A", "B", "C"))
  checkEquals(hed(pcom), HED)

  EED.nu <- setNames(rep(0.910003541733306, 3), c("A", "B", "C"))
  checkEquals(eed(pcom.nu), EED.nu)
  HED.nu <- setNames(rep(1.46459420051489, 3), c("A", "B", "C"))
  checkEquals(hed(pcom.nu), HED.nu)
}

test.aed <- function() {

  AED <- list(
    A = setNames(c(0.0569901020683333, 0.0678765523383333,
      0.0678765523383333, 0.0361490219733333, 0.0181818181818182),
      tipLabels(weeds)),
    B = setNames(c(0.0257789252355418, 0.0137074840755418,
      0.0137074840755418, 0.0806184329578947, 0.2), tipLabels(weeds)),
    C = setNames(c(0.0244355971816667, 0.0206535821316667,
      0.0206535821316667, 0.041876843225, 0.05), tipLabels(weeds)))
  checkEquals(aed(pcom), AED)

  EAED <- setNames(c(0.970281697458673, 0.880992381357733,
    0.984480140607184), c("A", "B", "C"))
  checkEquals(eaed(pcom), EAED)
  HAED <- setNames(c(4.46831234514660, 4.05711984871128,
    4.53369859222157), c("A", "B", "C"))
  checkEquals(haed(pcom), HAED)

  EAED.nu <- setNames(c(0.983944053257198, 0.838740167266669,
    0.968547587587551), c("A", "B", "C"))
  checkEquals(eaed(pcom.nu), EAED.nu)
  HAED.nu <- setNames(c(4.53122981874033, 3.86254121208713,
    4.46032647406888), c("A", "B", "C"))
  checkEquals(haed(pcom.nu), HAED.nu)
}

test.value <- function() {
  VALUE <- setNames(list(
    setNames(c(0.230660009755555, 0.274721498230706, 0.274721498230706,
      0.146308454598398, 0.0735885391846351), tipLabels(weeds)),
    setNames(c(0.0772258038456494, 0.0410634449172336, 0.0410634449172336,
      0.241508256572558, 0.599139049747326), tipLabels(weeds)),
    setNames(c(0.155028920627140, 0.131034348010884, 0.131034348010884,
      0.265682960648679, 0.317219422702413), tipLabels(weeds))
  ), c("A", "B", "C"))
  checkEquals(value(pcom), VALUE)
}

