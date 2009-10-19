# doRUnit.R --- Run RUnit tests
# slightly modified from fCalendar package
#------------------------------------------------------------------------

if(require("RUnit", quietly = TRUE)) {

  # Setup
  wd <- getwd()
  pkg <- sub("\\.Rcheck$", '', basename(dirname(wd)))

  # temporary fix to accommodate custom library location:
  .libPaths("~/R-dev")

  library(package=pkg, character.only=TRUE)

  if(Sys.getenv("RCMDCHECK") == "FALSE") {
    # Path to unit tests for standalone running under Makefile (not R CMD check)
    # PKG/tests/../inst/unitTests
    path <- file.path(getwd(), "..", "inst", "unitTests")
  } else {
    # Path to unit tests for R CMD check
    # PKG.Rcheck/tests/../PKG/unitTests
    path <- system.file(package=pkg, "unitTests")
  }

  stopifnot(file.exists(path), file.info(path.expand(path))$isdir)
  source(file.path(path, "runTests.R"), echo = TRUE)

}
