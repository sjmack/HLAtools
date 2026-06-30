# helper-alignments.R
#
# Alignment-consuming functions (compareSequences, motifMatch, alignmentSearch,
# uniSearch, multiSearch, customAlign, uniAlign, multiAlign, validatePositions,
# queryPositions, variantTable, posSort, validateAllele) reference a free
# variable `HLAalignments`. That object is normally built over the network via
# alignmentFull() and assigned to the global environment. It is NOT a package
# data object, so package functions resolve `HLAalignments` lexically through
# their enclosing environments: package namespace -> imports -> base -> the
# global environment. Because it lives in none of the package environments, the
# lookup falls through to the global environment.
#
# To test these functions offline we load the saved fixture (the recorded output
# of alignmentFull(loci=c("A","DRB1","DPB1"), alignType="all", version="3.64.0"))
# and make it visible to the consumers.

# Resolve the fixture path. Under a normal testthat run, test_path() locates the
# fixtures directory regardless of the working directory. When this helper is
# sourced outside that context (e.g. a manual smoke runner), test_path() errors,
# so fall back to searching the conventional fixtures locations.
.alignment_fixture_path <- tryCatch(
  testthat::test_path("fixtures", "alignmentFull_baseline.rds"),
  error = function(e) {
    candidates <- c(
      file.path("fixtures", "alignmentFull_baseline.rds"),
      file.path("tests", "testthat", "fixtures", "alignmentFull_baseline.rds")
    )
    hit <- candidates[file.exists(candidates)]
    if (length(hit) == 0) stop("Cannot locate alignmentFull_baseline.rds fixture.")
    hit[[1]]
  }
)

# Load the fixture once and stash it for the tests.
HLAalignments_fixture <- readRDS(.alignment_fixture_path)

# Make the fixture visible to the consumer functions.
#
# Primary mechanism: assign into the GLOBAL environment, where the free-variable
# lookup ultimately resolves.
assign("HLAalignments", HLAalignments_fixture, envir = globalenv())

# Belt-and-suspenders: also try to inject the object into the package namespace
# (and its imports environment) so the lookup resolves there even if some run
# context shadows the global assignment. Wrapped in tryCatch because the
# namespace may be locked; failure here is harmless given the global assignment.
tryCatch({
  ns <- asNamespace("HLAtools")
  # Unlock the binding if it already exists, then assign.
  if (exists("HLAalignments", envir = ns, inherits = FALSE)) {
    if (bindingIsLocked("HLAalignments", ns)) unlockBinding("HLAalignments", ns)
  }
  assign("HLAalignments", HLAalignments_fixture, envir = ns)
}, error = function(e) invisible(NULL))
