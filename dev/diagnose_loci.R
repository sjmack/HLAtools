# Identify which loci can / cannot be built per source, with the code currently
# in the working tree. alignmentFull(loci="all") crashes on release 3.64.0; this
# pinpoints the offending locus/loci so we can characterize the bug and choose a
# buildable locus set for the before/after timing.
.scriptArgs <- commandArgs(FALSE); .scriptFile <- sub("^--file=", "", .scriptArgs[grep("^--file=", .scriptArgs)])
pkg_root <- if (length(.scriptFile)) dirname(dirname(normalizePath(.scriptFile))) else normalizePath(getwd())  # repo root = parent of dev/
suppressMessages(pkgload::load_all(pkg_root, quiet = TRUE))
VERSION <- "3.64.0"

try_build <- function(locus, src) {
  out <- tryCatch({
    r <- suppressWarnings(buildAlignments(locus, src, version = VERSION))
    # a successful build returns a list whose [[1]] holds data frames; a
    # support-check failure returns a warning object instead
    if (inherits(r, "warning") || is.character(r)) "unsupported" else "OK"
  }, error = function(e) paste0("ERROR: ", conditionMessage(e)))
  out
}

cat("source=cDNA (nuc loci)\n")
for (lc in HLAgazetteer$nuc) cat(sprintf("  %-6s %s\n", lc, try_build(lc, "cDNA")))
cat("source=gDNA (gen loci)\n")
for (lc in HLAgazetteer$gen) cat(sprintf("  %-6s %s\n", lc, try_build(lc, "gDNA")))
cat("source=AA (prot loci)\n")
for (lc in HLAgazetteer$prot) cat(sprintf("  %-6s %s\n", lc, try_build(lc, "AA")))
cat("DONE\n")
