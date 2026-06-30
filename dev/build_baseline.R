# Build the characterization baseline for alignmentFull().
#
# This runs the CURRENT (un-optimized) alignmentFull() on a representative
# subset of loci and saves the exact result. The saved object serves two
# purposes:
#   1. Test fixture -- provides an HLAalignments object so the (otherwise
#      network-dependent) alignment-consuming functions can be tested offline.
#   2. Refactor safety net -- after we optimize alignmentFull/buildAlignments,
#      the new output must be identical() to this baseline.
#
# Version is PINNED to a concrete release (not "Latest") so the baseline is
# reproducible: IPD-IMGT/HLA release branches are immutable upstream.

.scriptArgs <- commandArgs(FALSE); .scriptFile <- sub("^--file=", "", .scriptArgs[grep("^--file=", .scriptArgs)])
pkg_root <- if (length(.scriptFile)) dirname(dirname(normalizePath(.scriptFile))) else normalizePath(getwd())  # repo root = parent of dev/
suppressMessages(pkgload::load_all(pkg_root, quiet = TRUE))

VERSION <- "3.64.0"
# Representative loci:
#   A     - class I (prot/nuc/gen), exercises the redundant cDNA download path
#   DRB1  - class II, the special DR handling in buildAlignments
#   DPB1  - class II, DP naming conventions
LOCI <- c("A", "DRB1", "DPB1")

fixtures <- file.path(pkg_root, "tests/testthat/fixtures")
dir.create(fixtures, recursive = TRUE, showWarnings = FALSE)

cat("Building baseline: loci =", paste(LOCI, collapse=", "),
    "| alignType = all | version =", VERSION, "\n")

t <- system.time({
  baseline <- alignmentFull(loci = LOCI, alignType = "all", version = VERSION)
})
cat("\n--- BUILD TIME (current code) ---\n")
print(t)

# Persist the baseline object and the metadata describing how it was built.
saveRDS(baseline, file.path(fixtures, "alignmentFull_baseline.rds"))
meta <- list(version = VERSION, loci = LOCI, alignType = "all",
             elapsed_sec = unname(t["elapsed"]),
             gazetteer_version = HLAgazetteer$version)
saveRDS(meta, file.path(fixtures, "alignmentFull_baseline_meta.rds"))

cat("\n--- BASELINE STRUCTURE ---\n")
cat("top-level names:", paste(names(baseline), collapse=", "), "\n")
for (typ in c("prot","codon","nuc","gen")) {
  loci_present <- names(baseline[[typ]])
  cat(sprintf("  %-6s: loci = %s\n", typ, paste(loci_present, collapse=", ")))
  for (lc in loci_present) {
    d <- baseline[[typ]][[lc]]
    cat(sprintf("      %-6s dim = %s x %s\n", lc,
                ifelse(is.null(nrow(d)), NA, nrow(d)),
                ifelse(is.null(ncol(d)), NA, ncol(d))))
  }
}
cat("version field:", baseline$version, "\n")
cat("\nSaved baseline + meta to", fixtures, "\n")
