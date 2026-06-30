# Compare the 40-locus (all-except-B) builds from original vs optimized code and
# assert identical() across every locus and alignment type -- the near-full-scale
# proof that the optimization preserves behavior.
.scriptArgs <- commandArgs(FALSE); .scriptFile <- sub("^--file=", "", .scriptArgs[grep("^--file=", .scriptArgs)])
pkg_root <- if (length(.scriptFile)) dirname(dirname(normalizePath(.scriptFile))) else normalizePath(getwd())  # repo root = parent of dev/
orig <- readRDS(file.path(pkg_root, "dev/most_loci_original.rds"))
opt  <- readRDS(file.path(pkg_root, "dev/most_loci_optimized.rds"))

cat("Comparing 40-locus builds (original vs optimized)...\n")
if (identical(orig, opt)) {
  cat("\nRESULT: IDENTICAL across all built loci and alignment types.\n")
  for (typ in c("prot","codon","nuc","gen"))
    cat(sprintf("  %-6s : %d loci\n", typ, length(opt[[typ]])))
  quit(status = 0)
} else {
  cat("\nRESULT: *** MISMATCH ***\n")
  for (typ in names(orig)) {
    if (!identical(orig[[typ]], opt[[typ]])) {
      cat(sprintf("  [%s] differs\n", typ))
      if (is.list(orig[[typ]])) for (lc in names(orig[[typ]]))
        if (!identical(orig[[typ]][[lc]], opt[[typ]][[lc]]))
          cat(sprintf("      -> locus %s differs\n", lc))
    }
  }
  quit(status = 1)
}
