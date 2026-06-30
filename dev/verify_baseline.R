# Refactor safety net: rebuild the baseline loci with the CURRENT source and
# assert the result is identical() to the saved baseline. Run after every
# optimization edit. Also reports build time so we can track the speedup.
#
# Exit code is 0 on identical match, 1 on any mismatch.

.scriptArgs <- commandArgs(FALSE); .scriptFile <- sub("^--file=", "", .scriptArgs[grep("^--file=", .scriptArgs)])
pkg_root <- if (length(.scriptFile)) dirname(dirname(normalizePath(.scriptFile))) else normalizePath(getwd())  # repo root = parent of dev/
suppressMessages(pkgload::load_all(pkg_root, quiet = TRUE))

fixtures <- file.path(pkg_root, "tests/testthat/fixtures")
baseline <- readRDS(file.path(fixtures, "alignmentFull_baseline.rds"))
meta     <- readRDS(file.path(fixtures, "alignmentFull_baseline_meta.rds"))

cat("Rebuilding loci =", paste(meta$loci, collapse=", "),
    "| version =", meta$version, "\n")
t <- system.time({
  rebuilt <- alignmentFull(loci = meta$loci, alignType = meta$alignType,
                           version = meta$version)
})
cat(sprintf("Rebuild elapsed: %.1fs  (baseline was %.1fs)  -> %.2fx\n",
            t["elapsed"], meta$elapsed_sec, meta$elapsed_sec / t["elapsed"]))

# identical() is the strict gate. If it fails, use all.equal to localize.
if (identical(rebuilt, baseline)) {
  cat("\nRESULT: IDENTICAL to baseline. Behavior preserved.\n")
  quit(status = 0)
} else {
  cat("\nRESULT: *** MISMATCH vs baseline ***\n")
  ae <- all.equal(rebuilt, baseline)
  cat("all.equal() differences:\n"); print(ae)
  # Drill into which alignType/locus differs.
  for (typ in names(baseline)) {
    if (!identical(rebuilt[[typ]], baseline[[typ]])) {
      cat(sprintf("  [%s] differs\n", typ))
      if (is.list(baseline[[typ]])) {
        for (lc in names(baseline[[typ]])) {
          if (!identical(rebuilt[[typ]][[lc]], baseline[[typ]][[lc]]))
            cat(sprintf("      -> locus %s differs\n", lc))
        }
      }
    }
  }
  quit(status = 1)
}
