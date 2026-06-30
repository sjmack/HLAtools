# Prove the HLA-B fix is a no-op for every other locus: the 40 non-B loci in the
# post-fix all-loci build must be identical() to the pre-fix 40-locus build.
.scriptArgs <- commandArgs(FALSE); .scriptFile <- sub("^--file=", "", .scriptArgs[grep("^--file=", .scriptArgs)])
pkg_root <- if (length(.scriptFile)) dirname(dirname(normalizePath(.scriptFile))) else normalizePath(getwd())  # repo root = parent of dev/
fixed <- readRDS(file.path(pkg_root, "dev/all_loci_optimized_fixed.rds"))  # all 41 loci, with fix
prefix <- readRDS(file.path(pkg_root, "dev/most_loci_optimized.rds"))       # 40 loci (no B), pre-fix

ok <- TRUE; checked <- 0
for (typ in c("prot","codon","nuc","gen")) {
  for (lc in names(prefix[[typ]])) {            # prefix has no B, so this skips B
    checked <- checked + 1
    if (!identical(fixed[[typ]][[lc]], prefix[[typ]][[lc]])) {
      ok <- FALSE
      cat(sprintf("  DIFF: %s / %s\n", typ, lc))
    }
  }
}
cat(sprintf("Compared %d non-B locus/type tables.\n", checked))
if (ok) cat("RESULT: all non-B tables IDENTICAL pre- vs post-fix. The fix only affects HLA-B.\n") else cat("RESULT: *** differences found ***\n")
# Also confirm B is now present in the fixed build
cat("B present in fixed build (nuc):", "B" %in% names(fixed$nuc),
    " | B cDNA dims:", paste(dim(fixed$nuc$B), collapse="x"), "\n")
quit(status = if (ok) 0 else 1)
