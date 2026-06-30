# Fast per-call regression harness for iterating on buildAlignments().
#
# Records / checks the output of a focused set of buildAlignments() calls that
# together exercise every CPU hotspot line, so we can verify behavior in ~50s
# instead of the ~140s full alignmentFull baseline.
#
# Usage:
#   Rscript verify_fast.R record   # build & save per-call baselines (run on clean code)
#   Rscript verify_fast.R          # rebuild & assert identical() to saved baselines
#
# Exit 0 = all identical, 1 = any mismatch.

.scriptArgs <- commandArgs(FALSE); .scriptFile <- sub("^--file=", "", .scriptArgs[grep("^--file=", .scriptArgs)])
pkg_root <- if (length(.scriptFile)) dirname(dirname(normalizePath(.scriptFile))) else normalizePath(getwd())  # repo root = parent of dev/
suppressMessages(pkgload::load_all(pkg_root, quiet = TRUE))

VERSION  <- "3.64.0"
fixtures <- file.path(pkg_root, "tests/testthat/fixtures")
fast_rds <- file.path(fixtures, "buildAlignments_fast_baseline.rds")

# Focused call set covering all hotspots:
#   A/AA    -> strsplit explosion (405), AA dash-distribute (606-607), pad (409)
#   A/cDNA  -> full cDNA path: 405,415, NA-fill (447/448/463/482), dash (616-617)
#   A/gDNA  -> gDNA path, negative DNA_start
#   DRB1/gDNA -> very wide (18527 cols): dash loop (616-617) stress + DRB handling
#   DPB1/cDNA -> DP naming path
calls <- list(
  c("A",    "AA"),
  c("A",    "cDNA"),
  c("A",    "gDNA"),
  c("DRB1", "gDNA"),
  c("DPB1", "cDNA")
)

run_calls <- function() {
  out <- vector("list", length(calls))
  for (k in seq_along(calls)) {
    lc <- calls[[k]][1]; src <- calls[[k]][2]
    out[[k]] <- suppressWarnings(buildAlignments(lc, src, version = VERSION))
  }
  names(out) <- vapply(calls, function(c) paste(c, collapse="/"), character(1))
  out
}

mode <- commandArgs(trailingOnly = TRUE)
if (length(mode) && mode[1] == "record") {
  cat("Recording per-call baselines for:",
      paste(names(setNames(calls, NULL)), collapse=", "), "\n")
  t <- system.time(base <- run_calls())
  saveRDS(base, fast_rds)
  cat(sprintf("Recorded %d calls in %.1fs -> %s\n", length(base), t["elapsed"], fast_rds))
  quit(status = 0)
}

base <- readRDS(fast_rds)
t <- system.time(now <- run_calls())
cat(sprintf("Fast rebuild: %.1fs\n", t["elapsed"]))
ok <- TRUE
for (nm in names(base)) {
  same <- identical(now[[nm]], base[[nm]])
  cat(sprintf("  %-12s %s\n", nm, if (same) "OK (identical)" else "*** MISMATCH ***"))
  if (!same) {
    ok <- FALSE
    cat("     all.equal: "); print(all.equal(now[[nm]], base[[nm]]))
  }
}
quit(status = if (ok) 0 else 1)
