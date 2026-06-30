#!/usr/bin/env Rscript
#
# IPD-IMGT/HLA release smoke-test.
#
# The package downloads alignment data live from the ANHIG/IMGTHLA GitHub repo,
# which IMGT updates roughly quarterly. A new release can introduce data the
# alignment parser chokes on (as HLA-B cDNA did: a trailing ragged block crashed
# buildAlignments(), and with it the default alignmentFull(loci = "all")). This
# script guards against that recurring: it builds the alignments against a given
# release and exits non-zero if anything fails, so a scheduled CI job turns red
# when a new release breaks the build.
#
# Strategy: run the real user-facing path first -- alignmentFull(loci = "all").
# If it succeeds, we are done. If it fails, run a per-locus / per-source sweep so
# the log names exactly which build(s) broke (alignmentFull aborts on the first
# error and would not tell you which locus).
#
# Usage:  Rscript dev/imgt_release_check.R [version]      (default: "Latest")
#
# Note: the set of loci tested comes from the bundled HLAgazetteer, so this
# verifies that *existing* loci still build under a new release. A release that
# adds an entirely new locus would also require rebuilding the gazetteer
# (buildGazetteer()/updateAll()); that is a separate concern from this check.

# Resolve the repo root from this script's own location (dev/<script>.R) so it
# runs from any working directory, in CI, and on any machine -- no hardcoded path.
.scriptArgs <- commandArgs(FALSE); .scriptFile <- sub("^--file=", "", .scriptArgs[grep("^--file=", .scriptArgs)])
pkg_root <- if (length(.scriptFile)) dirname(dirname(normalizePath(.scriptFile))) else normalizePath(getwd())  # repo root = parent of dev/
suppressMessages(pkgload::load_all(pkg_root, quiet = TRUE))

args    <- commandArgs(trailingOnly = TRUE)
version <- if (length(args) && nzchar(args[1])) args[1] else "Latest"

latest <- tryCatch(getLatestVersion(), error = function(e) NA_character_)
resolved <- if (identical(version, "Latest")) latest else version
cat("Latest IPD-IMGT/HLA release on ANHIG/IMGTHLA:", latest, "\n")
cat("Testing version: ", version, if (!identical(version, resolved)) paste0(" (", resolved, ")") else "", "\n\n", sep = "")

# Optional GitHub Actions step-summary output.
gh_summary <- function(...) {
  f <- Sys.getenv("GITHUB_STEP_SUMMARY")
  if (nzchar(f)) cat(..., file = f, append = TRUE)
}

## ---- 1. Happy path: the real user-facing call ----------------------------
cat("==> alignmentFull(loci = 'all', alignType = 'all', version = '", version, "')\n", sep = "")
ok <- TRUE; err <- NULL
t <- system.time(
  res <- tryCatch(
    suppressWarnings(suppressMessages(
      alignmentFull(loci = "all", alignType = "all", version = version))),
    error = function(e) { ok <<- FALSE; err <<- conditionMessage(e); NULL }
  )
)

if (ok) {
  msg <- sprintf("built prot=%d codon=%d nuc=%d gen=%d in %.0fs",
                 length(res$prot), length(res$codon),
                 length(res$nuc), length(res$gen), t["elapsed"])
  cat("\nSUCCESS:", msg, "against", resolved, "\n")
  gh_summary(sprintf("### IMGT release check: PASS\n\nVersion **%s** — %s\n", resolved, msg))
  quit(status = 0)
}

## ---- 2. Failure path: pinpoint the broken build(s) -----------------------
cat("\nalignmentFull() FAILED:", err, "\n")
cat("Running per-locus diagnostic to identify the offending build(s)...\n\n")

check_source <- function(loci, src) {
  do.call(rbind, lapply(loci, function(lc) {
    out <- tryCatch({
      r <- suppressWarnings(suppressMessages(buildAlignments(lc, src, version = version)))
      if (inherits(r, "warning") || is.character(r)) c("unsupported", "")
      else { tbl <- r[[1]][[1]]; c("OK", paste0(nrow(tbl), "x", ncol(tbl))) }
    }, error = function(e) c("ERROR", conditionMessage(e)))
    data.frame(locus = lc, source = src, status = out[1], detail = out[2],
               stringsAsFactors = FALSE)
  }))
}

all <- rbind(
  check_source(HLAgazetteer$nuc,  "cDNA"),
  check_source(HLAgazetteer$gen,  "gDNA"),
  check_source(HLAgazetteer$prot, "AA")
)

errs <- all[all$status == "ERROR", , drop = FALSE]
cat(sprintf("%d/%d builds OK, %d errored.\n\n",
            sum(all$status == "OK"), nrow(all), nrow(errs)))
cat("FAILURES:\n")
for (i in seq_len(nrow(errs)))
  cat(sprintf("  %-6s %-5s : %s\n", errs$locus[i], errs$source[i], errs$detail[i]))

gh_summary(sprintf("### IMGT release check: FAIL\n\nVersion **%s** — %d build(s) errored:\n\n",
                   resolved, nrow(errs)))
for (i in seq_len(nrow(errs)))
  gh_summary(sprintf("- `%s` / `%s`: %s\n", errs$locus[i], errs$source[i], errs$detail[i]))

quit(status = 1)
