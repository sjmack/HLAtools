# Line-level profiling of buildAlignments() to locate the CPU hotspots.
# The baseline build showed ~90% of time is CPU (user 124s of 139s elapsed),
# so we profile the parsing/data-manipulation inside buildAlignments, not the
# network. We pin a concrete version so downloads are cached-stable and the
# (small) network portion doesn't dominate the signal.

.scriptArgs <- commandArgs(FALSE); .scriptFile <- sub("^--file=", "", .scriptArgs[grep("^--file=", .scriptArgs)])
pkg_root <- if (length(.scriptFile)) dirname(dirname(normalizePath(.scriptFile))) else normalizePath(getwd())  # repo root = parent of dev/
suppressMessages(pkgload::load_all(pkg_root, quiet = TRUE))
VERSION <- "3.64.0"

# Pre-warm: do one throwaway build so the download for A is done before we
# start the CPU profiler (we want to measure parsing, not the first fetch).
invisible(suppressWarnings(buildAlignments("A", "AA", version = VERSION)))

prof <- tempfile(fileext = ".out")
Rprof(prof, line.profiling = TRUE, interval = 0.01)
# Heaviest representative calls: a wide class I cDNA build, a class I gDNA
# build, and the very-wide DRB1 gDNA build (662 x 18527).
invisible(suppressWarnings(buildAlignments("A",    "cDNA", version = VERSION)))
invisible(suppressWarnings(buildAlignments("A",    "gDNA", version = VERSION)))
invisible(suppressWarnings(buildAlignments("DRB1", "gDNA", version = VERSION)))
Rprof(NULL)

sr <- summaryRprof(prof, lines = "show")

cat("=== TOP BY SELF TIME (lines) ===\n")
print(utils::head(sr$by.self, 25))

cat("\n=== TOP BY TOTAL TIME (lines) ===\n")
print(utils::head(sr$by.total, 25))

cat("\n=== TOTAL SAMPLED TIME ===\n")
print(sr$sampling.time)
