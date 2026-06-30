# Standalone test runner that does NOT require installing the package.
# Loads HLAtools via pkgload::load_all(), then runs the testthat suite
# against tests/testthat/ and prints a summary.
#
# Run with:
#   Rscript dev/run_tests.R        (the script locates the repo root itself)

# Run skip_on_cran() tests too (the HLA-B network regression test). devtools
# sets this automatically; a bare Rscript run does not, so set it here. The
# network tests still self-skip when offline via skip_if_offline().
Sys.setenv(NOT_CRAN = "true")

# Resolve the repo root from this script's own location so this runs from any
# working directory, both locally and in CI -- no hardcoded path.
.scriptArgs <- commandArgs(FALSE); .scriptFile <- sub("^--file=", "", .scriptArgs[grep("^--file=", .scriptArgs)])
pkg_root <- if (length(.scriptFile)) dirname(dirname(normalizePath(.scriptFile))) else normalizePath(getwd())  # repo root = parent of dev/

# Load the package source so all (including non-exported) functions and
# bundled data objects are available to the tests.
suppressMessages(pkgload::load_all(pkg_root, quiet = TRUE))

library(testthat)

# Run every test-*.R file in the suite. stop_on_failure = FALSE so that the
# full tally is reported even if some tests fail.
res <- test_dir(
  file.path(pkg_root, "tests", "testthat"),
  reporter = "summary",
  stop_on_failure = FALSE
)

# Aggregate the per-file results into a single pass/fail/warn tally.
df <- as.data.frame(res)
total_pass    <- sum(df$passed,   na.rm = TRUE)
total_fail    <- sum(df$failed,   na.rm = TRUE)
total_warn    <- sum(df$warning,  na.rm = TRUE)
total_skip    <- sum(df$skipped,  na.rm = TRUE)
n_files       <- length(unique(df$file))

cat("\n================ TEST SUMMARY ================\n")
cat("Test files :", n_files, "\n")
cat("PASS       :", total_pass, "\n")
cat("FAIL       :", total_fail, "\n")
cat("WARN       :", total_warn, "\n")
cat("SKIP       :", total_skip, "\n")
cat("=============================================\n")

# Non-zero exit code when anything failed, so CI can detect it.
if (total_fail > 0) {
  quit(status = 1, save = "no")
}
