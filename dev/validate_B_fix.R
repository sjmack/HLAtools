# Validate the HLA-B cDNA fix: B should now build for all sources, the two
# formerly-ragged alleles should be normalized to the reference width, and the
# reference allele's sequence must be untouched.
.scriptArgs <- commandArgs(FALSE); .scriptFile <- sub("^--file=", "", .scriptArgs[grep("^--file=", .scriptArgs)])
pkg_root <- if (length(.scriptFile)) dirname(dirname(normalizePath(.scriptFile))) else normalizePath(getwd())  # repo root = parent of dev/
suppressMessages(pkgload::load_all(pkg_root, quiet = TRUE))
VERSION <- "3.64.0"

for (src in c("cDNA","gDNA","AA")) {
  cat("==== buildAlignments('B',", src, ") ====\n")
  r <- tryCatch(suppressWarnings(buildAlignments("B", src, version = VERSION)),
                error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NULL })
  if (is.null(r)) next
  tbl <- r[[1]][[ if (src=="cDNA") 2 else 1 ]]   # cDNA: take the nucleotide table
  cat("  class:", class(tbl)[1], " dim:", nrow(tbl), "x", ncol(tbl), "\n")
  # All rows must share the same width now (rectangular)
  cat("  ref allele (row 1):", paste(as.character(tbl[1, 1:4]), collapse=" / "), "\n")
  if (src == "cDNA") {
    # Locate the two formerly-ragged alleles and show their first/last few cells
    for (al in c("B*44:568Q","B*51:197")) {
      ix <- which(tbl$allele_name == al)
      if (length(ix)) {
        row <- as.character(tbl[ix, 5:ncol(tbl)])
        cat(sprintf("  %-12s present, %d position cells; tail: %s\n",
                    al, length(row), paste(tail(row, 12), collapse="")))
      } else cat("  ", al, "NOT FOUND\n")
    }
  }
}
cat("\nDONE\n")
