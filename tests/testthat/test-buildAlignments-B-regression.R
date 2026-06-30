# Regression test for the HLA-B cDNA build failure (IPD-IMGT/HLA 3.64.0).
#
# A few IMGT alleles (e.g. B*44:568Q, B*51:197, B*13:123Q, B*51:394Q) carry
# sequence in a trailing ("appendix") alignment block that the reference allele
# does not span. The terminal-block handling in buildAlignments() padded that
# block to the LAST appendix sequence rather than the LONGEST one, producing
# ragged rows; do.call(rbind, ...) then recycled them and the build later died
# with "number of items to replace is not a multiple of replacement length",
# which crashed the default alignmentFull(loci = "all").
#
# buildAlignments() (v1.9.0) now extends the terminal block to the longest
# appendix sequence, so the build succeeds AND the extra appendix positions are
# preserved for the alleles that carry them (rather than being truncated).
#
# These tests require network access to the ANHIG/IMGTHLA GitHub repo, so they
# are skipped on CRAN and when offline. The version is pinned for reproducibility
# (IMGT release branches are immutable), which makes the exact dimensions stable.

test_that("HLA-B cDNA builds despite ragged trailing-block alleles", {
  skip_on_cran()
  skip_if_offline()

  b <- suppressWarnings(suppressMessages(
    buildAlignments("B", "cDNA", version = "3.64.0")))

  # The cDNA build returns, under the locus, list(codon, cDNA, Version).
  cdna <- b[["B"]][["cDNA"]]
  expect_s3_class(cdna, "data.frame")

  # Rectangular table: 4 metadata columns + the reference-relative position
  # columns, including the preserved appendix insertion positions (...1089.NN).
  # The bug produced a width mismatch and aborted before this point.
  expect_equal(ncol(cdna), 1742L)   # 4 metadata + 1738 position columns (v3.64.0)
  expect_equal(nrow(cdna), 11110L)

  # The formerly-ragged alleles are present and conform to the table width.
  odd <- c("B*44:568Q", "B*51:197", "B*13:123Q", "B*51:394Q")
  expect_true(all(odd %in% cdna$allele_name))
  for (a in odd) expect_equal(sum(cdna$allele_name == a), 1L)

  # The reference allele (row 1) is unchanged.
  expect_equal(cdna$allele_name[1], "B*07:02:01:01")
})

test_that("HLA-B cDNA preserves appendix sequence for the trailing-block alleles", {
  skip_on_cran()
  skip_if_offline()

  cdna <- suppressWarnings(suppressMessages(
    buildAlignments("B", "cDNA", version = "3.64.0")))[["B"]][["cDNA"]]

  # Appendix insertion positions follow the final reference position (1089) and
  # are labeled 1089.NN. They must carry real (non-".") bases for the appendix
  # alleles -- the data the previous (truncating) behavior would have dropped.
  appendix_cols <- grep("^1089\\.", colnames(cdna), value = TRUE)
  expect_gt(length(appendix_cols), 24L)   # more than the truncated build retained

  carriers <- cdna$allele_name[rowSums(cdna[appendix_cols] != ".") > 0]
  expect_true(all(c("B*44:568Q", "B*51:197", "B*13:123Q", "B*51:394Q") %in% carriers))
})
