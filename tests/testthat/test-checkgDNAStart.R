# Characterization test for R/checkgDNAStart.R: checkgDNAstart().
# Reads the bundled HLAatlas object; fully offline.

test_that("checkgDNAstart returns a one-row data frame of first-feature positions", {
  out <- suppressMessages(checkgDNAstart())
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1L)
  # One column per gDNA-aligned locus.
  expect_equal(ncol(out), length(HLAtools::HLAatlas$gen))
  expect_identical(colnames(out), names(HLAtools::HLAatlas$gen))
})

test_that("checkgDNAstart reports the expected non-1 start positions", {
  out <- suppressMessages(checkgDNAstart())
  # Most loci start at position 1; a handful of pseudogenes do not. These are
  # the current values captured as the golden baseline.
  expect_equal(out$A, 1)
  expect_equal(out$N, 226)
  expect_equal(out$P, 476)
  expect_equal(out$R, 449)
  expect_equal(out$S, 223)
  expect_equal(out$T, 674)
  expect_equal(out$U, 52)
  expect_equal(out$W, 498)
})
