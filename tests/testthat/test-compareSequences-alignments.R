# Characterization tests for R/compareSequences.R against the saved
# alignmentFull_baseline.rds fixture (IPD-IMGT/HLA 3.64.0; loci A, DRB1, DPB1).
# The HLAalignments fixture is wired up by helper-alignments.R.
#
# These tests lock the CURRENT behaviour of compareSequences() so a future
# refactor that changes the comparison output is caught.

test_that("compareSequences returns a 2-row diff data frame for two A alleles", {
  out <- compareSequences("prot", c("A*01:01:01:01", "A*02:01:01:01"))
  # The result is a two-row data frame: one row per allele, plus an
  # 'allele_name' column followed by only the positions that differ.
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2L)
  # 'allele_name' is the first column; the remaining columns are differing
  # protein positions.
  expect_equal(colnames(out)[1], "allele_name")
  # Lock the exact set of differing positions and the allele identities.
  expect_equal(out$allele_name, c("A*01:01:01:01", "A*02:01:01:01"))
  expect_equal(
    colnames(out),
    c("allele_name", "-15", "44", "62", "66", "67", "74", "76", "77", "90",
      "95", "97", "105", "107", "114", "116", "127", "142", "145", "150",
      "152", "156", "158", "163", "166", "167", "184", "193", "194", "207",
      "253", "276", "294", "321")
  )
  # Spot-check a couple of known differing residues.
  expect_equal(out[["-15"]], c("L", "V"))
  expect_equal(out[["44"]], c("K", "R"))
})

test_that("compareSequences detects identical alleles", {
  # When both alleles are identical, a plain message string is returned.
  expect_identical(
    compareSequences("prot", c("A*01:01:01:01", "A*01:01:01:01")),
    "A*01:01:01:01 and A*01:01:01:01 are identical."
  )
})

test_that("compareSequences reports an allele not present in the alignment", {
  # validateAllele() fails for an unknown allele, so a descriptive string
  # (not a data frame) is returned.
  out <- suppressMessages(
    compareSequences("prot", c("A*99:99:99:99", "A*01:01:01:01"))
  )
  expect_identical(out, "A*99:99:99:99 is not found in the prot alignment for A.")
})

test_that("compareSequences compares alleles across two different loci", {
  # Cross-locus comparison is supported; only positions shared by both loci
  # (in the protein alignment) are compared.
  out <- compareSequences("prot", c("A*01:01:01:01", "DRB1*01:01:01:01"))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 2L)
  expect_equal(out$allele_name, c("A*01:01:01:01", "DRB1*01:01:01:01"))
  # Lock the number of differing shared positions (current behaviour: 244
  # columns including the leading allele_name column).
  expect_equal(ncol(out), 244L)
})

test_that("compareSequences errors when not given exactly two alleles", {
  expect_error(compareSequences("prot", "A*01:01:01:01"),
               "exactly two alleles")
  expect_error(
    compareSequences("prot", c("A*01:01:01:01", "A*02:01:01:01", "A*03:01:01:01")),
    "exactly two alleles"
  )
})

test_that("compareSequences errors when more than one alignType is given", {
  expect_error(
    compareSequences(c("prot", "nuc"), c("A*01:01:01:01", "A*02:01:01:01")),
    "only one 'alignType'"
  )
})
