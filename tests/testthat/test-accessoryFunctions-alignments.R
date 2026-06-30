# Characterization tests for the alignment-CONSUMING helpers in
# R/accessoryFunctions.R: posSort() and validateAllele(). The offline helpers
# (numFields, checkAlignType, ...) are covered in test-accessoryFunctions.R.
# The HLAalignments fixture is wired up by helper-alignments.R.

# ---- posSort --------------------------------------------------------------

test_that("posSort orders simple positions by alignment column order", {
  # Positions that exist in the DRB1 prot alignment are returned in sequence
  # (i.e. alignment column) order; input order is ignored. NOTE: posSort returns
  # the values with their INPUT type, so numeric input yields a numeric result.
  expect_equal(posSort(c(13, 11, 12), "prot", "DRB1"), c(11, 12, 13))
  # Character input yields a character result.
  expect_equal(posSort(c("13", "11", "12"), "prot", "DRB1"), c("11", "12", "13"))
})

test_that("posSort sorts text-formatted indel positions numerically", {
  # Indel columns such as 607.1/607.2/607.3 exist in the DRB1 nuc alignment and
  # must sort in sequence order, not lexicographic order.
  expect_equal(
    posSort(c("607.3", "607.1", "607.2"), "nuc", "DRB1"),
    c("607.1", "607.2", "607.3")
  )
})

test_that("posSort drops positions absent from the alignment", {
  # 999999 is not a column for DRB1 prot, so it is dropped from the result.
  # (Numeric input -> numeric result, see note above.)
  out <- posSort(c(11, 999999, 13), "prot", "DRB1")
  expect_equal(out, c(11, 13))
})

test_that("posSort errors for a locus absent from the alignment type", {
  expect_error(posSort(c(1, 2), "prot", "ZZZ"),
               "is not included among the prot alignments")
})

# ---- validateAllele -------------------------------------------------------

test_that("validateAllele is TRUE for a full-length allele in the fixture", {
  expect_true(validateAllele("A*01:01:01:01"))
})

test_that("validateAllele is TRUE for a two-field (trimmed) allele", {
  # Matched via the trimmed_allele column.
  expect_true(suppressMessages(validateAllele("A*01:01")))
})

test_that("validateAllele is FALSE for an absent allele with a message", {
  expect_false(suppressMessages(validateAllele("A*99:99:99:99")))
  expect_message(validateAllele("A*99:99:99:99"),
                 "is not found in version 3.64.0 alignments")
})

test_that("validateAllele rejects malformed names", {
  # No asterisk.
  expect_false(suppressMessages(validateAllele("A01:01")))
  expect_message(validateAllele("A01:01"), "No asterisk")
  # No colon.
  expect_false(suppressMessages(validateAllele("A*0101")))
  expect_message(validateAllele("A*0101"), "No colon")
})

test_that("validateAllele errors on an invalid locus (current behaviour)", {
  # CHARACTERIZATION of a current quirk: for a locus that fails validateLocus,
  # validateAllele() calls validateLocus(allele, "prot"). validateLocus() runs
  # checkSource("prot"), but "prot" is an alignType, not a source value, so
  # checkSource() raises an error rather than validateAllele() returning FALSE.
  # This locks the present (buggy) behaviour so a future fix is detected.
  expect_error(
    suppressMessages(validateAllele("ZZZ*01:01")),
    "are valid 'source' values"
  )
})
