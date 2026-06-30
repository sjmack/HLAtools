# Characterization tests for R/locusValidator.R: validateLocus(), multiLocusValidation().
# Uses the bundled HLAgazetteer.

test_that("validateLocus validates loci against the HLAgazetteer", {
  # Documented @examples.
  expect_true(validateLocus("DRB1", "AA"))
  expect_true(validateLocus("V", c("cDNA", "gDNA")))
  expect_true(validateLocus(c("E", "F", "G"), "gDNA"))
})

test_that("validateLocus returns FALSE with a message for invalid loci", {
  expect_false(suppressMessages(validateLocus("ZZZ", "gDNA")))
  expect_message(validateLocus("ZZZ", "gDNA"), "is not present in")
})

test_that("validateLocus treats DRB1/3/4/5 as always valid (skipped check)", {
  # The DRB1/DRB3/DRB4/DRB5 loci are short-circuited as valid regardless of source.
  expect_true(validateLocus("DRB3", "cDNA"))
  expect_true(validateLocus("DRB5", "gDNA"))
})

test_that("multiLocusValidation returns only valid loci", {
  # Documented @example: DQB8 is invalid and removed.
  expect_equal(suppressMessages(multiLocusValidation(c("DRB1", "DPB1", "DQB8"))),
               c("DRB1", "DPB1"))
  # 'D' and 'Q' are not gDNA loci and are removed.
  expect_equal(suppressMessages(multiLocusValidation(c("A", "B", "C", "D", "Q"))),
               c("A", "B", "C"))
})

test_that("multiLocusValidation emits messages for removed loci", {
  expect_message(multiLocusValidation(c("A", "D")), "is invalid in version")
})
