# Characterization tests for the alignment-consuming motifMatch() in
# R/motifMatch.R against alignmentFull_baseline.rds (3.64.0; A, DRB1, DPB1).
# validateMotif() (the offline helper) is already covered in test-motifMatch.R.
# The HLAalignments fixture is wired up by helper-alignments.R.

test_that("motifMatch returns full-length allele names sharing a single motif", {
  # A*-21M is the documented protein leader-peptide motif. Every A allele in
  # the fixture carries -21M, so the full A allele list is returned.
  out <- suppressMessages(motifMatch("A*-21M", "prot", full = TRUE))
  expect_type(out, "character")
  # Lock the current count of matching A alleles.
  expect_equal(length(out), 6459L)
  # Full-length names are returned (4 colon-delimited fields).
  expect_true(all(grepl("^A\\*", out)))
  expect_equal(out[1], "A*01:01:01:01")
})

test_that("motifMatch returns trimmed names when full = FALSE", {
  out <- suppressMessages(motifMatch("A*-21M", "prot", full = FALSE))
  expect_type(out, "character")
  # Same matched rows, but the two-field 'trimmed_allele' names are returned.
  expect_equal(length(out), 6459L)
  expect_equal(out[1], "A*01:01")
})

test_that("motifMatch intersects multiple variants in a motif", {
  # DRB1*11L~13F selects alleles carrying BOTH 11L and 13F (set intersection).
  out <- suppressMessages(motifMatch("DRB1*11L~13F", "prot", full = TRUE))
  expect_type(out, "character")
  expect_equal(length(out), 259L)
  expect_equal(out[1], "DRB1*01:01:01:01")
})

test_that("motifMatch returns NA with a message when no allele matches", {
  # A*1Z is a well-formed-but-absent motif (Z is allowed past validateMotif's
  # invalid-character screen? No -- 'Z' is screened, so the motif is invalid and
  # FALSE is returned). Use an absent-but-valid variant instead to hit the
  # not-found branch: pick a residue that does not occur at the position.
  # Position -21 of A is M in every allele; X never occurs there.
  out <- suppressMessages(motifMatch("A*-21X", "prot", full = TRUE))
  expect_true(is.na(out))
  expect_message(motifMatch("A*-21X", "prot"),
                 "was not found in any prot alleles")
})

test_that("motifMatch returns FALSE for an invalid locus", {
  # validateMotif() rejects the locus, so motifMatch() short-circuits to FALSE.
  expect_false(suppressMessages(motifMatch("ZZZ*1A", "prot")))
})

test_that("motifMatch errors when more than one alignType is supplied", {
  expect_error(motifMatch("A*-21M", c("prot", "nuc")),
               "single 'alignType'")
})
