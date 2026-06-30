# Characterization tests for the OFFLINE part of R/motifMatch.R: validateMotif().
# motifMatch() itself reads the un-bundled HLAalignments object and is excluded.
# validateMotif() only consults the HLAgazetteer (via validateLocus), so it is offline.

test_that("validateMotif accepts well-formed protein and genomic motifs", {
  # Documented @examples.
  expect_true(suppressMessages(validateMotif("A*-21M~2P", "prot")))
  expect_true(suppressMessages(validateMotif("A*196G~301A~3046T", "gen")))
})

test_that("validateMotif rejects invalid loci", {
  expect_false(suppressMessages(validateMotif("ZZZ*1A", "prot")))
})

test_that("validateMotif rejects motifs containing invalid characters", {
  # 'B' is an invalid amino-acid character for a protein motif body.
  expect_false(suppressMessages(validateMotif("A*1B", "prot")))
  expect_message(validateMotif("A*1B", "prot"), "contains invalid characters")
})
