# Characterization tests for offline helpers in R/accessoryFunctions.R:
# numFields(), checkAlignType(), checkSource(), typeToSource(), parseAlignmentHead(),
# addCodonLine(). Network-reading helpers (posSort, validateAllele) are excluded.

test_that("numFields counts colon-delimited fields", {
  expect_equal(numFields("HLA-A*01:01"), 2L)
  expect_equal(numFields("DRB1*04:03:01"), 3L)
  expect_equal(numFields("13:02:01:01"), 4L)
  # Digit-delimited names have no colons -> a single field.
  expect_equal(numFields("A*0101"), 1L)
})

test_that("checkAlignType keeps valid values and warns about invalid ones", {
  # 'gDNA' is not a valid alignType; it is dropped with a message.
  expect_message(out <- checkAlignType(c("nuc", "prot", "gDNA")),
                 "not a valid 'alignType'")
  expect_equal(out, c("nuc", "prot"))
  # All-valid input passes through silently.
  expect_equal(checkAlignType(c("nuc", "prot")), c("nuc", "prot"))
  # No valid values -> error.
  expect_error(checkAlignType(c("foo", "bar")), "are valid 'alignType'")
})

test_that("checkSource keeps valid values and warns about invalid ones", {
  expect_message(out <- checkSource(c("AA", "cDNA", "codon")),
                 "not a valid 'source'")
  expect_equal(out, c("AA", "cDNA"))
  expect_equal(checkSource(c("AA", "cDNA", "gDNA")), c("AA", "cDNA", "gDNA"))
  expect_error(checkSource(c("foo", "bar")), "are valid 'source'")
})

test_that("typeToSource converts between alignType and source", {
  # alignType -> source (nuc/codon both map to cDNA, deduplicated).
  expect_equal(typeToSource(c("nuc", "prot", "gen"), TRUE),
               c("cDNA", "AA", "gDNA"))
  # source -> alignType (cDNA expands to both nuc and codon).
  expect_equal(typeToSource(c("AA", "cDNA", "gDNA"), FALSE),
               c("prot", "nuc", "gen", "codon"))
})

test_that("parseAlignmentHead returns header-parsing guides for local versions", {
  # Documented @example: parseAlignmentHead("3.25.0").
  expect_equal(parseAlignmentHead("3.25.0"), c(2, 22))
  # field2 >= 32 -> line 3, skip 24.
  expect_equal(parseAlignmentHead("3.32.0"), c(3, 24))
  # field2 in c(25,27:31) -> line 2, skip 22.
  expect_equal(parseAlignmentHead("3.27.0"), c(2, 22))
  # field2 in c(0:24,26) -> line 2, skip 18.
  expect_equal(parseAlignmentHead("3.10.0"), c(2, 18))
})

test_that("addCodonLine inserts AA codon lines into a cDNA alignment", {
  # Build a minimal alignment vector where the 'afterLine' (8th) element value
  # appears twice, so two AA codon lines are inserted before each occurrence.
  cDNAalign <- c("h1", "h2", "h3", "h4", "h5", "h6", "h7",
                 "CDNAMARK", "x", "y", "CDNAMARK", "z")
  out <- addCodonLine(cDNAalign, firstPos = 1, afterLine = 8, codons = 25)
  # Two new " AA codon" lines should have been added.
  added <- grep("AA codon", out, value = TRUE)
  expect_equal(length(added), 2L)
  expect_true(all(grepl("^ AA codon", added)))
  # Output is longer than input by the number of inserted lines.
  expect_equal(length(out), length(cDNAalign) + 2L)
})
