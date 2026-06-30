# Characterization tests for R/SearchandBuildAlignments.R against
# alignmentFull_baseline.rds (3.64.0; A, DRB1, DPB1). The HLAalignments fixture
# is wired up by helper-alignments.R.
#
# Covers: alignmentSearch, uniSearch, multiSearch, customAlign, uniAlign,
# multiAlign, validatePositions, queryPositions, variantTable.

# ---- uniSearch / multiSearch / alignmentSearch ----------------------------

test_that("uniSearch returns the residue at a single position", {
  # Prefix on (default): position number prepended to the residue.
  expect_equal(
    suppressMessages(uniSearch("prot", "DRB1", "01:01:01:01", 13)),
    "13F"
  )
  # Prefix off: just the residue.
  expect_equal(
    suppressMessages(uniSearch("prot", "DRB1", "01:01:01:01", 13, prefix = FALSE)),
    "F"
  )
})

test_that("multiSearch concatenates residues across positions", {
  expect_equal(
    suppressMessages(multiSearch("prot", "DRB1", "01:01:01:01", c(11, 13))),
    "11L~13F"
  )
})

test_that("alignmentSearch returns a position-prefixed sequence block", {
  expect_equal(
    suppressMessages(alignmentSearch("prot", "DRB1*01:01:01:01", 11:13)),
    "11L~12K~13F"
  )
})

test_that("alignmentSearch resolves a 2-field (trimmed) allele name", {
  # DRB1*15:07:01 is a two-field name; the function falls back to the
  # trimmed_allele column and returns the first matching full-length allele.
  expect_equal(
    suppressMessages(alignmentSearch("nuc", "DRB1*15:07:01", 11:13)),
    "11T~12G~13A"
  )
})

test_that("alignmentSearch errors on an invalid locus", {
  expect_error(
    alignmentSearch("prot", "ZZZ*01:01", 11),
    "is not a valid HLA allele name"
  )
})

test_that("alignmentSearch reports when no valid positions are supplied", {
  # All positions absent from the alignment -> a 'There are no valid positions.'
  # message is emitted and NULL (the message() return value) is produced.
  expect_message(
    alignmentSearch("prot", "DRB1*01:01:01:01", 999999),
    "no valid positions"
  )
  # The return value is NULL.
  expect_null(suppressMessages(alignmentSearch("prot", "DRB1*01:01:01:01", 999999)))
})

# ---- validatePositions ----------------------------------------------------

test_that("validatePositions is TRUE for present positions, FALSE otherwise", {
  expect_true(validatePositions("prot", "DRB1", c(11, 13)))
  expect_false(suppressMessages(validatePositions("prot", "DRB1", 999999)))
  expect_message(validatePositions("prot", "DRB1", 999999),
                 "not found in the DRB1 prot alignment")
})

test_that("validatePositions is FALSE for an invalid locus", {
  expect_false(suppressMessages(validatePositions("prot", "ZZZ", 11)))
})

# ---- variantTable / queryPositions ---------------------------------------

test_that("variantTable summarises variant counts and frequencies", {
  vt <- variantTable("prot", "DRB1", 13)
  expect_s3_class(vt, "data.frame")
  expect_named(vt, c("Variant", "Count", "Frequency"))
  # Lock the variant set and a couple of key counts at DRB1 position 13.
  expect_setequal(
    vt$Variant,
    c("*", ".", "C", "D", "F", "G", "H", "K", "L", "P", "Q", "R", "S", "Y")
  )
  expect_equal(vt$Count[vt$Variant == "S"], 1703L)
  expect_equal(vt$Count[vt$Variant == "H"], 606L)
  # Frequencies sum to 1 (every allele contributes exactly one variant).
  expect_equal(sum(vt$Frequency), 1, tolerance = 1e-9)
})

test_that("variantTable collapses codon columns into whole codons", {
  # For 'codon' alignments the three codon sub-columns are pasted together.
  vt <- variantTable("codon", "A", 1)
  expect_s3_class(vt, "data.frame")
  expect_true(all(nchar(vt$Variant) == 3))
  # The reference codon GGC is the most common at A position 1.
  expect_equal(vt$Count[vt$Variant == "GGC"], 6468L)
})

test_that("queryPositions returns per-position variant vectors", {
  qp <- queryPositions("prot", "DRB1", c(11, 13))
  expect_type(qp, "list")
  expect_named(qp, c("DRB1_11", "DRB1_13"))
  expect_setequal(
    qp$DRB1_13,
    c("*", ".", "C", "D", "F", "G", "H", "K", "L", "P", "Q", "R", "S", "Y")
  )
})

test_that("queryPositions returns count/frequency tables when table = TRUE", {
  qpt <- queryPositions("prot", "DRB1", 13, table = TRUE)
  expect_type(qpt, "list")
  expect_named(qpt, "DRB1_13")
  expect_s3_class(qpt[[1]], "data.frame")
  expect_equal(dim(qpt[[1]]), c(14L, 3L))
})

test_that("queryPositions returns FALSE for invalid positions", {
  expect_false(suppressMessages(queryPositions("prot", "DRB1", 999999)))
})

# ---- customAlign / uniAlign / multiAlign ---------------------------------

test_that("uniAlign builds an Allele-by-position table at one position set", {
  ua <- suppressMessages(
    uniAlign("prot", c("DRB1*01:01:01:01", "DRB1*03:01:01:01"), 11:13)
  )
  expect_s3_class(ua, "data.frame")
  expect_equal(colnames(ua), c("Allele", "11", "12", "13"))
  expect_equal(ua$Allele, c("DRB1*01:01:01:01", "DRB1*03:01:01:01"))
  expect_equal(unname(unlist(ua[1, 2:4])), c("L", "K", "F"))
  expect_equal(unname(unlist(ua[2, 2:4])), c("S", "T", "S"))
})

test_that("uniAlign returns FALSE for invalid positions", {
  expect_false(suppressMessages(
    uniAlign("prot", c("DRB1*01:01:01:01"), 999999)
  ))
})

test_that("customAlign delegates to uniAlign for a single position vector", {
  ca <- suppressMessages(
    customAlign("prot", c("DRB1*01:01:01:01", "DRB1*03:01:01:01"), c(11, 12, 13))
  )
  expect_s3_class(ca, "data.frame")
  expect_equal(colnames(ca), c("Allele", "11", "12", "13"))
  expect_equal(unname(unlist(ca[1, 2:4])), c("L", "K", "F"))
})

test_that("multiAlign builds an interleaved per-allele position table", {
  ma <- suppressMessages(
    multiAlign("prot",
               c("DRB1*01:01:01:01", "DPB1*01:01:01:01"),
               list(c(11, 12), c(35, 36)))
  )
  expect_s3_class(ma, "data.frame")
  # First column is named for the first locus; rows interleave a positions row
  # and an allele-sequence row per allele.
  expect_equal(colnames(ma)[1], "DRB1")
  expect_equal(ma[[1]],
               c("DRB1*01:01:01:01", "DPB1", "DPB1*01:01:01:01"))
  # DRB1 residues at 11,12.
  expect_equal(unname(unlist(ma[1, 2:3])), c("L", "K"))
  # DPB1 residues at 35,36.
  expect_equal(unname(unlist(ma[3, 2:3])), c("Y", "A"))
})

test_that("customAlign delegates to multiAlign for a list of position vectors", {
  ca <- suppressMessages(
    customAlign("prot",
                c("DRB1*01:01:01:01", "DPB1*01:01:01:01"),
                list(c(11, 12), c(35, 36)))
  )
  expect_s3_class(ca, "data.frame")
  expect_equal(colnames(ca)[1], "DRB1")
  expect_equal(unname(unlist(ca[3, 2:3])), c("Y", "A"))
})
