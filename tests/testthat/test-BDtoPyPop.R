# Characterization tests for R/BDtoPyPop.R offline/pure paths: formatHead(),
# pypopHeaders(), convertAny(), and BDtoPyPop(..., save.file = FALSE).
# Uses the bundled sHLAdata BIGDAWG-style dataset.

test_that("formatHead suffixes paired locus columns with _1 and _2", {
  out <- formatHead(colnames(sHLAdata))
  expect_identical(
    out,
    c("Subject", "Status", "A_1", "A_2", "C_1", "C_2", "B_1", "B_2",
      "DRB1_1", "DRB1_2", "DQA1_1", "DQA1_2", "DQB1_1", "DQB1_2",
      "DPA1_1", "DPA1_2", "DPB1_1", "DPB1_2")
  )
  # First two columns (identifier + status) are left untouched.
  expect_identical(out[1:2], c("Subject", "Status"))
})

test_that("pypopHeaders produces PyPop-style unique locus headers", {
  out <- pypopHeaders(colnames(sHLAdata))
  expect_identical(
    out,
    c("Subject", "Status", "A_1", "A_2", "C_1", "C_2", "B_1", "B_2",
      "DRB1_1", "DRB1_2", "DQA1_1", "DQA1_2", "DQB1_1", "DQB1_2",
      "DPA1_1", "DPA1_2", "DPB1_1", "DPB1_2")
  )
  # Already-underscored headers pass through unchanged.
  underscored <- c("Subject", "Status", "A_1", "A_2")
  expect_identical(pypopHeaders(underscored), underscored)
})

test_that("convertAny replaces NA with the default placeholder", {
  out <- convertAny(sHLAdata)
  # No NA values remain.
  expect_false(any(is.na(out)))
  # The default placeholder '****' appears where NAs were.
  expect_true(any(out == "****", na.rm = TRUE))
  # Dimensions are preserved.
  expect_equal(dim(out), dim(sHLAdata))
})

test_that("convertAny supports custom change.from / change.to", {
  v <- c("x", "y", "x", "z")
  expect_identical(convertAny(v, change.from = "x", change.to = "Q"),
                   c("Q", "y", "Q", "z"))
})

test_that("BDtoPyPop returns a list of case/control frames when save.file = FALSE", {
  bd <- BDtoPyPop(sHLAdata, "BDHLA", FALSE)
  expect_type(bd, "list")
  expect_length(bd, 2L)
  # NOTE: the second list element name carries a spelling typo ('neagtive')
  # in the current implementation; encoded here as the actual behavior.
  expect_identical(names(bd), c("BDHLA.positive", "BDHLA.neagtive"))
  # The positive set holds Status == 1 rows; negative holds Status == 0 rows.
  expect_equal(dim(bd[["BDHLA.positive"]]), c(23L, 18L))
  expect_equal(dim(bd[["BDHLA.neagtive"]]), c(24L, 18L))
  # Columns were converted to PyPop headers.
  expect_identical(
    colnames(bd[["BDHLA.positive"]]),
    c("Subject", "Status", "A_1", "A_2", "C_1", "C_2", "B_1", "B_2",
      "DRB1_1", "DRB1_2", "DQA1_1", "DQA1_2", "DQB1_1", "DQB1_2",
      "DPA1_1", "DPA1_2", "DPB1_1", "DPB1_2")
  )
  # Every row in the positive set has Status 1; negative set Status 0.
  expect_true(all(bd[["BDHLA.positive"]][["Status"]] %in% c(1, "1")))
  expect_true(all(bd[["BDHLA.neagtive"]][["Status"]] %in% c(0, "0")))
})
