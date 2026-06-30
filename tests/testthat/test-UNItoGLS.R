# Characterization tests for R/UNItoGLS.R: validateUniformat(), UNtoGL(),
# UNItoGLS(), multiUNItoGLS(). Uses the bundled UNIFORMAT.example dataset.

test_that("validateUniformat detects permitted/forbidden characters", {
  expect_true(validateUniformat("A*02:01,A*03:01|A*02:01,A*03:02"))
  # "^" is not a UNIFORMAT operator -> FALSE with a message.
  expect_false(suppressMessages(validateUniformat("A^B")))
  expect_message(validateUniformat("A^B"), "not permitted in UNIFORMAT")
})

test_that("UNtoGL translates UNIFORMAT to GL Strings", {
  # Documented @example 1 (default pre = TRUE adds the HLA- prefix).
  expect_identical(
    UNtoGL("A*02:01,A*03:01|A*02:01,A*03:02|A*02:02,A*03:01|A*02:02,A*03:02"),
    "HLA-A*02:01/HLA-A*02:02+HLA-A*03:01/HLA-A*03:02"
  )
  # Documented @example 2 (pre = FALSE, no prefix).
  expect_identical(
    UNtoGL("A,B|A,C|D,B|D,C|E,B|E,C", pre = FALSE),
    "A/D/E+B/C"
  )
})

test_that("UNItoGLS wraps UNtoGL and rejects invalid strings", {
  expect_identical(
    suppressMessages(UNItoGLS("A*02:01,A*03:01|A*02:01,A*03:02|A*02:02,A*03:01|A*02:02,A*03:02")),
    "HLA-A*02:01/HLA-A*02:02+HLA-A*03:01/HLA-A*03:02"
  )
  expect_false(suppressMessages(UNItoGLS("A^B")))
})

test_that("multiUNItoGLS converts a vector of UNIFORMAT strings", {
  out <- suppressMessages(multiUNItoGLS(unname(as.vector(UNIFORMAT.example[2]))[[1]]))
  expect_type(out, "character")
  expect_equal(length(out), 20L)
  # The default pre = TRUE prefixes every allele with HLA-.
  expect_identical(
    out[1],
    "HLA-blank/HLA-A*02:01+HLA-A*02:01^HLA-blank/HLA-B*44:03+HLA-B*44:03^HLA-blank/HLA-DRB1*07:01+HLA-DRB1*07:01"
  )
  expect_identical(
    out[10],
    "HLA-A*24:02+HLA-A*01:01^HLA-B*15:01+HLA-B*13:02^HLA-DRB1*07:01+HLA-DRB1*11:01"
  )
})

test_that("multiUNItoGLS converts a data frame and preserves identifiers", {
  out <- suppressMessages(multiUNItoGLS(UNIFORMAT.example[1:2]))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 20L)
  expect_equal(out$sample.id[1], "hhrv_id195")
  # The genotype column is now GL String formatted (HLA- prefixes present).
  expect_true(all(grepl("HLA-", out$genotype)))
})
