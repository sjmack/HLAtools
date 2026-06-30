# Characterization tests for R/AlleleTrim.R: alleleTrim(), getField(), multiAlleleTrim().
# Expected values were captured from the current implementation (golden/regression net).

test_that("alleleTrim trims epoch 3 (colon-delimited) names by field", {
  # Documented @example: alleleTrim("A*03:01:01", resolution = 2) -> "A*03:01"
  expect_equal(alleleTrim("A*03:01:01", 2), "A*03:01")
  # Default version is 3.
  expect_equal(alleleTrim("A*02:01:01:02L", 2, version = 3, append = TRUE), "A*02:01L")
  # Requesting a resolution >= the current number of fields returns the input unchanged.
  expect_equal(alleleTrim("A*03:01", 4), "A*03:01")
})

test_that("alleleTrim trims epoch 2 (digit-delimited) names", {
  # Documented @example: alleleTrim("A*030101", resolution = 2, version = 2) -> "A*0301"
  expect_equal(alleleTrim("A*030101", 2, version = 2), "A*0301")
  # Documented @example: append the expression-variant suffix when truncating.
  expect_equal(alleleTrim("HLA-A*24020102L", 3, version = 2, append = TRUE), "HLA-A*240201L")
})

test_that("alleleTrim trims epoch 1 names and appends suffixes", {
  # Documented @example: alleleTrim("A*0303N", 1, version = 1, append = TRUE) -> "A*03N"
  expect_equal(alleleTrim("A*0303N", 1, version = 1, append = TRUE), "A*03N")
})

test_that("alleleTrim out-of-range resolution returns input with a message", {
  # resolution must be 1-4; otherwise the input is returned and a message is emitted.
  expect_message(out <- alleleTrim("A*03:01:01", 5), "Resolution must range")
  expect_equal(out, "A*03:01:01")
  # version must be 1-3; otherwise the input is returned and a message is emitted.
  expect_message(out2 <- alleleTrim("A*03:01:01", 2, version = 9), "Version must range")
  expect_equal(out2, "A*03:01:01")
})

test_that("getField trims colon-delimited names by field", {
  # Documented @example: getField("HLA-A*01:01:01:01", 3) -> "HLA-A*01:01:01"
  expect_equal(getField("HLA-A*01:01:01:01", 3), "HLA-A*01:01:01")
  # Documented @example: getField("DRB1*11:01:01:12N", 2, TRUE) -> "DRB1*11:01N"
  expect_equal(getField("DRB1*11:01:01:12N", 2, TRUE), "DRB1*11:01N")
  # res == 1 with append carries the suffix only.
  expect_equal(getField("DRB1*11:01:01:12N", 1, TRUE), "DRB1*11N")
})

test_that("getField edge cases", {
  # Requesting more fields than present returns the input unchanged.
  expect_equal(getField("A*01:01", 4), "A*01:01")
  # A name with no colon (length < 2 after split) is returned unchanged.
  expect_equal(getField("A*0101", 2), "A*0101")
  # A non-logical 'append' triggers an error.
  expect_error(getField("A*01:01:01", 2, append = "x"), "is not a logical")
})

test_that("multiAlleleTrim vectorizes alleleTrim across epochs", {
  alleles3 <- c("A*02:01:01:02L", "DRB1*08:07", "DQB1*04:02:01:16Q")
  expect_equal(multiAlleleTrim(alleles3, 2),
               c("A*02:01", "DRB1*08:07", "DQB1*04:02"))
  expect_equal(multiAlleleTrim(alleles3, 2, append = TRUE),
               c("A*02:01L", "DRB1*08:07", "DQB1*04:02Q"))

  alleles1 <- c("A*01010102L", "DRB1*1613N", "HLA-Cw*0322Q")
  # version 1, resolution 1, append.
  expect_equal(multiAlleleTrim(alleles1, 1, 2, TRUE),
               c("A*01L", "DRB1*16N", "HLA-Cw*03Q"))
  # version 2, resolution 2: names already at <= 2 fields are unchanged.
  expect_equal(multiAlleleTrim(alleles1, 2, 2),
               c("A*0101", "DRB1*1613N", "HLA-Cw*0322Q"))
})

test_that("multiAlleleTrim returns an unnamed vector", {
  out <- multiAlleleTrim(c("A*02:01:01", "B*07:02:01"), 2)
  expect_null(names(out))
})
