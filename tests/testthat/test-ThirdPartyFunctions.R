# Characterization tests for R/ThirdPartyFunctions.R: countSpaces().

test_that("countSpaces counts runs of spaces", {
  # Documented @examples.
  expect_equal(countSpaces("abc def"), 1)
  expect_equal(countSpaces("abc def  ghi   jkl mno"), c(1, 2, 3, 1))
})

test_that("countSpaces handles a string with no spaces", {
  # A token with no spaces yields an empty numeric vector.
  expect_equal(countSpaces("abc"), numeric())
})
